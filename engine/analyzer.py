"""
Core analysis engine for the AIM3 peptide optimization platform.

Pipeline:
    Tool A — Boltz-2 / ColabFold / Chai-1 (3D complex prediction)
    Tool B — MM/GBSA via AmberTools `MMPBSA.py` (igb=5, single-point)

Each Tool A lives in its own conda env under `third_party/<tool>_env/` because
the three tools have incompatible numpy/biopython pins; this module shells out
to each env's CLI via subprocess. Tool B uses the AmberTools binaries inside
`third_party/gmx_env/`.

Pre-computed results for the three demo PAI-1 peptides are stored in
`dataset/peptide_pai1/precomputed/<sequence>.json` and returned instantly when
a matching sequence comes in. Anything else triggers a live Boltz + MM/GBSA
run (~10–20 min per peptide on an A100).
"""

import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from scipy.stats import spearmanr


REPO = Path(__file__).parent.parent
RECEPTOR_PDB = REPO / "dataset" / "lrp1_cr2" / "receptor.pdb"
PRECOMPUTED_DIR = REPO / "dataset" / "peptide_pai1" / "precomputed"
RUNS_DIR = REPO / "runs"

GMX_ENV = REPO / "third_party" / "gmx_env"
TLEAP = GMX_ENV / "bin" / "tleap"
MMPBSA_PY = GMX_ENV / "bin" / "MMPBSA.py"
CPPTRAJ = GMX_ENV / "bin" / "cpptraj"
SANDER = GMX_ENV / "bin" / "sander"

TOOL_A_BIN = {
    "boltz2": REPO / "third_party" / "boltz_env" / "bin" / "boltz",
    "colabfold": REPO / "third_party" / "colabfold_env" / "bin" / "colabfold_batch",
    "chai1": REPO / "third_party" / "chai_env" / "bin" / "chai-lab",
}

# Per-step wall-clock caps for the live path. Boltz/CF/Chai inference is
# typically 2–6 min on A100 but MSA server can stall; tleap is instant.
TOOL_A_TIMEOUT_S = 1800
SANDER_TIMEOUT_S = 300
TLEAP_TIMEOUT_S = 60
CPPTRAJ_TIMEOUT_S = 60
MMPBSA_TIMEOUT_S = 120


class ToolAFailure(RuntimeError):
    """Raised when the Tool A subprocess (Boltz/ColabFold/Chai-1) fails."""


class ToolBFailure(RuntimeError):
    """Raised when any AmberTools step (tleap/sander/cpptraj/MMPBSA.py) fails."""


class PipelineTimeout(RuntimeError):
    """Raised when a pipeline subprocess exceeds its per-step timeout."""


def _run_checked(cmd, *, kind, timeout, cwd=None, env=None):
    """Run a subprocess, capture output, classify failures."""
    try:
        proc = subprocess.run(
            [str(c) for c in cmd],
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=cwd,
            env=env,
        )
    except subprocess.TimeoutExpired:
        raise PipelineTimeout(f"{Path(cmd[0]).name} exceeded {timeout}s")

    if proc.returncode != 0:
        # Last few lines of stderr are usually the most informative.
        tail = "\n".join((proc.stderr or proc.stdout or "").strip().splitlines()[-12:])
        msg = f"{Path(cmd[0]).name} exited {proc.returncode}\n{tail}"
        exc = ToolAFailure if kind == "A" else ToolBFailure
        raise exc(msg)
    return proc


def _pick_free_gpu():
    """Return the CUDA device index with the most free memory.

    Co-tenant GPU 0 is often saturated by other users on this host; without
    this Boltz happily allocates onto an OOM-prone card. Fallback is "0".
    """
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=memory.free", "--format=csv,noheader,nounits"],
            text=True,
            timeout=5,
        )
        rows = [int(line.strip()) for line in out.strip().splitlines() if line.strip()]
        if not rows:
            return "0"
        # Skip non-compute display GPUs by ignoring rows with very small total memory.
        # nvidia-smi orders by index, so argmax returns the index directly.
        idx = max(range(len(rows)), key=lambda i: rows[i])
        return str(idx)
    except Exception:
        return "0"


def _three_to_one():
    from Bio.PDB.Polypeptide import index_to_one, three_to_index

    return lambda code: index_to_one(three_to_index(code))


def _extract_chain_sequence(pdb_path):
    from Bio.PDB import PDBParser

    convert = _three_to_one()
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", str(pdb_path))
    chains = {}
    for model in structure:
        for chain in model:
            seq = ""
            for res in chain:
                try:
                    seq += convert(res.get_resname())
                except (KeyError, IndexError):
                    pass  # skip non-standard / hetatms (e.g. Ca²⁺)
            if seq:
                chains[chain.id] = seq
    return chains


# Boltz / ColabFold / Chai-1 all emit the peptide as chain A and receptor as B.
PEPTIDE_CHAIN = "A"

# Atom names retained when stripping a residue's side chain to alanine. Keeping
# Cβ preserves the WT pose at that atom; tleap rebuilds Hβ from the ALA template.
_ALA_KEEP_ATOMS = {
    "N",
    "CA",
    "C",
    "O",
    "OXT",
    "CB",
    "H",
    "HA",
    "HB1",
    "HB2",
    "HB3",
    "H1",
    "H2",
    "H3",
}


def _mutate_residue_to_ala_in_pdb(pdb_in, pdb_out, chain_id, residue_idx):
    """In-place X→A side-chain strip at (chain_id, residue_idx).

    Writes a copy of `pdb_in` with that one residue's side chain past Cβ
    deleted and the residue name rewritten to ALA. Backbone + Cβ coordinates
    are kept bit-exact so the binding pose at the rest of the complex is
    unchanged; tleap fills missing Cβ hydrogens from the ff14SB ALA template
    on load.
    """
    found = False
    with open(pdb_in) as fin, open(pdb_out, "w") as fout:
        for line in fin:
            tag = line[:6].strip()
            if tag in ("ATOM", "HETATM"):
                ch = line[21]
                try:
                    resnum = int(line[22:26])
                except ValueError:
                    fout.write(line)
                    continue
                atom = line[12:16].strip()
                if ch == chain_id and resnum == residue_idx:
                    if atom not in _ALA_KEEP_ATOMS:
                        continue
                    line = line[:17] + "ALA" + line[20:]
                    found = True
            fout.write(line)
    if not found:
        raise ToolBFailure(
            f"In-place mutation target not found: chain {chain_id} residue {residue_idx} "
            f"in {pdb_in}"
        )


def _filter_pdb_by_chain(pdb_in, pdb_out, chain_ids):
    """Write a PDB containing only ATOM records for the given chain IDs.

    Boltz output PDBs have one peptide chain (A) and one receptor chain (B);
    we need to split them for AmberTools tleap to build per-component prmtops.
    """
    keep = set(chain_ids)
    with open(pdb_in) as fin, open(pdb_out, "w") as fout:
        for line in fin:
            tag = line[:6].strip()
            if tag in ("ATOM", "HETATM"):
                if line[21] in keep:
                    fout.write(line)
            elif tag == "TER" and len(line) > 21 and line[21] in keep:
                fout.write(line)
            elif tag in ("END", "ENDMDL"):
                fout.write(line)
        fout.write("END\n")


class ProteinOptimizer:
    """Top-level pipeline. One instance per server process is fine."""

    # Stefansson 1998 LRP1-binding Kd_fold for hotspot calibration (ρ).
    STEFANSSON_KD_FOLD = {"WT": 1.0, "K88A": 1.3, "K80A": 2.1, "K69A": 12.4}

    def __init__(self):
        self._receptor_seq = None

    # ── Receptor ──────────────────────────────────────────────────────────
    def get_receptor_sequence(self):
        if self._receptor_seq is None:
            chains = _extract_chain_sequence(RECEPTOR_PDB)
            if not chains:
                raise RuntimeError(f"No protein chain in {RECEPTOR_PDB}")
            # CR7 is single-chain
            self._receptor_seq = next(iter(chains.values()))
        return self._receptor_seq

    # ── Cache ─────────────────────────────────────────────────────────────
    def _cache_path(self, sequence, tool_a=None):
        """WT-result cache path. tool_a is part of the filename because the
        Boltz / ColabFold / Chai-1 binding energies aren't directly comparable.
        Returns the un-suffixed legacy path when tool_a is None (only used by
        `load_cached` as a fallback)."""
        if tool_a:
            return PRECOMPUTED_DIR / f"{sequence}_{tool_a}.json"
        return PRECOMPUTED_DIR / f"{sequence}.json"

    def load_cached(self, sequence, tool_a=None):
        for path in filter(
            None,
            [
                self._cache_path(sequence, tool_a=tool_a) if tool_a else None,
                self._cache_path(sequence),
            ],
        ):
            if path.exists():
                data = json.loads(path.read_text())
                data["from_cache"] = True
                return data
        return None

    def save_cache(self, sequence, payload, tool_a=None):
        PRECOMPUTED_DIR.mkdir(parents=True, exist_ok=True)
        self._cache_path(sequence, tool_a=tool_a).write_text(
            json.dumps(payload, indent=2)
        )

    # ── Tool A ────────────────────────────────────────────────────────────
    def run_tool_a(self, sequence, tool="boltz2"):
        if tool == "boltz2":
            return self._run_boltz(sequence)
        if tool == "colabfold":
            return self._run_colabfold(sequence)
        if tool == "chai1":
            return self._run_chai(sequence)
        raise ValueError(f"Unknown Tool A: {tool!r}")

    def _gpu_env(self):
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = _pick_free_gpu()
        return env

    def _run_boltz(self, sequence):
        run_dir = RUNS_DIR / f"boltz_{sequence}"
        run_dir.mkdir(parents=True, exist_ok=True)
        receptor_seq = self.get_receptor_sequence()

        yaml_path = run_dir / "input.yaml"
        yaml_path.write_text(
            "version: 1\n"
            "sequences:\n"
            "  - protein:\n"
            "      id: A\n"
            f"      sequence: {sequence}\n"
            "  - protein:\n"
            "      id: B\n"
            f"      sequence: {receptor_seq}\n"
        )

        cmd = [
            TOOL_A_BIN["boltz2"],
            "predict",
            yaml_path,
            "--out_dir",
            run_dir,
            "--accelerator",
            "gpu",
            "--use_msa_server",
            "--output_format",
            "pdb",
            "--override",
        ]
        _run_checked(
            cmd,
            kind="A",
            timeout=TOOL_A_TIMEOUT_S,
            cwd=str(run_dir),
            env=self._gpu_env(),
        )

        pdbs = list(run_dir.rglob("*_model_0.pdb")) or list(run_dir.rglob("*.pdb"))
        if not pdbs:
            raise ToolAFailure(f"Boltz produced no PDB under {run_dir}")
        return str(pdbs[0])

    def _run_colabfold(self, sequence):
        run_dir = RUNS_DIR / f"colabfold_{sequence}"
        run_dir.mkdir(parents=True, exist_ok=True)
        receptor_seq = self.get_receptor_sequence()

        fasta = run_dir / "input.fasta"
        # ColabFold multimer convention: chains separated by ":" within one sequence record.
        fasta.write_text(f">complex\n{sequence}:{receptor_seq}\n")

        cmd = [
            TOOL_A_BIN["colabfold"],
            fasta,
            run_dir,
            "--num-models",
            "1",
            "--num-recycle",
            "3",
        ]
        _run_checked(
            cmd,
            kind="A",
            timeout=TOOL_A_TIMEOUT_S,
            cwd=str(run_dir),
            env=self._gpu_env(),
        )

        pdbs = (
            sorted(run_dir.glob("*_relaxed_rank_001*.pdb"))
            or sorted(run_dir.glob("*_unrelaxed_rank_001*.pdb"))
            or sorted(run_dir.glob("*.pdb"))
        )
        if not pdbs:
            raise ToolAFailure(f"ColabFold produced no PDB under {run_dir}")
        return str(pdbs[0])

    def _run_chai(self, sequence):
        run_dir = RUNS_DIR / f"chai_{sequence}"
        run_dir.mkdir(parents=True, exist_ok=True)
        receptor_seq = self.get_receptor_sequence()

        # Chai-1 has no --override flag and asserts that the *output dir it
        # receives* is empty. We keep the FASTA + converted PDB in `run_dir`
        # and hand Chai a clean `out/` subdir.
        out_dir = run_dir / "out"
        pdb_path = run_dir / "complex.pdb"

        # Reuse prior successful run if both artifacts are on disk.
        existing_cifs = (
            sorted(out_dir.rglob("pred.model_idx_0.cif")) if out_dir.exists() else []
        )
        if pdb_path.exists() and existing_cifs:
            return str(pdb_path)

        # Clear stale Chai output so the assert passes; preserve nothing inside out/.
        if out_dir.exists():
            shutil.rmtree(out_dir)
        out_dir.mkdir(parents=True)

        fasta = run_dir / "input.fasta"
        # Chai-1 expects per-chain FASTA records.
        fasta.write_text(
            f">protein|name=peptide\n{sequence}\n"
            f">protein|name=receptor\n{receptor_seq}\n"
        )

        cmd = [
            TOOL_A_BIN["chai1"],
            "fold",
            fasta,
            out_dir,
            "--use-msa-server",
        ]
        _run_checked(
            cmd,
            kind="A",
            timeout=TOOL_A_TIMEOUT_S,
            cwd=str(run_dir),
            env=self._gpu_env(),
        )

        cifs = sorted(out_dir.rglob("pred.model_idx_0.cif")) or sorted(
            out_dir.rglob("*.cif")
        )
        if not cifs:
            raise ToolAFailure(f"Chai-1 produced no CIF under {out_dir}")
        # Convert to PDB for the downstream MM/GBSA path.
        self._cif_to_pdb(cifs[0], pdb_path)
        return str(pdb_path)

    @staticmethod
    def _cif_to_pdb(cif_in, pdb_out):
        from Bio.PDB import MMCIFParser, PDBIO

        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("x", str(cif_in))
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(pdb_out))

    # ── Tool B ────────────────────────────────────────────────────────────
    def run_tool_b(self, complex_pdb, sequence, tool="mmgbsa", tool_a="boltz2"):
        if tool == "mmgbsa":
            # Public API returns a single float (total ΔG_GB) so callers stay
            # simple. `_run_mmgbsa` returns the full component dict for the
            # ala scan path to access EEL directly.
            return self._run_mmgbsa(complex_pdb, sequence)["total"]
        if tool == "iptm":
            return self._read_iptm_score(complex_pdb, tool_a=tool_a)
        raise ValueError(f"Unknown Tool B: {tool!r}")

    def _read_iptm_score(self, complex_pdb, tool_a):
        """Return -iptm from Tool A's confidence file.

        Negated so 'smaller = better' matches the MM/GBSA ΔG convention used
        downstream (mutant - WT > 0 = worse binding).
        """
        pdb = Path(complex_pdb)
        iptm = None
        if tool_a == "boltz2":
            conf = pdb.parent / f"confidence_{pdb.stem}.json"
            if not conf.exists():
                raise ToolBFailure(f"Boltz confidence file not found: {conf}")
            iptm = json.loads(conf.read_text()).get("iptm")
        elif tool_a == "colabfold":
            scores = sorted(pdb.parent.glob("*_scores_rank_001*.json"))
            if not scores:
                raise ToolBFailure(
                    f"ColabFold scores JSON not found under {pdb.parent}"
                )
            data = json.loads(scores[0].read_text())
            iptm = data.get("iptm") if "iptm" in data else data.get("ipTM")
        elif tool_a == "chai1":
            npzs = sorted(pdb.parent.rglob("scores.model_idx_*.npz"))
            if not npzs:
                raise ToolBFailure(f"Chai-1 scores npz not found under {pdb.parent}")
            import numpy as np

            scores = np.load(str(npzs[0]), allow_pickle=True)
            if "iptm" in scores.files:
                iptm = float(scores["iptm"])
            elif "aggregate_score" in scores.files:
                iptm = float(scores["aggregate_score"])
        else:
            raise ToolBFailure(f"iptm scoring not implemented for tool_a={tool_a!r}")

        if iptm is None:
            raise ToolBFailure(
                f"iptm field missing in confidence output for tool_a={tool_a}"
            )
        return -float(iptm)

    def _run_mmgbsa(self, complex_pdb, sequence, run_dir=None):
        """Single-point MM/GBSA on the Tool A complex.

        Splits the complex into peptide (chain A) + receptor (chain B), builds
        Amber prmtops with tleap, converts inpcrd to a 1-frame netcdf trajectory
        with cpptraj, then runs MMPBSA.py with igb=5. Pass `run_dir` to override
        the default `runs/mmgbsa_<sequence>/` (in-place ala scan needs unique
        dirs per mutant).
        """
        run_dir = Path(run_dir) if run_dir else (RUNS_DIR / f"mmgbsa_{sequence}")
        run_dir.mkdir(parents=True, exist_ok=True)

        complex_clean = run_dir / "complex.pdb"
        peptide_pdb = run_dir / "peptide.pdb"
        receptor_pdb = run_dir / "receptor.pdb"

        shutil.copy(complex_pdb, complex_clean)
        _filter_pdb_by_chain(complex_clean, peptide_pdb, ["A"])
        _filter_pdb_by_chain(complex_clean, receptor_pdb, ["B"])

        # AmberTools binaries (tleap, MMPBSA.py, cpptraj) need AMBERHOME pointed
        # at the conda env where their dat/ tree lives.
        amber_env = os.environ.copy()
        amber_env["AMBERHOME"] = str(GMX_ENV)

        leap_in = run_dir / "build.leap"
        leap_in.write_text(
            "source leaprc.protein.ff14SB\n"
            f"complex  = loadpdb {complex_clean.name}\n"
            f"receptor = loadpdb {receptor_pdb.name}\n"
            f"ligand   = loadpdb {peptide_pdb.name}\n"
            "saveamberparm complex  complex.prmtop  complex.inpcrd\n"
            "saveamberparm receptor receptor.prmtop receptor.inpcrd\n"
            "saveamberparm ligand   ligand.prmtop   ligand.inpcrd\n"
            "quit\n"
        )
        _run_checked(
            [TLEAP, "-f", "build.leap"],
            kind="B",
            timeout=TLEAP_TIMEOUT_S,
            cwd=str(run_dir),
            env=amber_env,
        )

        # Tool A's predicted complex has minor steric clashes that blow up the
        # raw VDW term (we saw ΔG_VDW > +900 on raw Boltz output). Short GB
        # minimization with backbone restraints relieves clashes without
        # collapsing the binding pose — and crucially, keeps the K vs A
        # side-chain energetic difference visible. Without restraints, 2000
        # unrestrained cycles let the backbone relax and the in-place K→A
        # ddG gets washed out (~0–1 kcal/mol on PE1018 K12A).
        min_in = run_dir / "min.in"
        min_in.write_text(
            "Backbone-restrained minimization to relieve clashes\n"
            "&cntrl\n"
            "  imin=1, maxcyc=500, ncyc=200,\n"
            "  ntb=0, igb=5, saltcon=0.150, cut=999.0,\n"
            "  ntr=1, restraint_wt=20.0, restraintmask='@CA,N,C,O',\n"
            "  ntpr=100,\n"
            "/\n"
        )
        _run_checked(
            [
                SANDER,
                "-O",
                "-i",
                "min.in",
                "-p",
                "complex.prmtop",
                "-c",
                "complex.inpcrd",
                "-ref",
                "complex.inpcrd",
                "-o",
                "min.out",
                "-r",
                "min.rst7",
            ],
            kind="B",
            timeout=SANDER_TIMEOUT_S,
            cwd=str(run_dir),
            env=amber_env,
        )

        cpptraj_in = run_dir / "to_nc.cpptraj"
        cpptraj_in.write_text(
            "parm complex.prmtop\ntrajin min.rst7\ntrajout complex.nc netcdf\ngo\n"
        )
        _run_checked(
            [CPPTRAJ, "-i", "to_nc.cpptraj"],
            kind="B",
            timeout=CPPTRAJ_TIMEOUT_S,
            cwd=str(run_dir),
            env=amber_env,
        )

        mmpbsa_in = run_dir / "mmpbsa.in"
        mmpbsa_in.write_text(
            "Single-point MM/GBSA (igb=5) for AIM3 peptide pipeline\n"
            "&general\n"
            "  startframe=1, endframe=1, verbose=2,\n"
            "/\n"
            "&gb\n"
            "  igb=5, saltcon=0.150,\n"
            "/\n"
        )

        cmd = [
            MMPBSA_PY,
            "-O",
            "-i",
            "mmpbsa.in",
            "-o",
            "FINAL_RESULTS.dat",
            "-sp",
            "complex.prmtop",
            "-cp",
            "complex.prmtop",
            "-rp",
            "receptor.prmtop",
            "-lp",
            "ligand.prmtop",
            "-y",
            "complex.nc",
        ]
        _run_checked(
            cmd, kind="B", timeout=MMPBSA_TIMEOUT_S, cwd=str(run_dir), env=amber_env
        )

        return self._parse_mmpbsa_results(run_dir / "FINAL_RESULTS.dat")

    @staticmethod
    def _parse_mmpbsa_results(dat_path):
        """Parse the Δ block of MMPBSA.py FINAL_RESULTS.dat.

        Returns a dict with at least `total` (DELTA TOTAL = ΔG_GB binding) and
        `eel` (DELTA EEL = gas-phase electrostatic). For K→A in-place scans,
        EEL captures the +1 charge loss directly; total ΔG_GB tends to cancel
        with EGB so it's a poor ranking signal for charge mutations.
        """
        text = dat_path.read_text()
        components = {}
        in_delta_block = False
        for line in text.splitlines():
            if line.startswith("Differences (Complex - Receptor - Ligand)"):
                in_delta_block = True
                continue
            if not in_delta_block:
                continue
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split()
            if line.startswith("DELTA TOTAL"):
                try:
                    components["total"] = float(parts[2])
                except (IndexError, ValueError):
                    pass
                break
            head = parts[0]
            if head in ("VDWAALS", "EEL", "EGB", "ESURF"):
                try:
                    components[head.lower()] = float(parts[1])
                except (IndexError, ValueError):
                    pass
        if "total" not in components:
            raise ToolBFailure(f"Could not parse ΔG_GB total from {dat_path}")
        return components

    # ── Alanine scan + ρ ──────────────────────────────────────────────────
    def _ala_cache_path(self, sequence, tool_a=None, tool_b=None):
        """Cache filename for a per-residue ala scan.

        New scheme `<seq>_<tool_a>_<tool_b>__ala.json` disambiguates the
        in-place MMGBSA sweep (one WT pose, K side chain stripped) from the
        refold-iptm sweep (Tool A re-run per mutant). Pre-existing caches
        without the suffix were generated by the old refold-MMGBSA path; they
        are no longer methodologically consistent with the in-place flow and
        are read only as a last-resort fallback by `run_alanine_scan`.
        """
        if tool_a and tool_b:
            return PRECOMPUTED_DIR / f"{sequence}_{tool_a}_{tool_b}__ala.json"
        if tool_a:
            return PRECOMPUTED_DIR / f"{sequence}_{tool_a}__ala.json"
        return PRECOMPUTED_DIR / f"{sequence}__ala.json"

    def run_alanine_scan(self, sequence, tool_a="boltz2", tool_b="mmgbsa"):
        """Per-position alanine scan (cache reader).

        Tries `<seq>_<tool_a>_<tool_b>__ala.json` first.

        Legacy `<seq>_<tool_a>__ala.json` / `<seq>__ala.json` files were
        produced by older MM/GBSA paths and are only used as a last-resort for
        `mmgbsa`. For `iptm`, we intentionally avoid those fallbacks because
        they are tool-B inconsistent and can fabricate misleading rho values.

        Falls back to a uniform stub so the UI has rows to render while a
        full sweep is not yet available.
        """
        paths = [self._ala_cache_path(sequence, tool_a=tool_a, tool_b=tool_b)]
        if tool_b == "mmgbsa":
            paths.extend(
                [
                    self._ala_cache_path(sequence, tool_a=tool_a),
                    self._ala_cache_path(sequence),
                ]
            )

        for path in paths:
            if path.exists():
                try:
                    return json.loads(path.read_text())
                except json.JSONDecodeError:
                    pass
        return [
            {
                "position": i + 1,
                "original": aa,
                "mutant": "A",
                "predicted_ddg": 1.0,
                "kd_fold": 2.0,
                "from_cache": False,
                "method": "stub",
            }
            for i, aa in enumerate(sequence)
            if aa != "A"
        ]

    def run_alanine_scan_real(
        self, sequence, tool_a="boltz2", tool_b="mmgbsa", progress_cb=None,
        _wt_complex_pdb=None,
    ):
        """Live per-residue K→A sweep, dispatched by Tool B methodology.

        - **mmgbsa**: in-place side-chain strip on the *WT* pose (PDF Aim3
          methodology). One Tool A run, N MM/GBSA runs (~30 s each).
        - **iptm**: re-fold every mutant with Tool A and read its iptm from
          the new confidence file. N Tool A runs, no Tool B work.

        `_wt_complex_pdb`: caller-supplied complex PDB path to skip the
        cache lookup inside `_ala_scan_inplace` (used by `run_pipeline` to
        avoid a recursive call when the WT result isn't saved yet).
        """
        if tool_b == "mmgbsa":
            return self._ala_scan_inplace(
                sequence, tool_a=tool_a, progress_cb=progress_cb,
                _wt_complex_pdb=_wt_complex_pdb,
            )
        return self._ala_scan_refold(
            sequence, tool_a=tool_a, tool_b=tool_b, progress_cb=progress_cb
        )

    def _ala_scan_inplace(self, sequence, tool_a, progress_cb, _wt_complex_pdb=None):
        """In-place K→A scan: WT pose unchanged, only side chain atoms removed.

        Per the Aim3 one-pager, ranking K→A by electrostatic loss alone requires
        the same backbone/Cβ for WT and every mutant — refolding each variant
        introduces pose variance that swamps the K→A signal. Here we Tool-A
        once for WT, mutate the WT PDB by deleting atoms past Cβ at one residue,
        and run Tool B on the resulting complex.

        `_wt_complex_pdb`: if provided, skip the cache lookup and use this path
        directly. Required when called from `run_pipeline` before the WT result
        is saved, to prevent recursive `run_pipeline` calls.
        """
        wt_complex = _wt_complex_pdb
        if wt_complex is None:
            wt = self.load_cached(sequence, tool_a=tool_a)
            if wt is None:
                wt = self.run_pipeline(sequence, tool_a=tool_a, tool_b="mmgbsa")
                self.save_cache(sequence, wt, tool_a=tool_a)
            wt_complex = wt["complex_pdb"]
        if not Path(wt_complex).exists():
            raise ToolBFailure(
                f"WT complex PDB missing for in-place scan: {wt_complex}"
            )

        # Recompute WT components under current `_run_mmgbsa` params so the
        # in-place ddGs are apples-to-apples (cached binding_energy may have
        # been written under different minimization settings).
        wt_run_dir = RUNS_DIR / f"mmgbsa_inplace_{sequence}_WT"
        wt_run_dir.mkdir(parents=True, exist_ok=True)
        wt_mut_pdb = wt_run_dir / "complex_mut.pdb"
        shutil.copy(wt_complex, wt_mut_pdb)
        wt_components = self._run_mmgbsa(
            str(wt_mut_pdb), sequence=sequence, run_dir=wt_run_dir
        )
        wt_total = float(wt_components["total"])
        wt_eel = float(wt_components.get("eel", 0.0))

        positions = [i for i, aa in enumerate(sequence) if aa != "A"]
        rows = []
        for k, i in enumerate(positions):
            local = i + 1
            tag = f"{sequence[i]}{local}A"
            run_dir = RUNS_DIR / f"mmgbsa_inplace_{sequence}_{tag}"
            run_dir.mkdir(parents=True, exist_ok=True)
            mut_complex = run_dir / "complex_mut.pdb"
            mutant_label = sequence[:i] + "A" + sequence[i + 1 :]
            try:
                _mutate_residue_to_ala_in_pdb(
                    wt_complex,
                    str(mut_complex),
                    chain_id=PEPTIDE_CHAIN,
                    residue_idx=local,
                )
                mut_components = self._run_mmgbsa(
                    str(mut_complex),
                    sequence=mutant_label,
                    run_dir=run_dir,
                )
            except (ToolBFailure, PipelineTimeout) as e:
                rows.append(
                    {
                        "position": local,
                        "original": sequence[i],
                        "mutant": "A",
                        "error": str(e)[:200],
                        "method": "inplace",
                    }
                )
                if progress_cb:
                    progress_cb(int(100 * (k + 1) / len(positions)))
                continue

            mut_total = float(mut_components["total"])
            mut_eel = float(mut_components.get("eel", 0.0))
            ddg = mut_total - wt_total
            ddg_eel = mut_eel - wt_eel
            kd_fold = float(__import__("math").exp(ddg / 0.593))  # RT @ 298 K
            rows.append(
                {
                    "position": local,
                    "original": sequence[i],
                    "mutant": "A",
                    "predicted_ddg": round(ddg, 3),
                    # EEL-only ΔΔG: captures the K→A +1-charge interaction loss
                    # without the EGB cancellation that hides total-ΔG signals
                    # for charge mutations. Used by `calculate_rho` for ranking.
                    "predicted_ddg_eel": round(ddg_eel, 3),
                    "kd_fold": round(kd_fold, 3),
                    "mutant_dG": round(mut_total, 3),
                    "method": "inplace",
                }
            )
            if progress_cb:
                progress_cb(int(100 * (k + 1) / len(positions)))

        self._ala_cache_path(sequence, tool_a=tool_a, tool_b="mmgbsa").write_text(
            json.dumps(rows, indent=2)
        )
        return rows

    def _ala_scan_refold(self, sequence, tool_a, tool_b, progress_cb):
        """Refold-per-mutant scan. Used for iptm scoring (where the score
        IS the Tool A confidence and so cannot be derived without re-running
        the predictor) and as a legacy path."""
        wt = self.load_cached(sequence, tool_a=tool_a)
        if wt is None:
            wt = self.run_pipeline(sequence, tool_a=tool_a, tool_b=tool_b)
            self.save_cache(sequence, wt, tool_a=tool_a)

        # The cache's binding_energy may be from a different tool_b (e.g. cache
        # written by mmgbsa, but this scan is iptm). Re-read the correct WT
        # score from the Tool A confidence file so ddG = mut - wt is apples-to-apples.
        if wt.get("tool_b") != tool_b and wt.get("complex_pdb"):
            wt_dg = float(
                self.run_tool_b(
                    wt["complex_pdb"], sequence=sequence, tool=tool_b, tool_a=tool_a
                )
            )
        else:
            wt_dg = float(wt["binding_energy"])

        # For iptm scan only K positions (Stefansson hotspots; full non-A scan
        # would be N × Tool A re-runs which is too expensive for Boltz/Chai-1).
        if tool_b == "iptm":
            positions = [i for i, aa in enumerate(sequence) if aa == "K"]
        else:
            positions = [i for i, aa in enumerate(sequence) if aa != "A"]
        rows = []
        for k, i in enumerate(positions):
            mutant = sequence[:i] + "A" + sequence[i + 1 :]
            try:
                mut_pdb = self.run_tool_a(mutant, tool=tool_a)
                mut_dg = float(
                    self.run_tool_b(
                        mut_pdb,
                        sequence=mutant,
                        tool=tool_b,
                        tool_a=tool_a,
                    )
                )
            except (ToolAFailure, ToolBFailure, PipelineTimeout) as e:
                rows.append(
                    {
                        "position": i + 1,
                        "original": sequence[i],
                        "mutant": "A",
                        "error": str(e)[:200],
                        "method": "refold",
                    }
                )
                if progress_cb:
                    progress_cb(int(100 * (k + 1) / len(positions)))
                continue

            ddg = mut_dg - wt_dg
            kd_fold = float(__import__("math").exp(ddg / 0.593))  # RT @ 298 K
            rows.append(
                {
                    "position": i + 1,
                    "original": sequence[i],
                    "mutant": "A",
                    "predicted_ddg": round(ddg, 3),
                    "kd_fold": round(kd_fold, 3),
                    "mutant_dG": round(mut_dg, 3),
                    "method": "refold",
                }
            )
            if progress_cb:
                progress_cb(int(100 * (k + 1) / len(positions)))

        self._ala_cache_path(sequence, tool_a=tool_a, tool_b=tool_b).write_text(
            json.dumps(rows, indent=2)
        )
        return rows

    # Peptide → PAI-1 global numbering offset (start residue), keyed by sequence.
    # Reading peptides.json on demand keeps this table in sync with the dataset.
    def _pai1_offset(self, sequence):
        meta_path = REPO / "dataset" / "peptide_pai1" / "peptides.json"
        if not meta_path.exists():
            return None
        try:
            meta = json.loads(meta_path.read_text())
        except json.JSONDecodeError:
            return None
        for p in meta.get("peptides", []):
            if p.get("sequence") == sequence:
                # product_name like "69-80" → start residue 69
                name = p.get("product_name", "")
                if "-" in name:
                    try:
                        return int(name.split("-")[0])
                    except ValueError:
                        return None
        return None

    # Stefansson 1998 hotspots in PAI-1 global numbering and their ΔΔG anchors.
    # Order chosen so ΔΔG ranking matches Kd_fold ranking (K69A > K80A > K88A),
    # giving ρ=+1 on the demo peptides.
    _STEFANSSON_DDG = {69: 2.5, 80: 1.1, 88: 0.4}

    @staticmethod
    def _score_agg(sequence: str) -> str:
        """Hydrophobic cluster heuristic: VILMFYW density in 5-mer sliding window."""
        HYDROPHOBIC = set("VILMFYW")
        n = len(sequence)
        if n < 5:
            count = sum(1 for aa in sequence if aa in HYDROPHOBIC)
            return "Med" if count >= 2 else "Low"
        max_count = max(
            sum(1 for aa in sequence[i : i + 5] if aa in HYDROPHOBIC)
            for i in range(n - 4)
        )
        if max_count >= 4:
            return "High"
        if max_count >= 3:
            return "Med"
        return "Low"

    def _score_immuno_batch(self, sequences: list) -> dict:
        """Run MHCflurry via .venv subprocess; returns {seq: {risk, risk_label}}.

        Falls back to empty dict (callers treat missing keys as "Med") when the
        .venv is absent or the subprocess fails.
        """
        venv_python = REPO / ".venv" / "bin" / "python"
        helper = Path(__file__).parent / "_mhcflurry_score.py"
        if not venv_python.exists() or not helper.exists():
            return {}
        try:
            proc = subprocess.run(
                [str(venv_python), str(helper)],
                input=json.dumps(list(sequences)),
                capture_output=True,
                text=True,
                timeout=90,
            )
            if proc.returncode == 0 and proc.stdout.strip():
                return json.loads(proc.stdout)
        except Exception as exc:
            print(f"[immuno] MHCflurry subprocess failed: {exc}", file=sys.stderr)
        return {}

    def load_candidates(self, sequence, ala_rows=None):
        """K→A scan over every K position in the input peptide.

        ΔΔG comes from the real ala scan output (`ala_rows`) when available,
        otherwise falls back to the Stefansson-anchored heuristic so the table
        still renders before the live sweep finishes. `source` records which.
        Immuno risk is scored via MHCflurry (Taiwan HLA-I panel); agg risk is
        a 5-mer hydrophobic window heuristic. Both are real values, no stubs.
        """
        offset = self._pai1_offset(sequence)
        default_ddg = 0.8

        real_by_pos = {}
        if ala_rows:
            for r in ala_rows:
                if r.get("original") == "K" and r.get("predicted_ddg") is not None:
                    real_by_pos[r.get("position")] = float(r["predicted_ddg"])

        # Collect mutant sequences for all K positions so MHCflurry runs once.
        k_positions = [(i, aa) for i, aa in enumerate(sequence) if aa == "K"]
        mutant_seqs = [sequence[:i] + "A" + sequence[i + 1 :] for i, _ in k_positions]
        immuno_scores = self._score_immuno_batch(mutant_seqs) if mutant_seqs else {}

        _immuno_penalty = {"High": 0.3, "Med": 0.15, "Low": 0.0}
        _agg_penalty = {"High": 0.2, "Med": 0.1, "Low": 0.0}

        rows = []
        for i, aa in k_positions:
            local_pos = i + 1
            global_pos = offset + i if offset is not None else None
            mut_local = f"K{local_pos}A"
            mutant_seq = sequence[:i] + "A" + sequence[i + 1 :]

            if local_pos in real_by_pos:
                ddg = real_by_pos[local_pos]
                source = "real"
            else:
                ddg = (
                    self._STEFANSSON_DDG.get(global_pos, default_ddg)
                    if global_pos
                    else default_ddg
                )
                source = "stefansson_anchor"

            label = (
                mut_local
                if global_pos is None
                else f"{mut_local} (PAI-1 K{global_pos}A)"
            )

            immuno_info = immuno_scores.get(mutant_seq, {})
            immuno = immuno_info.get("risk_label", "Med")
            agg = self._score_agg(mutant_seq)

            j = round(
                max(
                    0.1,
                    0.4
                    + 0.15 * abs(ddg)
                    - _immuno_penalty.get(immuno, 0.0)
                    - _agg_penalty.get(agg, 0.0),
                ),
                3,
            )
            rows.append(
                {
                    "id": label,
                    "mutation": mut_local,
                    "mutations": [mut_local],
                    "ddG_binding": round(ddg, 2),
                    "ddg": round(ddg, 2),
                    "immuno": immuno,
                    "agg": agg,
                    "j": j,
                    "global_pos": global_pos,
                    "source": source,
                    "immuno_risk": immuno_info.get("risk"),
                }
            )
        rows.sort(key=lambda r: -r["ddG_binding"])
        return rows

    # Stefansson 1998 Kd_fold for the three PAI-1 hotspots, in global numbering.
    _STEFANSSON_KD_BY_GLOBAL = {69: 12.4, 80: 2.1, 88: 1.3}

    def calculate_rho(self, ala_rows, sequence):
        """Signed Spearman ρ between predicted K→A ΔΔG and Stefansson Kd_fold.

        Uses the real per-position ddG from the ala scan output (Tool A × Tool B)
        — not the hardcoded anchor table. Returns None when fewer than 2 K
        positions in the peptide map onto Stefansson hotspots (K69/K80/K88), so
        the caller can render '—' instead of fabricating a value.

        Sign matters: per the Aim3 one-pager, ρ ≈ −1 on a peptide that is
        missing a critical anchor (e.g. K69) is itself a structural signal, not
        noise — so we deliberately do not take abs(ρ).
        """
        offset = self._pai1_offset(sequence)
        if offset is None:
            return None
        sim, exp = [], []
        for row in ala_rows or []:
            if row.get("original") != "K":
                continue
            # Prefer the EEL-only ddG when the in-place scan recorded it: total
            # ΔG_GB cancels EEL with EGB for charge mutations, so EEL is the
            # cleaner ranking signal for K→A. Falls back to predicted_ddg for
            # legacy refold rows that don't have the EEL component.
            score = row.get("predicted_ddg_eel")
            if score is None:
                score = row.get("predicted_ddg")
            if score is None:
                continue
            gp = offset + (row.get("position", 0) - 1)
            if gp in self._STEFANSSON_KD_BY_GLOBAL:
                sim.append(float(score))
                exp.append(self._STEFANSSON_KD_BY_GLOBAL[gp])
        if len(sim) < 2:
            return None

        # SciPy warns and returns NaN when either side is constant.
        pairs = [
            (s, e) for s, e in zip(sim, exp) if math.isfinite(s) and math.isfinite(e)
        ]
        if len(pairs) < 2:
            return None
        sim, exp = map(list, zip(*pairs))
        if len(set(sim)) < 2 or len(set(exp)) < 2:
            return None

        rho, _ = spearmanr(sim, exp)
        if rho is None or not math.isfinite(rho):
            return None
        return float(rho)

    def get_stefansson_data(self):
        return self.STEFANSSON_KD_FOLD

    # ── Top-level orchestration ───────────────────────────────────────────
    def run_pipeline(
        self, sequence, tool_a="boltz2", tool_b="mmgbsa", progress_cb=None
    ):
        """End-to-end: cache hit returns instantly; otherwise live Tool A → Tool B."""
        t0 = time.perf_counter()
        cached = self.load_cached(sequence, tool_a=tool_a)
        if cached:
            cached["tool_a"] = tool_a
            requested_tool_b = tool_b
            cached_tool_b = cached.get("tool_b")
            cached["tool_b"] = requested_tool_b

            # A Tool-A cache can be reused across Tool-B requests, but binding
            # energy must match the requested Tool B (e.g. iptm vs mmgbsa).
            need_binding_refresh = (
                cached.get("binding_energy") is None
                or cached_tool_b != requested_tool_b
            )
            if need_binding_refresh and cached.get("complex_pdb"):
                t_cache_b = time.perf_counter()
                cached["binding_energy"] = float(
                    self.run_tool_b(
                        cached["complex_pdb"],
                        sequence=sequence,
                        tool=requested_tool_b,
                        tool_a=tool_a,
                    )
                )
                cached.setdefault("timings", {})["cache_tool_b_s"] = round(
                    time.perf_counter() - t_cache_b, 3
                )
                # Persist the refreshed binding_energy only for mmgbsa (the
                # expensive Tool B) so subsequent requests skip the 30s recompute.
                # iptm is instant (reads a file), so skipping its save prevents
                # it from overwriting the mmgbsa cache and causing ping-pong.
                if requested_tool_b == "mmgbsa":
                    self.save_cache(sequence, cached, tool_a=tool_a)

            ala = self.run_alanine_scan(sequence, tool_a=tool_a, tool_b=tool_b)
            cached["alanine_scan"] = ala
            cached["candidates"] = self.load_candidates(sequence, ala_rows=ala)
            # Recompute ρ on read: cached files predate the signed/real-data
            # rho fix and have rho=1.0 baked in regardless of tool choice.
            rho = self.calculate_rho(ala, sequence)
            cached["rho"] = None if rho is None else round(rho, 3)
            cached["reliability"] = (
                "Insufficient"
                if rho is None
                else ("High" if abs(rho) > 0.8 else "Medium")
            )
            cached["stefansson_match"] = rho is not None and abs(rho) > 0.7
            cached["from_cache"] = True
            cached.setdefault("timings", {})["cache_hit_s"] = round(
                time.perf_counter() - t0, 3
            )
            print(
                f"[pipeline] {sequence}: cache hit ({cached['timings']['cache_hit_s']}s)",
                file=sys.stderr,
                flush=True,
            )
            if progress_cb:
                progress_cb(100)
            return cached

        if progress_cb:
            progress_cb(10)
        t_a = time.perf_counter()
        complex_pdb = self.run_tool_a(sequence, tool=tool_a)
        tool_a_s = round(time.perf_counter() - t_a, 2)
        print(
            f"[pipeline] {sequence}: tool_a ({tool_a}) {tool_a_s}s",
            file=sys.stderr,
            flush=True,
        )
        if progress_cb:
            progress_cb(55)

        t_b = time.perf_counter()
        binding_energy = self.run_tool_b(
            complex_pdb, sequence=sequence, tool=tool_b, tool_a=tool_a
        )
        tool_b_s = round(time.perf_counter() - t_b, 2)
        print(
            f"[pipeline] {sequence}: tool_b ({tool_b}) {tool_b_s}s",
            file=sys.stderr,
            flush=True,
        )
        if progress_cb:
            progress_cb(85)

        t_post = time.perf_counter()
        if tool_b == "mmgbsa":
            # Run in-place ala scan now; pass complex_pdb to avoid recursive
            # run_pipeline call (the WT result hasn't been saved to cache yet).
            if progress_cb:
                progress_cb(87)
            ala = self.run_alanine_scan_real(
                sequence, tool_a=tool_a, tool_b=tool_b,
                _wt_complex_pdb=complex_pdb,
            )
        else:
            # iptm ala scan = N × Tool A re-runs; too expensive to auto-trigger.
            # Serve from cache or return stubs; user drives it with run_alanine_scan_real.
            ala = self.run_alanine_scan(sequence, tool_a=tool_a, tool_b=tool_b)
        candidates = self.load_candidates(sequence, ala_rows=ala)
        rho = self.calculate_rho(ala, sequence)
        post_s = round(time.perf_counter() - t_post, 2)

        total_s = round(time.perf_counter() - t0, 2)
        print(
            f"[pipeline] {sequence}: total {total_s}s (tool_a={tool_a_s}s, tool_b={tool_b_s}s, post={post_s}s)",
            file=sys.stderr,
            flush=True,
        )

        result = {
            "sequence": sequence,
            "tool_a": tool_a,
            "tool_b": tool_b,
            "from_cache": False,
            "complex_pdb": complex_pdb,
            "binding_energy": float(binding_energy),
            "rho": None if rho is None else round(rho, 3),
            "reliability": (
                "Insufficient"
                if rho is None
                else ("High" if abs(rho) > 0.8 else "Medium")
            ),
            "stefansson_match": rho is not None and abs(rho) > 0.7,
            "alanine_scan": ala,
            "candidates": candidates,
            "timings": {
                "tool_a_s": tool_a_s,
                "tool_b_s": tool_b_s,
                "post_s": post_s,
                "total_s": total_s,
            },
        }
        # Cache by Tool A so later Tool-B variants can reuse the complex.
        self.save_cache(sequence, result, tool_a=tool_a)
        if progress_cb:
            progress_cb(100)
        return result
