"""
Build the LRP1 receptor structure for peptide docking (子計畫 3 Year-1 S2).

Source
------
PDB 1J8E — LRP1 Complement-Type Repeat 7 (CR7), 1.85 Å X-ray, 44 residues
(Simonovic et al., Biochemistry 2001, 40:15127). Bound Ca²⁺ is retained —
it is structurally essential for the CR fold β-hairpin.

Why CR7 alone (and not a CR tandem)?
-----------------------------------
PAI-1 binds LRP1 across CR cluster II (CR3-CR10). Single-domain crystals
dominate the available PDB entries; no public co-crystal of PAI-1 + LRP1 exists.
For Year-1 baseline we use CR7 (the best-resolution LRP1 CR crystal and one of
the known PAI-1 contact repeats) as a minimal receptor. A multi-CR tandem or
AlphaFold2-predicted CR3-CR10 model is a Year-2 upgrade path documented in
dataset/lrp1_cr2/metadata.json.

Pipeline
--------
1. PDBFixer: strip waters, keep Ca²⁺, add missing atoms + hydrogens.
2. OpenMM Amber14 + OBC2 implicit solvent, L-BFGS minimization.
3. Heavy-atom RMSD vs. crystal must stay < 2.0 Å (sanity check against
   over-minimization that would flatten the binding surface).

Outputs
-------
dataset/lrp1_cr2/receptor.pdb      Cleaned + minimized, single chain + Ca²⁺
dataset/lrp1_cr2/metadata.json     Provenance, minimization stats, schema
"""

import json
from pathlib import Path

from pdbfixer import PDBFixer
from openmm import LangevinIntegrator, Platform, unit, CustomExternalForce  # noqa: F401
from openmm.app import ForceField, Modeller, PDBFile, Simulation, HBonds

ROOT = Path(__file__).parent.parent
OUT_DIR = ROOT / "dataset" / "lrp1_cr2"
RAW_PDB = OUT_DIR / "raw" / "1J8E.pdb"
OUT_PDB = OUT_DIR / "receptor.pdb"
META_PATH = OUT_DIR / "metadata.json"

RMSD_LIMIT_ANG = 2.0
MIN_TOLERANCE_KJ_MOL = 10.0
MAX_MIN_ITER = 2000


def _heavy_atom_rmsd(coords_a, coords_b):
    import numpy as np
    a = np.asarray(coords_a, dtype=float)
    b = np.asarray(coords_b, dtype=float)
    if a.shape != b.shape:
        raise ValueError(f"Atom count mismatch: {a.shape} vs {b.shape}")
    return float(((a - b) ** 2).sum(axis=1).mean() ** 0.5)


def _collect_heavy_coords(topology, positions_nm):
    import numpy as np
    coords = []
    pos = np.asarray(positions_nm.value_in_unit(unit.angstrom))
    for atom in topology.atoms():
        if atom.element is not None and atom.element.symbol != "H":
            coords.append(pos[atom.index])
    return coords


def build():
    if not RAW_PDB.exists():
        raise FileNotFoundError(f"Expected raw PDB at {RAW_PDB}")

    print(f"[build_lrp1_receptor] fixing {RAW_PDB.name}")
    fixer = PDBFixer(filename=str(RAW_PDB))
    fixer.removeHeterogens(keepWater=False)  # strips water; Ca²⁺ re-added below
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)

    # Re-inject Ca²⁺ (removeHeterogens strips it; it's structural, not a ligand).
    import numpy as np
    ca_line = None
    for line in RAW_PDB.read_text().splitlines():
        if line.startswith("HETATM") and line[17:20].strip() == "CA" and line[12:16].strip() == "CA":
            ca_line = line
            break
    if ca_line is None:
        raise RuntimeError("No Ca²⁺ found in raw 1J8E — unexpected")
    ca_xyz_ang = [float(ca_line[30:38]), float(ca_line[38:46]), float(ca_line[46:54])]

    modeller = Modeller(fixer.topology, fixer.positions)
    # Append Ca²⁺ as a separate chain/residue with its own topology entry
    from openmm.app import Element
    ca_chain = modeller.topology.addChain(id="C")
    ca_res = modeller.topology.addResidue("CA", ca_chain)
    ca_atom = modeller.topology.addAtom("CA", Element.getByAtomicNumber(20), ca_res)  # noqa: F841
    positions_nm = list(modeller.positions.value_in_unit(unit.nanometer))
    positions_nm.append((ca_xyz_ang[0] / 10.0, ca_xyz_ang[1] / 10.0, ca_xyz_ang[2] / 10.0))
    modeller.positions = positions_nm * unit.nanometer

    print(f"[build_lrp1_receptor] topology: {modeller.topology.getNumResidues()} residues, "
          f"{modeller.topology.getNumAtoms()} atoms")

    # Amber14 protein + TIP3P ion parameters (needed for Ca²⁺) + OBC2 implicit solvent.
    # amber14/tip3p.xml only contributes water + ion templates; with no HOH residues
    # in the topology, only the Ca²⁺ template gets matched.
    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml", "implicit/obc2.xml")
    system = forcefield.createSystem(modeller.topology, constraints=HBonds,
                                     nonbondedCutoff=2.0 * unit.nanometer)

    integrator = LangevinIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond,
                                    0.002 * unit.picoseconds)
    try:
        platform = Platform.getPlatformByName("CUDA")
    except Exception:
        try:
            platform = Platform.getPlatformByName("OpenCL")
        except Exception:
            platform = Platform.getPlatformByName("CPU")
    print(f"[build_lrp1_receptor] platform = {platform.getName()}")

    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    pre_state = simulation.context.getState(getEnergy=True, getPositions=True)
    pre_energy = pre_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    pre_coords = _collect_heavy_coords(modeller.topology, pre_state.getPositions())
    print(f"[build_lrp1_receptor] pre-min energy  = {pre_energy:.2f} kJ/mol")

    simulation.minimizeEnergy(tolerance=MIN_TOLERANCE_KJ_MOL * unit.kilojoule_per_mole / unit.nanometer,
                              maxIterations=MAX_MIN_ITER)

    post_state = simulation.context.getState(getEnergy=True, getPositions=True)
    post_energy = post_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    post_coords = _collect_heavy_coords(modeller.topology, post_state.getPositions())
    rmsd_ang = _heavy_atom_rmsd(pre_coords, post_coords)
    print(f"[build_lrp1_receptor] post-min energy = {post_energy:.2f} kJ/mol "
          f"(Δ = {post_energy - pre_energy:+.2f})")
    print(f"[build_lrp1_receptor] heavy-atom RMSD = {rmsd_ang:.3f} Å")

    if rmsd_ang > RMSD_LIMIT_ANG:
        raise RuntimeError(
            f"RMSD {rmsd_ang:.3f} Å exceeds {RMSD_LIMIT_ANG} Å cap — "
            "minimization distorted the fold; review forcefield / Ca²⁺ coordination"
        )

    with OUT_PDB.open("w") as f:
        PDBFile.writeFile(modeller.topology, post_state.getPositions(), f, keepIds=True)

    metadata = {
        "receptor_id":   "lrp1_cr7",
        "source_pdb":    "1J8E",
        "pdb_title":     "Crystal structure of ligand-binding repeat CR7 from LRP",
        "resolution_A":  1.85,
        "chain":         "A",
        "n_residues":    modeller.topology.getNumResidues() - 1,  # minus Ca²⁺ pseudo-residue
        "n_atoms":       modeller.topology.getNumAtoms(),
        "kept_ions":     ["CA (Ca²⁺, structural, β-hairpin coordination)"],
        "removed":       ["HOH (crystallographic water)"],
        "forcefield":    ["amber14-all.xml", "amber14/tip3p.xml (ions only)", "implicit/obc2.xml"],
        "platform":      platform.getName(),
        "minimization": {
            "pre_energy_kJ_mol":  round(pre_energy, 2),
            "post_energy_kJ_mol": round(post_energy, 2),
            "tolerance_kJ_mol_nm": MIN_TOLERANCE_KJ_MOL,
            "max_iterations":     MAX_MIN_ITER,
            "heavy_atom_rmsd_A":  round(rmsd_ang, 3),
            "rmsd_cap_A":         RMSD_LIMIT_ANG,
        },
        "note": (
            "Year-1 baseline uses the single CR7 domain (1J8E) as minimal "
            "receptor. PAI-1 binds LRP1 across CR cluster II (CR3-CR10); a "
            "multi-CR tandem or AlphaFold2-predicted CR3-CR10 model is the "
            "Year-2 upgrade path."
        ),
        "sources": {
            "raw_pdb":      "raw/1J8E.pdb",
            "citation":     "Simonovic et al., Biochemistry 2001, 40:15127. DOI:10.1021/bi015688m",
        },
    }
    META_PATH.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[build_lrp1_receptor] wrote {OUT_PDB.relative_to(ROOT)} + {META_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    build()
