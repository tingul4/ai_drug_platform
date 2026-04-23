"""
Parse PDBbind v2020 R1 INDEX files into a SQLite affinity lookup table.

Input  : dataset/PDBbind_v2020_R1/index/index/INDEX_general_{PL,PP}.2020R1.lst
Output : `pdbbind_affinity` table in db/skempi.db

Why we need it:
  The analyzer's ΔΔG model is relative — ΔpKd is just a shift. To report an
  ABSOLUTE pKd for each candidate (as PDF Aim3 Summary Table requires), we
  need a WT reference. PDBbind gives us experimental Kd/Ki/IC50 for 21,835
  PDB complexes (19,037 PL + 2,798 PP). When our analyzer matches a SKEMPI
  PDB to the user's sequence, we look up that PDB's WT affinity and anchor
  the shift to a real number.

Line format (both PL and PP lists):
    <PDB>  <resolution>  <year>  <affinity>   // <notes...>
  e.g.  1fc2  2.80  1981  Kd=22.5nM     // 1fc2.pdf (C|D) ...
        5tln  2.30  1982  Ki=0.43uM     // 5tln.pdf (BAN) ...
        4cts  2.90  1984  Kd<10uM       // 4cts.pdf (OAA) ...

Affinity tokens: Kd / Ki / IC50, relation ∈ {=, <, >, ~, <=, >=},
unit ∈ {fM, pM, nM, uM, mM, M}.
"""

import re
import sqlite3
import math
from pathlib import Path

DB_PATH  = Path(__file__).parent / "skempi.db"
IDX_DIR  = Path(__file__).parent.parent / "dataset" / "PDBbind_v2020_R1" / "index" / "index"

UNITS = {"fM": 1e-15, "pM": 1e-12, "nM": 1e-9,
         "uM": 1e-6,  "mM": 1e-3,  "M":  1.0}

# 1fc2  2.80  1981  Kd=22.5nM     // 1fc2.pdf (C|D) ...
LINE_RE = re.compile(
    r"^(?P<pdb>[0-9a-z]{4})\s+"
    r"(?P<res>\S+)\s+"
    r"(?P<year>\d{4})\s+"
    r"(?P<kind>Kd|Ki|IC50)"
    r"(?P<rel><=|>=|<|>|~|=)"
    r"(?P<val>[0-9.]+)"
    r"(?P<unit>fM|pM|nM|uM|mM|M)"
    r"\s*//\s*(?P<notes>.*)$",
    re.IGNORECASE,
)


def parse_line(line: str):
    line = line.strip()
    if not line or line.startswith("#"):
        return None
    m = LINE_RE.match(line)
    if not m:
        return None
    d = m.groupdict()
    try:
        val_m = float(d["val"]) * UNITS[d["unit"]]
        if val_m <= 0:
            return None
        pkd = -math.log10(val_m)
    except Exception:
        return None
    return {
        "pdb_id":         d["pdb"].upper(),
        "kind":           d["kind"].capitalize() if d["kind"].lower() != "ic50" else "IC50",
        "relation":       d["rel"],
        "value_molar":    val_m,
        "pKd":            round(pkd, 3),
        "resolution":     None if d["res"].upper() == "NMR" else (
                             float(d["res"]) if _isfloat(d["res"]) else None),
        "year":           int(d["year"]),
        "notes":          d["notes"][:200],
    }


def _isfloat(x):
    try: float(x); return True
    except: return False


def load_index_file(path: Path, complex_type: str):
    out = []
    if not path.exists():
        print(f"[PDBbind] missing {path}")
        return out
    for ln in path.read_text(errors="ignore").splitlines():
        rec = parse_line(ln)
        if rec:
            rec["complex_type"] = complex_type   # 'PL' or 'PP'
            out.append(rec)
    print(f"[PDBbind] {path.name}: parsed {len(out)} entries")
    return out


def main():
    conn = sqlite3.connect(DB_PATH)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS pdbbind_affinity (
            pdb_id        TEXT,
            complex_type  TEXT,      -- 'PL' or 'PP'
            kind          TEXT,      -- 'Kd' / 'Ki' / 'IC50'
            relation      TEXT,      -- '=', '<', '>', '~', '<=', '>='
            value_molar   REAL,      -- affinity in mol/L
            pKd           REAL,      -- −log10(value_molar)
            resolution    REAL,
            year          INTEGER,
            notes         TEXT,
            PRIMARY KEY (pdb_id, complex_type, kind)
        )
    """)
    conn.execute("CREATE INDEX IF NOT EXISTS idx_pdbbind_pdb ON pdbbind_affinity(pdb_id)")
    # Clean slate so repeat builds stay idempotent
    conn.execute("DELETE FROM pdbbind_affinity")

    records  = []
    records += load_index_file(IDX_DIR / "INDEX_general_PP.2020R1.lst", "PP")
    records += load_index_file(IDX_DIR / "INDEX_general_PL.2020R1.lst", "PL")

    conn.executemany("""
        INSERT OR REPLACE INTO pdbbind_affinity
          (pdb_id, complex_type, kind, relation, value_molar, pKd,
           resolution, year, notes)
        VALUES (:pdb_id, :complex_type, :kind, :relation, :value_molar, :pKd,
                :resolution, :year, :notes)
    """, records)
    conn.commit()

    # Coverage vs SKEMPI
    n_total    = conn.execute("SELECT COUNT(*) FROM pdbbind_affinity").fetchone()[0]
    n_pp       = conn.execute("SELECT COUNT(*) FROM pdbbind_affinity WHERE complex_type='PP'").fetchone()[0]
    n_pl       = conn.execute("SELECT COUNT(*) FROM pdbbind_affinity WHERE complex_type='PL'").fetchone()[0]
    n_skempi_covered = conn.execute("""
        SELECT COUNT(DISTINCT s.pdb_id)
        FROM skempi_mutations s
        JOIN pdbbind_affinity p ON s.pdb_id = p.pdb_id
    """).fetchone()[0]
    n_skempi_total   = conn.execute("SELECT COUNT(DISTINCT pdb_id) FROM skempi_mutations").fetchone()[0]

    print(f"[PDBbind] total entries      = {n_total}  (PP={n_pp}, PL={n_pl})")
    print(f"[PDBbind] SKEMPI PDBs w/ Kd  = {n_skempi_covered} / {n_skempi_total}")
    conn.close()


if __name__ == "__main__":
    main()
