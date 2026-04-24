"""
Build the SQLite database from real SKEMPI 2.0 data + PDB structures.
Schema:
  - skempi_mutations     : 7085 real experimental records
  - pdb_structures       : metadata for 345 PDB structures
  - pdb_residues         : per-residue info parsed from PDB + mapping files
  - analysis_sessions    : user analysis runs
  - candidates           : generated optimized candidates
  - optimization_runs    : pareto optimization results
"""

import sqlite3, csv, math, os, json, hashlib
from pathlib import Path
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

DB_PATH = Path(__file__).parent / "skempi.db"
PROJECT_ROOT = Path(__file__).parent.parent
SKEMPI_DIR = PROJECT_ROOT / "dataset" / "skempi_v2"
CSV_PATH = SKEMPI_DIR / "skempi_v2.csv"
PDB_DIR = SKEMPI_DIR / "PDBs"

# ── helpers ──────────────────────────────────────────────────────────────────
THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "MSE": "M",
    "HSD": "H",
    "HSE": "H",
    "HSP": "H",
}


def safe_float(v):
    try:
        f = float(v)
        return None if (math.isnan(f) or math.isinf(f)) else f
    except:
        return None


def calc_ddG(kd_mut, kd_wt, T=298.15):
    """ΔΔG = RT·ln(Kd_mut/Kd_wt)  [kcal/mol]"""
    try:
        RT = 1.987e-3 * T  # kcal/mol
        if kd_mut and kd_wt and kd_mut > 0 and kd_wt > 0:
            return RT * math.log(kd_mut / kd_wt)
    except:
        pass
    return None


def parse_sequence_from_pdb(pdb_path):
    """Extract chain sequences from PDB file."""
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("s", str(pdb_path))
        chains = {}
        for model in structure:
            for chain in model:
                seq = ""
                for res in chain:
                    if res.id[0] == " ":  # standard residues only
                        aa = THREE_TO_ONE.get(res.resname, "X")
                        seq += aa
                if seq:
                    chains[chain.id] = seq
            break  # first model only
        return chains
    except Exception as e:
        return {}


def parse_mapping(mapping_path):
    """Parse .mapping file → {chain: [(resname, chain, pdb_num, seq_num), ...]}"""
    chains = {}
    try:
        with open(mapping_path) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 4:
                    resname, chain, pdb_num, seq_num = (
                        parts[0],
                        parts[1],
                        int(parts[2]),
                        int(parts[3]),
                    )
                    chains.setdefault(chain, []).append(
                        (resname, chain, pdb_num, seq_num)
                    )
    except:
        pass
    return chains


# ── create schema ─────────────────────────────────────────────────────────────
def create_schema(conn):
    conn.executescript("""
    CREATE TABLE IF NOT EXISTS skempi_mutations (
        id              INTEGER PRIMARY KEY,
        pdb_complex     TEXT,
        pdb_id          TEXT,
        chain1          TEXT,
        chain2          TEXT,
        mutation_pdb    TEXT,
        mutation_clean  TEXT,
        location        TEXT,
        protein1        TEXT,
        protein2        TEXT,
        affinity_mut    REAL,
        affinity_wt     REAL,
        ddG_kcal        REAL,
        kon_mut         REAL,
        kon_wt          REAL,
        koff_mut        REAL,
        koff_wt         REAL,
        koff_ratio      REAL,
        dH_mut          REAL,
        dH_wt           REAL,
        temperature     REAL,
        method          TEXT,
        reference       TEXT,
        skempi_version  INTEGER
    );

    CREATE TABLE IF NOT EXISTS pdb_structures (
        pdb_id          TEXT PRIMARY KEY,
        chains          TEXT,      -- JSON list of chain IDs
        sequences       TEXT,      -- JSON {chain: sequence}
        n_residues      INTEGER,
        n_mutations     INTEGER,
        proteins        TEXT,      -- JSON list of protein names from SKEMPI
        created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    );

    CREATE TABLE IF NOT EXISTS pdb_residues (
        id              INTEGER PRIMARY KEY,
        pdb_id          TEXT,
        chain           TEXT,
        pdb_resnum      INTEGER,
        seq_num         INTEGER,
        residue_aa      TEXT,
        is_interface    INTEGER DEFAULT 0,
        FOREIGN KEY(pdb_id) REFERENCES pdb_structures(pdb_id)
    );

    CREATE TABLE IF NOT EXISTS analysis_sessions (
        id              INTEGER PRIMARY KEY,
        session_id      TEXT UNIQUE,
        created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
        input_sequence  TEXT,
        input_name      TEXT,
        input_target    TEXT,
        sequence_hash   TEXT,
        matched_pdb     TEXT,
        n_candidates    INTEGER,
        status          TEXT DEFAULT 'pending',
        runtime_sec     REAL,
        result_summary  TEXT   -- JSON
    );

    CREATE TABLE IF NOT EXISTS candidates (
        id              INTEGER PRIMARY KEY,
        session_id      TEXT,
        rank            INTEGER,
        variant_id      TEXT,
        mutations       TEXT,    -- JSON list of mutation strings
        ddG_binding     REAL,
        ddG_stability   REAL,
        pKd_shift       REAL,    -- ΔpKd derived from ΔΔG: −ΔΔG/(RT·ln10); ΔpIC50 ≈ ΔpKd under 1:1 competitive
        pKd_abs         REAL,    -- PDBbind-anchored absolute pKd (if WT Kd known)
        koff_relative   REAL,
        residence_time  REAL,
        immunogenicity  REAL,
        aggregation     REAL,
        solubility      REAL,
        pI              REAL,
        mw_kDa          REAL,
        hbonds          INTEGER,
        pareto_rank     INTEGER,
        score_composite REAL,
        skempi_support  INTEGER,  -- num supporting SKEMPI records
        created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
        FOREIGN KEY(session_id) REFERENCES analysis_sessions(session_id)
    );

    CREATE TABLE IF NOT EXISTS optimization_runs (
        id              INTEGER PRIMARY KEY,
        session_id      TEXT,
        round_num       INTEGER,
        n_evaluated     INTEGER,
        best_ddG        REAL,
        best_koff       REAL,
        model_r2        REAL,
        notes           TEXT,
        created_at      TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    );

    CREATE INDEX IF NOT EXISTS idx_mutations_pdb  ON skempi_mutations(pdb_id);
    CREATE INDEX IF NOT EXISTS idx_mutations_ddG  ON skempi_mutations(ddG_kcal);
    CREATE INDEX IF NOT EXISTS idx_candidates_session ON candidates(session_id);
    CREATE INDEX IF NOT EXISTS idx_residues_pdb   ON pdb_residues(pdb_id, chain);
    """)
    conn.commit()


# ── load SKEMPI CSV ────────────────────────────────────────────────────────────
def load_skempi(conn):
    print("Loading SKEMPI 2.0 data...")
    rows_inserted = 0
    with open(CSV_PATH, newline="") as f:
        reader = csv.DictReader(f, delimiter=";")
        batch = []
        for row in reader:
            pdb_complex = row["#Pdb"]
            parts = pdb_complex.split("_")
            pdb_id = parts[0]
            chain1 = parts[1] if len(parts) > 1 else ""
            chain2 = parts[2] if len(parts) > 2 else ""

            kd_mut = safe_float(row["Affinity_mut_parsed"])
            kd_wt = safe_float(row["Affinity_wt_parsed"])
            ddG = calc_ddG(kd_mut, kd_wt, safe_float(row["Temperature"]) or 298.15)

            koff_mut = safe_float(row["koff_mut_parsed"])
            koff_wt = safe_float(row["koff_wt_parsed"])
            koff_ratio = None
            if koff_mut and koff_wt and koff_wt > 0:
                koff_ratio = koff_mut / koff_wt

            batch.append(
                (
                    pdb_complex,
                    pdb_id,
                    chain1,
                    chain2,
                    row["Mutation(s)_PDB"],
                    row["Mutation(s)_cleaned"],
                    row.get("iMutation_Location(s)", ""),
                    row.get("Protein 1", ""),
                    row.get("Protein 2", ""),
                    kd_mut,
                    kd_wt,
                    ddG,
                    safe_float(row["kon_mut_parsed"]),
                    safe_float(row["kon_wt_parsed"]),
                    koff_mut,
                    koff_wt,
                    koff_ratio,
                    safe_float(row.get("dH_mut (kcal mol^(-1))")),
                    safe_float(row.get("dH_wt (kcal mol^(-1))")),
                    safe_float(row["Temperature"]),
                    row.get("Method", ""),
                    row.get("Reference", ""),
                    int(row.get("SKEMPI version", 1) or 1),
                )
            )

    conn.executemany(
        """
        INSERT INTO skempi_mutations
          (pdb_complex, pdb_id, chain1, chain2,
           mutation_pdb, mutation_clean, location,
           protein1, protein2,
           affinity_mut, affinity_wt, ddG_kcal,
           kon_mut, kon_wt, koff_mut, koff_wt, koff_ratio,
           dH_mut, dH_wt, temperature, method, reference, skempi_version)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
    """,
        batch,
    )
    conn.commit()
    rows_inserted = len(batch)
    print(f"  Inserted {rows_inserted} SKEMPI records")
    return rows_inserted


# ── load PDB structures ────────────────────────────────────────────────────────
def load_pdb_structures(conn):
    print("Parsing PDB structures...")
    pdb_files = sorted(PDB_DIR.glob("*.pdb"))
    print(f"  Found {len(pdb_files)} PDB files")

    # Get protein names from SKEMPI for each PDB
    pdb_proteins = {}
    cur = conn.execute(
        "SELECT DISTINCT pdb_id, protein1, protein2 FROM skempi_mutations"
    )
    for row in cur:
        pid = row[0]
        pdb_proteins.setdefault(pid, set())
        if row[1]:
            pdb_proteins[pid].add(row[1])
        if row[2]:
            pdb_proteins[pid].add(row[2])

    inserted = 0
    residue_batch = []

    for pdb_path in pdb_files:
        pdb_id = pdb_path.stem
        mapping_path = pdb_path.with_suffix(".mapping")

        sequences = parse_sequence_from_pdb(pdb_path)
        mapping = parse_mapping(mapping_path) if mapping_path.exists() else {}

        if not sequences:
            continue

        chains = list(sequences.keys())
        n_residues = sum(len(s) for s in sequences.values())
        proteins = list(pdb_proteins.get(pdb_id, set()))

        # Count mutations in SKEMPI
        n_muts = conn.execute(
            "SELECT COUNT(*) FROM skempi_mutations WHERE pdb_id=?", (pdb_id,)
        ).fetchone()[0]

        conn.execute(
            """
            INSERT OR REPLACE INTO pdb_structures
              (pdb_id, chains, sequences, n_residues, n_mutations, proteins)
            VALUES (?,?,?,?,?,?)
        """,
            (
                pdb_id,
                json.dumps(chains),
                json.dumps(sequences),
                n_residues,
                n_muts,
                json.dumps(proteins),
            ),
        )

        # Build residue table from mapping
        for chain, residues in mapping.items():
            for resname, ch, pdb_num, seq_num in residues:
                aa = THREE_TO_ONE.get(resname, "X")
                residue_batch.append((pdb_id, ch, pdb_num, seq_num, aa))

        inserted += 1
        if inserted % 50 == 0:
            print(f"  Parsed {inserted}/{len(pdb_files)} structures...")

    # Batch insert residues
    conn.executemany(
        """
        INSERT INTO pdb_residues (pdb_id, chain, pdb_resnum, seq_num, residue_aa)
        VALUES (?,?,?,?,?)
    """,
        residue_batch,
    )
    conn.commit()
    print(f"  Inserted {inserted} PDB structures, {len(residue_batch)} residue records")


# ── mark interface residues using SKEMPI location field ───────────────────────
def mark_interface_residues(conn):
    """Use SKEMPI location (COR/RIM/INT/SUR) to flag interface residues."""
    print("Marking interface residues from SKEMPI location data...")
    cur = conn.execute("""
        SELECT DISTINCT pdb_id, mutation_clean, location
        FROM skempi_mutations
        WHERE location IN ('COR','RIM','INT','SUP')
    """)
    interface_set = set()
    for row in cur:
        pdb_id, mut_clean, loc = row
        if loc in ("COR", "RIM"):
            # Parse mutation string e.g. "LA38G" → chain=L, pos=38
            for mut in mut_clean.split(","):
                mut = mut.strip()
                if len(mut) >= 3:
                    chain = mut[1]
                    pos_str = "".join(c for c in mut[2:-1] if c.isdigit())
                    if pos_str:
                        interface_set.add((pdb_id, chain, int(pos_str)))

    batch = list(interface_set)
    conn.executemany(
        """
        UPDATE pdb_residues SET is_interface=1
        WHERE pdb_id=? AND chain=? AND seq_num=?
    """,
        batch,
    )
    conn.commit()
    print(f"  Marked {len(batch)} interface residue positions")


# ── statistics ────────────────────────────────────────────────────────────────
def print_stats(conn):
    stats = {
        "skempi_total": conn.execute(
            "SELECT COUNT(*) FROM skempi_mutations"
        ).fetchone()[0],
        "with_ddG": conn.execute(
            "SELECT COUNT(*) FROM skempi_mutations WHERE ddG_kcal IS NOT NULL"
        ).fetchone()[0],
        "with_koff": conn.execute(
            "SELECT COUNT(*) FROM skempi_mutations WHERE koff_ratio IS NOT NULL"
        ).fetchone()[0],
        "pdb_structures": conn.execute(
            "SELECT COUNT(*) FROM pdb_structures"
        ).fetchone()[0],
        "residues": conn.execute("SELECT COUNT(*) FROM pdb_residues").fetchone()[0],
        "interface_res": conn.execute(
            "SELECT COUNT(*) FROM pdb_residues WHERE is_interface=1"
        ).fetchone()[0],
        "hotspots": conn.execute(
            "SELECT COUNT(*) FROM skempi_mutations WHERE ddG_kcal > 1.0"
        ).fetchone()[0],
    }
    print("\n=== Database Statistics ===")
    for k, v in stats.items():
        print(f"  {k:20s}: {v:,}")
    return stats


# ── main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    if DB_PATH.exists():
        DB_PATH.unlink()
        print(f"Removed old database")

    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA cache_size=10000")

    create_schema(conn)
    load_skempi(conn)
    load_pdb_structures(conn)
    mark_interface_residues(conn)
    stats = print_stats(conn)
    conn.close()

    print(
        f"\n✅ Database built at {DB_PATH} ({DB_PATH.stat().st_size / 1024 / 1024:.1f} MB)"
    )
