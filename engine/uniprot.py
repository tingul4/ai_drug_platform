"""
UniProt REST client with local file cache.

Usage:
    from engine.uniprot import fetch_by_id
    entry = fetch_by_id("P01116")
    #  -> {"id", "primary_accession", "sequence", "length",
    #      "protein_name", "gene", "organism",
    #      "features": [{"type","description","start","end"}, ...],
    #      "pdb_xrefs": ["7RPZ", ...], "alphafold_id": "AF-P01116-F1"}

We pull only the fields that actually feed into the optimizer:
  - canonical sequence  (primary input for analyzer.run_optimization)
  - feature track       (later: gate alanine-scan away from active sites / PTMs)
  - PDB / AlphaFold IDs (later: prefer these over k-mer search)
"""

from pathlib import Path
import json
import time
import requests

UNIPROT_URL   = "https://rest.uniprot.org/uniprotkb/{acc}.json"
CACHE_DIR     = Path(__file__).parent.parent / "db" / "uniprot_cache"
REQUEST_TIMEOUT = 15  # seconds

CACHE_DIR.mkdir(parents=True, exist_ok=True)

# Feature types we surface to the frontend / optimizer.
# UniProt uses these controlled-vocab strings in the JSON payload.
RELEVANT_FEATURES = {
    "Active site", "Binding site", "Site",
    "Domain", "Region", "Motif",
    "Modified residue",          # phosphorylation etc.
    "Glycosylation",
    "Disulfide bond",
    "Signal", "Transit peptide", "Propeptide",
    "Natural variant",
}


class UniProtError(Exception):
    pass


def _cache_path(acc: str) -> Path:
    return CACHE_DIR / f"{acc.upper()}.json"


def _http_get(acc: str) -> dict:
    url = UNIPROT_URL.format(acc=acc)
    r = requests.get(url, timeout=REQUEST_TIMEOUT,
                     headers={"Accept": "application/json"})
    if r.status_code == 404:
        raise UniProtError(f"UniProt accession '{acc}' not found")
    if r.status_code != 200:
        raise UniProtError(f"UniProt HTTP {r.status_code} for '{acc}'")
    return r.json()


def _parse(raw: dict) -> dict:
    """Reduce the UniProt JSON to the fields we actually use."""
    acc = raw.get("primaryAccession", "")
    seq_block = raw.get("sequence", {}) or {}
    sequence  = seq_block.get("value", "")
    length    = seq_block.get("length", len(sequence))

    # Protein name: recommendedName.fullName.value, fall back to submittedNames
    prot = raw.get("proteinDescription", {}) or {}
    rec  = prot.get("recommendedName", {}) or {}
    full = (rec.get("fullName") or {}).get("value")
    if not full and prot.get("submissionNames"):
        full = ((prot["submissionNames"][0] or {}).get("fullName") or {}).get("value")
    protein_name = full or acc

    # Gene name
    genes = raw.get("genes", []) or []
    gene = None
    if genes:
        gene = ((genes[0].get("geneName") or {}).get("value"))

    # Organism
    org = raw.get("organism", {}) or {}
    organism = org.get("scientificName") or ""

    # Features — keep only the types the optimizer cares about
    features = []
    for f in raw.get("features", []) or []:
        ftype = f.get("type")
        if ftype not in RELEVANT_FEATURES:
            continue
        loc = f.get("location", {}) or {}
        start = (loc.get("start") or {}).get("value")
        end   = (loc.get("end")   or {}).get("value")
        features.append({
            "type":        ftype,
            "description": f.get("description", ""),
            "start":       start,
            "end":         end,
        })

    # Cross-references: PDB + AlphaFold
    pdb_xrefs = []
    alphafold_id = None
    for x in raw.get("uniProtKBCrossReferences", []) or []:
        db = x.get("database")
        if db == "PDB":
            pid = x.get("id")
            if pid: pdb_xrefs.append(pid)
        elif db == "AlphaFoldDB":
            alphafold_id = x.get("id")

    return {
        "id":                acc,
        "primary_accession": acc,
        "sequence":          sequence,
        "length":            length,
        "protein_name":      protein_name,
        "gene":              gene,
        "organism":          organism,
        "features":          features,
        "pdb_xrefs":         pdb_xrefs[:25],   # truncate for payload size
        "alphafold_id":      alphafold_id,
        "source":            "UniProtKB",
    }


def fetch_by_id(acc: str, use_cache: bool = True) -> dict:
    """Fetch a UniProt entry by accession (e.g. 'P01116'). Cached on disk."""
    acc = (acc or "").strip().upper()
    if not acc:
        raise UniProtError("Empty accession")
    # Minimal sanity check — UniProt accessions are [OPQ][0-9][A-Z0-9]{3}[0-9]
    # or [A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}. We just enforce alnum-only.
    if not acc.isalnum() or len(acc) < 6 or len(acc) > 10:
        raise UniProtError(f"Malformed accession: {acc}")

    cp = _cache_path(acc)
    if use_cache and cp.exists():
        try:
            return json.loads(cp.read_text())
        except Exception:
            pass  # fall through to refetch

    raw    = _http_get(acc)
    parsed = _parse(raw)
    parsed["_fetched_at"] = int(time.time())
    cp.write_text(json.dumps(parsed, ensure_ascii=False))
    return parsed


# A curated demo list used by the frontend and the demo doc.
# Picked because each has strong SKEMPI coverage OR is a headline drug target.
DEMO_TARGETS = [
    {"acc": "P01241", "label": "Human Growth Hormone (hGH)",
     "note": "251 SKEMPI records (hGH / hGHbp)"},
    {"acc": "P62593", "label": "TEM-1 β-lactamase",
     "note": "277 SKEMPI records vs BLIP; antibiotic resistance"},
    {"acc": "P01563", "label": "Interferon α-2",
     "note": "220 SKEMPI records; cytokine-receptor interface"},
    {"acc": "P08246", "label": "Human leukocyte elastase",
     "note": "259 SKEMPI records; serine-protease / inhibitor"},
    {"acc": "P01116", "label": "KRAS (G12D POC target)",
     "note": "matches small-molecule POC (PDB 7RPZ)"},
]
