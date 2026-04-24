"""
Flask REST API for the AI-Physics Drug Optimization Platform.

Endpoints:
  GET  /api/stats              - database statistics
  GET  /api/history            - past analysis sessions
  POST /api/analyze            - run full optimization pipeline
  GET  /api/session/<id>       - fetch saved session results
  GET  /api/skempi/search      - search SKEMPI records by PDB/protein
  GET  /api/pdb/<pdb_id>       - get PDB structure info
"""

import sys, json, sqlite3, traceback, threading
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS

from engine.analyzer import ProteinOptimizer
from engine import smallmol
from engine import uniprot

DB_PATH      = Path(__file__).parent.parent / "db" / "skempi.db"
FRONTEND_DIR = Path(__file__).parent.parent / "frontend"
PAI1_DIR     = Path(__file__).parent.parent / "dataset" / "peptide_pai1"
PAI1_POOLS   = {"69-80", "76-88", "69-88", "merged"}

app = Flask(__name__, static_folder=str(FRONTEND_DIR))
CORS(app)

# Single shared optimizer instance, initialized eagerly at module import so the
# first HTTP request doesn't pay the ~1 s cold-start cost. A lock guards any
# lazy re-init path to prevent concurrent workers from duplicating the load.
_optimizer = None
_optimizer_lock = threading.Lock()

def get_optimizer():
    global _optimizer
    if _optimizer is not None:
        return _optimizer
    with _optimizer_lock:
        if _optimizer is None:
            print("[API] Initialising ProteinOptimizer...")
            _optimizer = ProteinOptimizer()
            print("[API] Ready.")
    return _optimizer


# ── health / stats ────────────────────────────────────────────────────────────
@app.route("/api/stats")
def api_stats():
    try:
        opt = get_optimizer()
        stats = opt.get_db_stats()
        stats["iedb_binders"]     = opt.iedb.n_binder
        stats["iedb_nonbinders"]  = opt.iedb.n_nonbinder
        stats["iedb_source"]      = opt.iedb.source
        return jsonify({"ok": True, "stats": stats})
    except Exception as e:
        return jsonify({"ok": False, "error": str(e)}), 500


@app.route("/api/immuno/score", methods=["POST"])
def api_immuno_score():
    """IEDB-PSSM MHC-II T-cell epitope risk for a single sequence."""
    try:
        data = request.get_json(force=True)
        seq = (data or {}).get("sequence", "").strip().upper()
        seq = "".join(a for a in seq if a.isalpha())
        if len(seq) < 9:
            return jsonify({"ok": False, "error": "Sequence must be ≥9 residues"}), 400
        opt = get_optimizer()
        out = {"sequence_len": len(seq)}
        r = opt.iedb.score_sequence(seq)
        if r:
            out["pssm_risk"]        = r["risk"]
            out["pssm_top_hits"]    = r["top_hits"]
            out["pssm_n_hot"]       = r["n_hot_windows"]
            out["pssm_n_windows"]   = r["n_windows"]
            out["pssm_max_score"]   = r["max_score"]
            out["pssm_source"]      = opt.iedb.source
        return jsonify({"ok": True, "result": _clean(out)})
    except Exception as e:
        traceback.print_exc()
        return jsonify({"ok": False, "error": str(e)}), 500


# ── session history ───────────────────────────────────────────────────────────
@app.route("/api/history")
def api_history():
    try:
        opt   = get_optimizer()
        limit = int(request.args.get("limit", 20))
        hist  = opt.get_session_history(limit)
        for h in hist:
            try:
                h["result_summary"] = json.loads(h.get("result_summary") or "{}")
            except: pass
        return jsonify({"ok": True, "sessions": hist})
    except Exception as e:
        return jsonify({"ok": False, "error": str(e)}), 500


# ── main analysis endpoint ────────────────────────────────────────────────────
@app.route("/api/analyze", methods=["POST"])
def api_analyze():
    try:
        data = request.get_json(force=True)
        if not data or "sequence" not in data:
            return jsonify({"ok": False, "error": "Missing 'sequence' field"}), 400

        sequence    = data["sequence"].strip()
        name        = data.get("name", "Query Protein")[:80]
        target      = data.get("target", "Unknown Target")[:80]
        strategy    = data.get("strategy", "mixed")
        max_variants = min(int(data.get("max_variants", 80)), 200)

        if strategy not in ("conservative", "aggressive", "mixed"):
            strategy = "mixed"

        opt    = get_optimizer()
        result = opt.run_optimization(sequence, name=name, target=target,
                                      max_variants=max_variants, strategy=strategy)

        if "error" in result:
            return jsonify({"ok": False, "error": result["error"]}), 400

        # Serialise numpy types
        result = _clean(result)
        return jsonify({"ok": True, "result": result})

    except Exception as e:
        traceback.print_exc()
        return jsonify({"ok": False, "error": str(e)}), 500


# ── fetch saved session ───────────────────────────────────────────────────────
@app.route("/api/session/<session_id>")
def api_session(session_id):
    try:
        conn = sqlite3.connect(DB_PATH)
        conn.row_factory = sqlite3.Row

        sess = conn.execute(
            "SELECT * FROM analysis_sessions WHERE session_id=?", (session_id,)
        ).fetchone()
        if not sess:
            return jsonify({"ok": False, "error": "Session not found"}), 404

        cands = conn.execute("""
            SELECT * FROM candidates WHERE session_id=? ORDER BY rank
        """, (session_id,)).fetchall()

        result = dict(sess)
        result["candidates"] = [dict(c) for c in cands]
        for c in result["candidates"]:
            try: c["mutations"] = json.loads(c.get("mutations") or "[]")
            except: c["mutations"] = []
        try: result["result_summary"] = json.loads(result.get("result_summary") or "{}")
        except: pass
        conn.close()
        return jsonify({"ok": True, "session": result})
    except Exception as e:
        return jsonify({"ok": False, "error": str(e)}), 500


# ── SKEMPI search ─────────────────────────────────────────────────────────────
@app.route("/api/skempi/search")
def api_skempi_search():
    try:
        pdb      = request.args.get("pdb", "").upper()
        protein  = request.args.get("protein", "")
        mut_like = request.args.get("mutation", "")
        limit    = min(int(request.args.get("limit", 50)), 200)

        conn = sqlite3.connect(DB_PATH)
        conn.row_factory = sqlite3.Row

        clauses, params = [], []
        if pdb:
            clauses.append("pdb_id = ?"); params.append(pdb)
        if protein:
            clauses.append("(protein1 LIKE ? OR protein2 LIKE ?)")
            params += [f"%{protein}%", f"%{protein}%"]
        if mut_like:
            clauses.append("mutation_clean LIKE ?"); params.append(f"%{mut_like}%")
        clauses.append("ddG_kcal IS NOT NULL")

        where = " AND ".join(clauses)
        rows = conn.execute(f"""
            SELECT pdb_id, mutation_clean, ddG_kcal, koff_ratio,
                   protein1, protein2, temperature, location
            FROM skempi_mutations WHERE {where}
            ORDER BY ABS(ddG_kcal) DESC LIMIT ?
        """, params + [limit]).fetchall()
        conn.close()

        return jsonify({"ok": True, "records": [dict(r) for r in rows], "total": len(rows)})
    except Exception as e:
        return jsonify({"ok": False, "error": str(e)}), 500


# ── PDB info ──────────────────────────────────────────────────────────────────
@app.route("/api/pdb/<pdb_id>")
def api_pdb(pdb_id):
    try:
        pdb_id = pdb_id.upper()
        conn = sqlite3.connect(DB_PATH)
        conn.row_factory = sqlite3.Row

        struct = conn.execute(
            "SELECT * FROM pdb_structures WHERE pdb_id=?", (pdb_id,)
        ).fetchone()
        if not struct:
            return jsonify({"ok": False, "error": f"PDB {pdb_id} not found"}), 404

        muts = conn.execute("""
            SELECT mutation_clean, ddG_kcal, koff_ratio, location
            FROM skempi_mutations
            WHERE pdb_id=? AND ddG_kcal IS NOT NULL
            ORDER BY ABS(ddG_kcal) DESC LIMIT 20
        """, (pdb_id,)).fetchall()

        result = dict(struct)
        try: result["sequences"] = json.loads(result.get("sequences") or "{}")
        except: pass
        try: result["chains"] = json.loads(result.get("chains") or "[]")
        except: pass
        try: result["proteins"] = json.loads(result.get("proteins") or "[]")
        except: pass
        result["top_mutations"] = [dict(m) for m in muts]
        conn.close()
        return jsonify({"ok": True, "pdb": result})
    except Exception as e:
        return jsonify({"ok": False, "error": str(e)}), 500


# ── SKEMPI ddG distribution (for charts) ─────────────────────────────────────
@app.route("/api/skempi/distribution")
def api_skempi_dist():
    try:
        conn = sqlite3.connect(DB_PATH)
        rows = conn.execute(
            "SELECT ddG_kcal, koff_ratio, location FROM skempi_mutations WHERE ddG_kcal IS NOT NULL"
        ).fetchall()
        ddGs     = [r[0] for r in rows]
        koffs    = [r[1] for r in rows if r[1] is not None]
        locs     = {}
        for r in rows:
            if r[2]: locs[r[2]] = locs.get(r[2], 0) + 1
        conn.close()
        return jsonify({
            "ok":    True,
            "ddG":   ddGs,
            "koffs": koffs,
            "locs":  locs,
            "total": len(ddGs)
        })
    except Exception as e:
        return jsonify({"ok": False, "error": str(e)}), 500


# ── UniProt lookup ────────────────────────────────────────────────────────────
@app.route("/api/uniprot/<acc>")
def api_uniprot(acc):
    """Fetch a UniProt entry (sequence + features + PDB/AF xrefs).

    Cached on disk under db/uniprot_cache/. Appends ``?refresh=1`` to bypass.
    """
    try:
        use_cache = request.args.get("refresh", "0") != "1"
        entry = uniprot.fetch_by_id(acc, use_cache=use_cache)
        return jsonify({"ok": True, "entry": _clean(entry)})
    except uniprot.UniProtError as e:
        return jsonify({"ok": False, "error": str(e)}), 404
    except Exception as e:
        traceback.print_exc()
        return jsonify({"ok": False, "error": str(e)}), 500


@app.route("/api/uniprot/demo")
def api_uniprot_demo():
    """Return the curated demo-target list shown in the frontend dropdown."""
    return jsonify({"ok": True, "targets": uniprot.DEMO_TARGETS})


# ── small-molecule (KRAS G12D POC) ────────────────────────────────────────────
@app.route("/api/smallmol/poc")
def api_smallmol_poc():
    """Serve cached small-molecule screening results for the KRAS G12D POC."""
    try:
        return jsonify({"ok": True, "result": _clean(smallmol.load_poc_results())})
    except Exception as e:
        traceback.print_exc()
        return jsonify({"ok": False, "error": str(e)}), 500


@app.route("/api/smallmol/plot")
def api_smallmol_plot():
    """Serve the Pareto plot PNG."""
    p = smallmol.poc_plot_path()
    if not p.exists():
        return jsonify({"ok": False, "error": "plot not found"}), 404
    return send_from_directory(str(p.parent), p.name)


# ── PAI-1 peptide demo (Year-1 baseline: 3 seeds + candidate pools) ──────────
def _pai1_read_json(rel_path: str):
    p = PAI1_DIR / rel_path
    if not p.exists():
        return None
    with p.open("r", encoding="utf-8") as fh:
        return json.load(fh)


@app.route("/api/peptide_pai1/seeds")
def api_pai1_seeds():
    """List the three PAI-1 mimicking seed peptides from Peptide info.zip."""
    data = _pai1_read_json("peptides.json")
    if data is None:
        return jsonify({"ok": False, "error": "peptides.json not found"}), 404
    return jsonify({"ok": True, "data": data})


@app.route("/api/peptide_pai1/summary")
def api_pai1_summary():
    """Aggregate pool statistics (runtime, J distribution, penalty counts)."""
    data = _pai1_read_json("candidates/summary.json")
    if data is None:
        return jsonify({"ok": False, "error": "summary.json not found"}), 404
    return jsonify({"ok": True, "data": data})


@app.route("/api/peptide_pai1/candidates/<pool>")
def api_pai1_candidates(pool):
    """Return the variant pool (full JSON) for one of: 69-80, 76-88, 69-88, merged."""
    if pool not in PAI1_POOLS:
        return jsonify({"ok": False, "error": f"Unknown pool '{pool}'. "
                        f"Valid: {sorted(PAI1_POOLS)}"}), 400
    data = _pai1_read_json(f"candidates/{pool}/variants.json")
    if data is None:
        return jsonify({"ok": False, "error": f"{pool}/variants.json not found"}), 404

    # Allow the client to cap how many variants come across the wire.
    try:
        limit = int(request.args.get("limit", 0))
    except ValueError:
        limit = 0
    if limit > 0 and isinstance(data.get("variants"), list):
        data = dict(data)
        data["variants"] = data["variants"][:limit]
        data["_truncated_to"] = limit
    return jsonify({"ok": True, "data": _clean(data)})


@app.route("/api/peptide_pai1/benchmark")
def api_pai1_benchmark():
    """Alanine-scan consistency benchmark vs Stefansson 2004."""
    data = _pai1_read_json("benchmark/ala_consistency.json")
    if data is None:
        return jsonify({"ok": False, "error": "benchmark not found"}), 404
    return jsonify({"ok": True, "data": data})


@app.route("/api/peptide_pai1/modifiable_sites/<name>")
def api_pai1_sites(name):
    """Per-residue modifiable-site list (allowed mutations + SKEMPI evidence)."""
    if name not in {"69-80", "76-88", "69-88"}:
        return jsonify({"ok": False, "error": f"Unknown seed '{name}'"}), 400
    data = _pai1_read_json(f"modifiable_sites/{name}.json")
    if data is None:
        return jsonify({"ok": False, "error": f"sites for {name} not found"}), 404
    return jsonify({"ok": True, "data": _clean(data)})


# ── serve frontend ────────────────────────────────────────────────────────────
@app.route("/")
@app.route("/<path:path>")
def serve_frontend(path="index.html"):
    try:
        return send_from_directory(str(FRONTEND_DIR), path)
    except:
        return send_from_directory(str(FRONTEND_DIR), "index.html")


# ── helpers ───────────────────────────────────────────────────────────────────
def _clean(obj):
    """Recursively convert numpy scalars to Python primitives."""
    import numpy as np
    if isinstance(obj, dict):
        return {k: _clean(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_clean(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


if __name__ == "__main__":
    import os
    port = int(os.environ.get("PORT", 7860))
    # Pre-warm models before binding the socket so the first HTTP request is fast
    get_optimizer()
    print(f"[API] Starting server on port {port}...")
    app.run(host="0.0.0.0", port=port, debug=False, threaded=True)
