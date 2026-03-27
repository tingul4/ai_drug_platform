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

import sys, json, sqlite3, traceback
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS

from engine.analyzer import ProteinOptimizer

DB_PATH      = Path(__file__).parent.parent / "db" / "skempi.db"
FRONTEND_DIR = Path(__file__).parent.parent / "frontend"

app = Flask(__name__, static_folder=str(FRONTEND_DIR))
CORS(app)

# Single shared optimizer instance (loads SKEMPI stats once)
_optimizer = None

def get_optimizer():
    global _optimizer
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
        return jsonify({"ok": True, "stats": stats})
    except Exception as e:
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
    print(f"[API] Starting server on port {port}...")
    app.run(host="0.0.0.0", port=port, debug=False, threaded=True)
