import sys, traceback, threading, uuid
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))

from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS

from engine.analyzer import (
    ProteinOptimizer,
    ToolAFailure,
    ToolBFailure,
    PipelineTimeout,
)
from engine import uniprot

app = Flask(__name__)
CORS(app)

jobs = {}
jobs_lock = threading.Lock()

optimizer = ProteinOptimizer()


def background_worker(job_id, sequence, tool_a, tool_b):
    def progress_cb(pct):
        with jobs_lock:
            jobs[job_id]["progress"] = int(pct)

    try:
        with jobs_lock:
            jobs[job_id]["status"] = "processing"
            jobs[job_id]["progress"] = 0

        results = optimizer.run_pipeline(
            sequence,
            tool_a=tool_a,
            tool_b=tool_b,
            progress_cb=progress_cb,
        )

        with jobs_lock:
            jobs[job_id]["status"] = "completed"
            jobs[job_id]["progress"] = 100
            jobs[job_id]["results"] = results
            jobs[job_id]["completed_at"] = datetime.now().isoformat()

    except (ToolAFailure, ToolBFailure, PipelineTimeout, Exception) as e:
        if isinstance(e, ToolAFailure):
            kind = "tool_a"
        elif isinstance(e, ToolBFailure):
            kind = "tool_b"
        elif isinstance(e, PipelineTimeout):
            kind = "timeout"
        else:
            kind = "internal"
        with jobs_lock:
            jobs[job_id]["status"] = "error"
            jobs[job_id]["error"] = str(e)
            jobs[job_id]["error_kind"] = kind
        traceback.print_exc()


@app.route("/api/jobs", methods=["POST"])
def create_job():
    data = request.get_json()
    sequence = data.get("sequence")
    tool_a = data.get("tool_a", "boltz2")
    tool_b = data.get("tool_b", "mmgbsa")

    if not sequence:
        return jsonify({"error": "Sequence is required"}), 400

    job_id = str(uuid.uuid4())
    with jobs_lock:
        jobs[job_id] = {
            "id": job_id,
            "status": "pending",
            "progress": 0,
            "created_at": datetime.now().isoformat(),
            "params": {"sequence": sequence, "tool_a": tool_a, "tool_b": tool_b},
        }

    thread = threading.Thread(target=background_worker, args=(job_id, sequence, tool_a, tool_b))
    thread.start()

    return jsonify({"job_id": job_id}), 202


@app.route("/api/jobs/<job_id>", methods=["GET"])
def get_job(job_id):
    with jobs_lock:
        job = jobs.get(job_id)
    if not job:
        return jsonify({"error": "Job not found"}), 404
    return jsonify(job)


@app.route("/api/uniprot/demo")
def api_uniprot_demo():
    return jsonify({"ok": True, "targets": uniprot.DEMO_TARGETS})


@app.route("/", defaults={"path": ""})
@app.route("/<path:path>")
def serve(path):
    frontend_dist = Path(__file__).parent.parent / "frontend" / "dist"
    if path != "" and (frontend_dist / path).exists():
        return send_from_directory(str(frontend_dist), path)
    return send_from_directory(str(frontend_dist), "index.html")


if __name__ == "__main__":
    import os
    app.run(host="0.0.0.0", port=int(os.environ.get("PORT", 7860)), debug=False)
