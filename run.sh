#!/bin/bash
# ================================================================
# AI-Physics 混合策略全維度藥物優化平台  |  Startup Script
# ================================================================
# Usage:  ./run.sh          (default port 7860)
#         ./run.sh 8080     (custom port)
# ================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PORT="${1:-7860}"

# Check Python
if ! command -v python3 &>/dev/null; then
    echo "[ERROR] python3 not found. Please install Python 3.9+."
    exit 1
fi

# Check / install dependencies
echo "[INIT] Checking dependencies..."
python3 -c "import flask, flask_cors, Bio, numpy, scipy" 2>/dev/null || {
    echo "[INIT] Installing required packages..."
    pip3 install flask flask-cors biopython numpy scipy --break-system-packages -q
}

# Verify database
DB="$SCRIPT_DIR/db/skempi.db"
if [ ! -f "$DB" ]; then
    echo "[ERROR] Database not found: $DB"
    echo "        Please rebuild with:  python3 db/build_database.py"
    exit 1
fi

echo "[INIT] Database OK: $(du -sh "$DB" | cut -f1)"

# Kill any previous server still holding the port (e.g. a failed Ctrl+C)
if command -v lsof &>/dev/null; then
    STALE=$(lsof -ti tcp:$PORT 2>/dev/null || true)
elif command -v ss &>/dev/null; then
    STALE=$(ss -tlnp 2>/dev/null | awk -v p=":$PORT" '$4 ~ p {print}' | grep -oP 'pid=\K[0-9]+' | head -1)
fi
if [ -n "$STALE" ]; then
    echo "[INIT] Killing stale process on port $PORT (PID $STALE)..."
    kill -TERM "$STALE" 2>/dev/null || true
    sleep 1
    kill -KILL "$STALE" 2>/dev/null || true
fi

echo "[INIT] Starting server on http://localhost:$PORT ..."
echo "[INFO] Press Ctrl+C to stop."
echo ""

cd "$SCRIPT_DIR"
PORT=$PORT python3 api/server.py
