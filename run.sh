#!/bin/bash
# ================================================================
# AIM3 胜肽藥物優化平台  |  Startup Script
# ================================================================
# Usage:  ./run.sh          (default port 7860)
#         ./run.sh 8080     (custom port)
# ================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PORT="${1:-7860}"
PY="$SCRIPT_DIR/third_party/gmx_env/bin/python"

if [ ! -x "$PY" ]; then
    echo "[ERROR] gmx_env not found at $PY"
    echo "        See README.md §環境建置 to bootstrap third_party/gmx_env."
    exit 1
fi

# Rebuild frontend if dist is missing or older than any source file
DIST="$SCRIPT_DIR/frontend/dist/index.html"
NEEDS_BUILD=0
if [ ! -f "$DIST" ]; then
    NEEDS_BUILD=1
elif [ -n "$(find "$SCRIPT_DIR/frontend/src" -newer "$DIST" -print -quit 2>/dev/null)" ]; then
    NEEDS_BUILD=1
fi
if [ "$NEEDS_BUILD" = "1" ]; then
    echo "[INIT] Frontend dist is stale; running npm run build..."
    (cd "$SCRIPT_DIR/frontend" && npm run build)
fi

# Kill any previous server still holding the port
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
PORT=$PORT "$PY" api/server.py
