#!/usr/bin/env bash
set -euo pipefail

echo "NOTE: reproduce.sh runs the legacy REM4 asymmetric-XY benchmark, not the current Spec v2.2 Phase D validation chain." >&2
echo "For current reproduction guidance, see REPRODUCING.md." >&2

# Use .venv/bin/python if available (preferred), fallback to python3
if [ -x .venv/bin/python ]; then
    PYTHON=".venv/bin/python"
elif command -v python3 &>/dev/null; then
    PYTHON="python3"
else
    PYTHON="python"
fi

exec "$PYTHON" src/rem4_numerical.py \
  --outdir outputs \
  --mode all \
  --state-mode ground \
  --j12 1.5 \
  --j23 0.6 \
  --h 0.2 \
  --lambda 0.2 \
  --lambda-max 1.0 \
  --grid-r 80 \
  --grid-l 80
