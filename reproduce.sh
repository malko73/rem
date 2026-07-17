#!/usr/bin/env bash
set -euo pipefail
python src/rem4_numerical.py \
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
