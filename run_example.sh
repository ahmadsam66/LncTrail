#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python "${SCRIPT_DIR}/lnctrail.py" \
  --lncrna "${SCRIPT_DIR}/examples/lncrna.fa" \
  --targets "${SCRIPT_DIR}/examples/target_transcripts.fa" \
  --outdir "${SCRIPT_DIR}/demo_output" \
  --engine auto \
  --threads 4 \
  --batch-size 1
