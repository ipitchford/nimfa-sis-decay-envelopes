#!/usr/bin/env bash
# Run every checker (each exits with status 1 on any failure, also under python -O), then rebuild the PDF.
# Usage: uv run --python 3.11 --with-requirements requirements.txt bash build.sh
set -euo pipefail
cd -- "$(dirname -- "$0")"
export MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
python3 -I -O -B exact_checks.py --out exact_checks.json          # seconds
python3 -I -O -B numerical_checks.py --out numerical_checks.json  # seconds
python3 -I -O -B eq20_check.py --out eq20_check.json              # about 2 minutes
python3 -I -O -B illustrations.py --outdir .                      # about 5 minutes
bash build_pdf.sh
