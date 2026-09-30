#!/usr/bin/env bash
# Runs once, right after the codespace is created.
set -euo pipefail

echo ""
echo "Setting up the Bioinformatics Teaching Lab..."
echo ""

python data/generate_mock_data.py --outdir data/mock
bash scripts/check_env.sh || true

cat <<'MSG'

Setup complete.

Next steps:
  1. Open README.md and start at Lesson 1
  2. Or jump straight in:  nextflow run main.nf

MSG
