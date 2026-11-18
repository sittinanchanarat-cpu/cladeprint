#!/usr/bin/env bash
# Overnight driver: pure lineages, top up low-coverage/failed ones from reserves (up to 3 rounds),
# then hybrids (smallest first), then score. Every step skips finished runs, so it can be rerun.
cd "$(dirname "$0")"
SET=pure bash run_external.sh
for round in 1 2 3; do
  python topup.py || break
  SHEET=topup.tsv bash run_external.sh
done
python topup.py
SET=hybrid bash run_external.sh
python score_external.py > scores_summary.txt 2>&1
echo "overnight finished $(date)"
