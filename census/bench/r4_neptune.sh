#!/usr/bin/env bash
# R4: Neptune per target (defaults, 4 processes: 8 x 3.5 GB starved the R2 runs); log ~/r4.log
export PATH=$HOME/.local/bin:$PATH MAMBA_ROOT_PREFIX=$HOME/micromamba
eval "$(micromamba shell hook -s bash)"; micromamba activate bench
cd ~/cladeprint/census/out/bench
tail -n +2 targets.tsv | cut -f2 | while read t; do
  [ -s neptune/out_$t/consolidated/consolidated.fasta ] && continue
  /usr/bin/time -v -o neptune/$t.time neptune -i $(cat neptune/$t.inc) -e $(cat neptune/$t.exc) \
    -o neptune/out_$t -p 4 > neptune/$t.log 2>&1
  rc=$?
  echo "$(date -u) neptune $t exit $rc $(grep -E 'Elapsed|Maximum resident' neptune/$t.time | tr -s ' ' | tr '\n' ' ')" >> ~/r4.log
  du -sh neptune/out_$t >> ~/r4.log
  find neptune/out_$t -mindepth 1 -maxdepth 1 -type d ! -name consolidated -exec rm -rf {} +
done
echo "$(date -u) NEPTUNE DONE" >> ~/r4.log
