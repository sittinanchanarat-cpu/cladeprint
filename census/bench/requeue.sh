#!/usr/bin/env bash
# after the SKA2 OOMs (18 Oct): finish exact_other1 and the last Kover CV targets, then the solo benchmark
( bash ~/r2_job.sh exact_other1 ) &
( bash ~/r4_kover_cv.sh ) &
wait
echo "$(date -u) REQUEUE DONE" >> ~/r4.log
bash ~/solo_bench.sh
