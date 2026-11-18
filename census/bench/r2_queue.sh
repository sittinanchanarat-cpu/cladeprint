#!/usr/bin/env bash
main="VNII VNB VNBI VNBII VNIb VNIc VNIa VNIa_C1 VNIa_C2 VNIa_C3 other1 other2"
exact="exact_full $(for t in $main; do printf 'exact_%s ' $t; done)"
echo "$(date -u) QUEUE $main $exact" >> ~/r2.log
printf '%s\n' $main $exact | xargs -P 2 -I{} bash ~/r2_job.sh {}
echo "$(date -u) R2 QUEUE DONE" >> ~/r2.log
