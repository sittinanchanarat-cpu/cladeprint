#!/usr/bin/env bash
# R4 Neptune fairness test (18 Oct): one target (VNII) alone with -p 2 after Neptune -p 8 used the
# whole 64 GB at signature extraction and hung the instance. Memory capped at 40 GB in a cgroup
# (an overrun kills Neptune, not the machine) and 6 h wall time. usage: neptune_p2_test.sh <pid to wait for>
tail --pid=$1 -f /dev/null
B=~/cladeprint/census/out/bench; S=$B/solo; t=C.neoformans_VNII
cd $B && rm -rf neptune/out_$t
s=$(date -u +%s)
sudo systemd-run --scope --quiet -p MemoryMax=40G -p MemorySwapMax=0 \
  sudo -u ubuntu bash -c "export PATH=\$HOME/micromamba/envs/bench/bin:\$PATH; cd $B; \
  timeout 6h /usr/bin/time -v -o $S/neptune_p2_$t.time neptune -i \$(cat neptune/$t.inc) -e \$(cat neptune/$t.exc) \
  -o neptune/out_$t -p 2 > neptune/${t}_p2.log 2>&1"
rc=$?
echo "$(date -u) solo neptune_p2 $t exit $rc wall $(( $(date -u +%s) - s )) s $(grep -E 'Elapsed|Maximum resident' $S/neptune_p2_$t.time 2>/dev/null | tr -s ' ' | tr '\n' ' ') signatures $(grep -c '>' neptune/out_$t/consolidated/consolidated.fasta 2>/dev/null)" >> ~/r4.log
du -sh neptune/out_$t >> ~/r4.log
find neptune/out_$t -mindepth 1 -maxdepth 1 ! -name consolidated -exec rm -rf {} +
