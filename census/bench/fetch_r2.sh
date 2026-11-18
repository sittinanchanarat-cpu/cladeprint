#!/usr/bin/env bash
# copy finished R2 holdouts home (small tables only)
cd ~/Projects/cladeprint
for t in $(ssh -i ~/.ssh/cladeprint.pem ubuntu@3.236.75.1 'grep " END .* exit 0" ~/r2.log | awk "{print \$(NF-2)}"'); do
  mkdir -p census/out/r2/$t
  [ -s census/out/r2/$t/heldout_calls.tsv ] && continue
  for f in heldout.tsv heldout_calls.tsv qc.tsv loo.tsv; do scp -q -i ~/.ssh/cladeprint.pem ubuntu@3.236.75.1:cladeprint/census/out/r2/$t/$f census/out/r2/$t/ 2>/dev/null; done
  scp -q -i ~/.ssh/cladeprint.pem ubuntu@3.236.75.1:cladeprint/census/out/r2/$t/bundle/groups.tsv census/out/r2/$t/groups.tsv
done
