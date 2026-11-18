#!/usr/bin/env bash
# Extra markers for exact leave-one-out sensitivity: k-mers held by exactly need-1 target isolates
# (need = ceil(0.9 n)), passed through the same out-group filters as the real markers.
#
# Holding out a target isolate X leaves the out-group unchanged and lowers the threshold to
# ceil(0.9 (n-1)) = need-1 (true for every target here). So the markers built without X are the
# current markers plus those extras that X lacks. Outputs groups/T.x.ed2.u64.npy; filter stats go
# to loo_extra.tsv (kept out of near_miss.tsv / edit_miss.tsv).
#
# Usage: micromamba run -n cladeprint bash loo_extra.sh
set -euo pipefail
cd "$(dirname "$0")"
declare -A N=([VNII]=37 [VNBI]=14 [VNBII]=8 [VNI]=49 [VNIa]=24 [VNIb]=12 [VNIc]=13)
targets=(VNIc VNIa VNIb VNBII VNBI VNI VNII)

for t in "${targets[@]}"; do
  need=$(python -c "import math; print(math.ceil(0.9 * ${N[$t]}))")
  x=groups/$t.x
  if [[ ! -s $x.strict.kmc_suf ]]; then
    kmc_tools -t2 transform groups/$t.in -ci$((need - 1)) -cx$((need - 1)) reduce $x.cand >/dev/null 2>&1
    kmc_tools -t2 simple $x.cand groups/$t.out kmers_subtract $x.strict >/dev/null 2>&1
    rm -f $x.cand.kmc_pre $x.cand.kmc_suf
  fi
  ln -sf $t.out.u64.npy $x.out.u64.npy  # same out-group, so reuse its cached k-mer array
done

for f in near_miss.tsv edit_miss.tsv; do cp $f $f.bak; done
xs=("${targets[@]/%/.x}")
python near_miss.py "${xs[@]}"
python edit_miss.py "${xs[@]}"
python - <<'EOF'
import pandas as pd
nm = pd.read_csv("near_miss.tsv", sep="\t").set_index("target")
ed = pd.read_csv("edit_miss.tsv", sep="\t").set_index("target")
x = [t for t in nm.index if t.endswith(".x")]
res = nm.loc[x, ["strict", "survive"]].rename(columns={"survive": "hamming_survive"})
res["edit_survive"] = ed.loc[x, "survive"]
res.to_csv("loo_extra.tsv", sep="\t")
print(res.to_string())
EOF
for f in near_miss.tsv edit_miss.tsv; do mv $f.bak $f; done
