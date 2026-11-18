"""Why cladeprint's VNI sublineage markers score high on the near-miss check (Fig. 3A; R4 subset).

bench_eval.py counts a marker as a near miss when ANY one out-group genome holds it within 2 edits. cladeprint
tolerates one out-group run per marker (max_out = 1; such k-mers also do not block in the near-miss filters),
so a marker held by one out-group genome is allowed by design. Per cladeprint target and split this prints:
  all      share within 2 edits of >= 1 out-group genome (= near2 in eval/summary.tsv)
  no_cand  the same without the VNI mislabel/intermediate candidates of census/r3_fig.py
  ge2      share within 2 edits of >= 2 out-group genomes
and the out-group genomes that hold the most markers. Writes census/out/bench/eval/nearmiss.tsv.
Usage: python census/bench_nearmiss.py
"""
import os

import numpy as np
import pandas as pd

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
EV = "census/out/bench/eval"
CAND = {"SRR822918", "SRR836604", "SRR798180", "SRR1944858", "SRR642856", "ERR2624237", "SRR1063253", "ERR842663"}
LEAVES = {"C.neoformans/VNI/VNIa": ["VNIa"], "C.neoformans/VNI/VNIb": ["VNIb"], "C.neoformans/VNI/VNIc": ["VNIc"],
          "C.neoformans/VNB/VNBI": ["VNBI"], "C.neoformans/VNB/VNBII": ["VNBII"],
          "C.neoformans/VNI": ["VNIa", "VNIb", "VNIc"], "C.neoformans/VNB": ["VNBI", "VNBII"],
          "C.neoformans/VNII": ["VNII"]}

idx = pd.read_csv(f"{EV}/index.tsv", sep="\t")
g = np.array(open(f"{EV}/genomes.txt").read().split())
d = np.fromfile(f"{EV}/dist.u8", dtype=np.uint8).reshape(len(g), -1)
b = pd.read_csv("census/out/bench/samples.tsv", sep="\t", dtype=str).set_index("sample").loc[g]
rows = []
for t, leaf in LEAVES.items():
    r = idx[(idx.tool == "cp") & (idx.target == t)].iloc[0]
    m = d[:, r.start:r.end]
    for split in ("train", "test"):
        sel = (~b.leaf.isin(leaf) & (b.split == split)).values
        near = m[sel] <= 2
        share = pd.Series(near.mean(axis=1), index=g[sel]).sort_values(ascending=False)
        rows.append(dict(target=t, split=split, all=near.any(0).mean(),
                         no_cand=(m[sel & ~np.isin(g, list(CAND))] <= 2).any(0).mean(),
                         ge2=(near.sum(0) >= 2).mean(),
                         top=", ".join(f"{k}({b.loc[k, 'leaf']}) {v:.2f}" for k, v in share.head(2).items())))
out = pd.DataFrame(rows)
out.to_csv(f"{EV}/nearmiss.tsv", sep="\t", index=False, float_format="%.4f")
pd.set_option("display.width", 200)
print(out.round(3).to_string(index=False))
