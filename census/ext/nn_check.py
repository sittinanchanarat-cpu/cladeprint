"""Lineage calls on runs with no lineage label vs the nearest literature-labelled run by SNP distance.

Runs labelled by the markers (label_source = propagated, census/propagate_labels.py) were in no
lineage's training set. Their marker-based lineage/sublineage is compared with the label of the
closest literature-labelled run in census/out/clusters_neo/pairs.tsv (split k-mer SNPs/Mb).
Writes census/ext/nn_check.tsv. Usage: python census/ext/nn_check.py
"""
import os

import pandas as pd

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
s = pd.read_csv("census/out/samples_propagated.tsv", sep="\t", dtype=str).fillna("").set_index("sample")
lit = s[s.label_source == "literature"]
prop = s[s.label_source == "propagated"]
p = pd.read_csv("census/out/clusters_neo/pairs.tsv", sep="\t", usecols=["sample_a", "sample_b", "snps_per_mb"])
both = pd.concat([p.rename(columns={"sample_a": "q", "sample_b": "r"}),
                  p.rename(columns={"sample_b": "q", "sample_a": "r"})])
both = both[both.q.isin(prop.index) & both.r.isin(lit.index)]
nn = both.loc[both.groupby("q").snps_per_mb.idxmin()].set_index("q")
out = prop[["lineage", "sublineage"]].join(nn, how="inner")
out["nn_lineage"] = out.r.map(lit.lineage)
out["nn_sublineage"] = out.r.map(lit.sublineage)
out.index.name = "sample"
out.to_csv("census/ext/nn_check.tsv", sep="\t")
for lv in ("lineage", "sublineage"):
    o = out[(out[lv] != "") & (out[f"nn_{lv}"] != "")]
    print(f"{lv}: {(o[lv] == o[f'nn_{lv}']).sum()}/{len(o)} agree; NN distance median "
          f"{o.snps_per_mb.median():.0f} SNPs/Mb, p90 {o.snps_per_mb.quantile(0.9):.0f}, > 500: {(o.snps_per_mb > 500).sum()}")
print(out.lineage.value_counts().to_dict())
