"""Evaluate `cladeprint classify` calls on the external test set against the published labels.

Input: a classify TSV (--long) for the external runs, made with a bundle built from the training
isolates (e.g. kmers/tmp/runtest/out/bundle). Pure isolates: is the lineage call right, and the
sublineage call for VNIa/b/c? Hybrids: is one parent called and the other flagged MIXED?
Low-coverage isolates (external/stats.tsv low_cov) are reported separately, as in score_external.py.
Usage: micromamba run -n cladeprint python external/classify_eval.py CALLS.tsv
"""
import os, sys
import pandas as pd

os.chdir(os.path.dirname(os.path.abspath(__file__)))
calls = pd.read_csv(sys.argv[1], sep="\t").set_index("sample")
calls.columns = ["call_" + c for c in calls.columns]
lab = pd.read_csv("samples.tsv", sep="\t").drop_duplicates("run").set_index("run")
st = pd.read_csv("stats.tsv", sep="\t").drop_duplicates("run", keep="last").set_index("run")
d = calls.join(lab[["isolate", "lineage"]]).join(st[["low_cov", "peak_cov"]])
SUB = ("VNIa", "VNIb", "VNIc")
major = lambda l: "VNI" if l in SUB else l
leaf = lambda x: x.split("/")[-1] if isinstance(x, str) else "NA"

pure = d[~d.lineage.str.contains("/")].copy()
pure["lin_ok"] = pure.call_lineage == pure.lineage.map(major)
pure["sub_ok"] = [leaf(c) == l if l in SUB else True for c, l in zip(pure.call_sublineage, pure.lineage)]
for name, part in (("usable coverage", pure[pure.low_cov == 0]), ("low coverage", pure[pure.low_cov == 1])):
    print(f"Pure isolates, {name}: {len(part)}; lineage right {part.lin_ok.sum()}, "
          f"sublineage right {part.sub_ok.sum()} (of {part.lineage.isin(SUB).sum()} VNIa/b/c + others n/a); "
          f"confidence {part.call_lineage_confidence.value_counts().to_dict()}")
bad = pure[~(pure.lin_ok & pure.sub_ok)]
if len(bad):
    print(bad[["isolate", "lineage", "call_lineage", "call_lineage_score", "call_sublineage", "call_flags", "low_cov"]].to_string())

hyb = d[d.lineage.str.contains("/")]
print(f"\nHybrids and diploids: {len(hyb)}")
cols = ["isolate", "lineage", "call_lineage", "call_lineage_score", "call_lineage_confidence", "call_sublineage", "call_flags"]
print(hyb[cols].sort_values("lineage").to_string())

# UNEXPLAINED (bundle format 2): the --kmc-db databases carry no stats, so low coverage comes from
# external/stats.tsv here; classify does not flag low-coverage samples
if "call_unexplained" in d:
    import tomllib
    flag = tomllib.load(open(os.path.join(os.path.dirname(sys.argv[1]), "out/bundle/bundle.toml"), "rb"))["unexplained_flag"]
    d["flagged"] = (d.call_unexplained > flag) & (d.low_cov == 0)
    print(f"\nUNEXPLAINED threshold {flag:.3f}")
    for name, part in (("pure, usable coverage", pure.index), ("hybrids", hyb.index)):
        x = d.loc[part]
        print(f"{name}: {len(x)}, unexplained median {x.call_unexplained.median():.3f}, max {x.call_unexplained.max():.3f}, flagged {x.flagged.sum()}")
    print(d.loc[hyb.index][["isolate", "lineage", "call_unexplained", "flagged", "low_cov"]].sort_values("lineage").to_string())
