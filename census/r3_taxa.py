"""R3: taxon set for the VNIa phylogeny (VNIa note, org-roam; census/out/vnia_lit_components_600.tsv).

One run per cluster at 10 SNPs/Mb within VNIa (literature or propagated label), preferring a
literature-labelled run and then the highest genome-peak coverage; every literature VNIa run in the
small complexes (C4-C6 and the singletons, at 600 SNPs/Mb) so their placement is seen directly;
N_SIB runs each of VNIb and VNIc (one per 10-SNP cluster, seeded draw) and N_ROOT VNB runs to root;
and the runs flagged as mislabel or intermediate candidates by the census checks.
Only runs that pass QC in the v2 run.

Usage: python census/r3_taxa.py census/out   -> census/out/r3/taxa.tsv
"""

import sys
from pathlib import Path

import pandas as pd

N_SIB, N_ROOT, SEED = 30, 6, 20261018
CANDIDATES = ["SRR822918", "SRR836604", "SRR798180", "SRR1944858", "SRR642856", "ERR2624237",
              "SRR1063253"]


def main(out):
    out = Path(out)
    lit = pd.read_csv(out / "samples.tsv", sep="\t", dtype=str).set_index("sample")
    prop = pd.read_csv(out / "prop/samples.tsv", sep="\t", dtype=str).set_index("sample")
    qc = pd.read_csv(out / "qc.tsv", sep="\t", dtype=str).set_index("sample")
    c10 = pd.read_csv(out / "clusters_10_50/clusters.tsv", sep="\t", dtype=str).set_index("sample")
    comp = pd.read_csv(out / "vnia_lit_components_600.tsv", sep="\t", index_col=0).iloc[:, 0]
    d = prop[(qc.reindex(prop.index).train == "true") & (prop.species == "C.neoformans")].copy()
    d["c10"] = c10.reindex(d.index).cluster_10
    d["peak"] = pd.to_numeric(qc.reindex(d.index).peak)
    d["lit"] = lit.reindex(d.index).sublineage
    d["is_lit"] = d.lit.notna()
    picks = {}

    v = d[d.sublineage == "VNIa"].sort_values(["is_lit", "peak"], ascending=False)
    for s in v.drop_duplicates("c10").index:
        picks[s] = "VNIa rep"
    for s, c in comp.items():
        if c not in ("C1", "C2", "C3") and s in d.index:
            picks.setdefault(s, f"VNIa {c}")
    for sub, n in [("VNIb", N_SIB), ("VNIc", N_SIB)]:
        x = d[d.sublineage == sub].sample(frac=1, random_state=SEED).sort_values("is_lit", ascending=False, kind="stable")
        for s in x.drop_duplicates("c10").head(n).index:
            picks.setdefault(s, sub)
    r = d[(d.lineage == "VNB") & d.is_lit].sample(frac=1, random_state=SEED).drop_duplicates("c10").head(N_ROOT)
    for s in r.index:
        picks.setdefault(s, "VNB root")
    for s in CANDIDATES:
        picks[s] = picks.get(s, "") + " candidate"
    t = pd.DataFrame({"sample": list(picks), "role": [v.strip() for v in picks.values()]})
    t["label_lit"] = t["sample"].map(lit.sublineage)
    t["label_prop"] = t["sample"].map(prop.sublineage)
    t["complex600"] = t["sample"].map(comp)
    t["country"] = t["sample"].map(lit.country)
    (out / "r3").mkdir(exist_ok=True)
    t.to_csv(out / "r3/taxa.tsv", sep="\t", index=False)
    print(t.role.value_counts().to_string())
    print("total", len(t))


if __name__ == "__main__":
    main(sys.argv[1])
