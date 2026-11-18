"""R3 focused tree taxa (~90): bootstrap support where the R3 claims need it, at a size SKA2 can
merge (the 530-taxon merge needed > 60 GB).

10 literature-labelled runs per Ashton 2019 VNIa sub-clade (from the R3 taxon set, seeded draw);
every run that the nearest-neighbour placement left ambiguous or outside the named sub-clades,
ERR842663 and the census candidates; 5 VNIb and 5 VNIc (plus the candidates labelled with them);
3 VNB to root (the alignment has no distance limit, unlike the split-k-mer pair table).

Usage: python census/r3_focus.py census/out   -> census/out/r3/focus_taxa.tsv
"""

import sys
from pathlib import Path

import pandas as pd

PER_SUB, PER_SIB, N_ROOT, SEED = 10, 5, 3, 20261018


def main(out):
    out = Path(out)
    taxa = pd.read_csv(out / "r3/taxa.tsv", sep="\t", dtype=str).set_index("sample")
    L = pd.read_csv("labels/literature_labels.tsv", sep="\t", dtype=str)
    L["sub"] = L.label_raw.str.extract(r"(VNIa-[0-9A-Z]+)", expand=False)
    lab = L[L["sub"].notna()].assign(run=L.runs.str.split(";")).explode("run").set_index("run")["sub"]
    lab = lab[~lab.index.duplicated()]
    nn = pd.read_csv(out / "r3/vnia_nn_subclade.tsv", sep="\t", index_col=0)
    picks = {}
    for sub in sorted(lab.unique()):
        pool = [s for s in taxa.index if lab.get(s) == sub and s != "ERR842663"]
        for s in pd.Series(pool).sample(min(PER_SUB, len(pool)), random_state=SEED):
            picks[s] = sub
    odd = list(nn.index[(nn.d1 > 300) | (nn.d2 < 3 * nn.d1)]) + ["ERR842663"]
    odd += list(taxa.index[taxa.role.str.contains("candidate|small", regex=True)])
    for s in odd:
        if s in taxa.index:
            picks.setdefault(s, "odd: " + str(taxa.role.get(s)))
    for sib in ["VNIb", "VNIc"]:
        pool = [s for s in taxa.index[taxa.role == sib]]
        for s in pd.Series(pool).sample(PER_SIB, random_state=SEED):
            picks.setdefault(s, sib)
    for s in taxa.index[taxa.role == "VNB root"][:N_ROOT]:
        picks[s] = "VNB root"
    t = pd.DataFrame({"sample": list(picks), "group": list(picks.values())})
    t["ashton"] = t["sample"].map(lab)
    t["label_prop"] = t["sample"].map(taxa.label_prop)
    t.to_csv(out / "r3/focus_taxa.tsv", sep="\t", index=False)
    print(t.group.str.replace(r"odd: .*", "odd", regex=True).value_counts().to_string())
    print("total", len(t))


if __name__ == "__main__":
    main(sys.argv[1])
