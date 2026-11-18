"""Leave-one-clade-out (R2) and exact-baseline (R4) run directories for the census.

Each holdout is a full `cladeprint run` on the literature-label sheet (v2, the primary bundle)
with every member of one clade set to train = no; the held-out runs are then classified against
the bundle built without them.

Membership has to include runs that carry no literature label for the clade, or they leak it back
into training at the parent level (48 unlabelled C. neoformans runs are VNII, and still train the
species markers). A run is a member if
  - its literature or propagated label (census/out/prop/samples.tsv) is the clade, or
  - it has no label at that level and its nearest run that has one (SNP distance,
    clusters_neo/pairs.tsv) is in the clade, within NN_MAX SNPs/Mb.
VNIa complexes: single linkage at 600 SNPs/Mb over literature VNIa (the structure in the VNIa note),
extended to other VNIa-labelled or unlabelled runs whose nearest literature VNIa run is in the
complex. Other species: the largest and second-largest cluster at 1,000 SNPs/Mb per species
(clusters_other/), held out together in two runs (one clade per species; they do not interact
except through the out-group, which can only add cross-reacting markers, i.e. it is conservative).

Usage: python census/holdouts.py census/out     (writes census/out/r2/<tag>/ ...)
"""

import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

NN_MAX = 500  # SNPs/Mb: nearest labelled neighbour must be this close to transfer a label
VNIA_LINK = 600
OTHER_CUT = "cluster_1000"

CONFIG = """samples = "samples.tsv"
output = "."
levels = ["species", "lineage", "sublineage"]

[markers]
max_out = 1
near_miss = {near_miss}

[resources]
threads = {threads}
memory_gb = 8
"""


def nearest(pairs, query, labelled):
    """For each run in `query`: (nearest run in `labelled`, distance), from the pair table."""
    p = pairs[
        (pairs.a.isin(query) & pairs.b.isin(labelled))
        | (pairs.b.isin(query) & pairs.a.isin(labelled))
    ]
    q = pd.concat(
        [
            p[p.a.isin(query)].rename(columns={"a": "q", "b": "n"}),
            p[p.b.isin(query)].rename(columns={"b": "q", "a": "n"}),
        ]
    )
    q = q[q.n.isin(labelled)]
    best = q.sort_values("d").drop_duplicates("q").set_index("q")
    return best


def components(pairs, members, cut):
    parent = {s: s for s in members}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    e = pairs[pairs.a.isin(members) & pairs.b.isin(members) & (pairs.d <= cut)]
    for a, b in zip(e.a, e.b):
        parent[find(a)] = find(b)
    comp = defaultdict(list)
    for s in members:
        comp[find(s)].append(s)
    return sorted(comp.values(), key=len, reverse=True)


def main(out):
    out = Path(out)
    lit = pd.read_csv(out / "samples.tsv", sep="\t", dtype=str)
    prop = pd.read_csv(out / "prop/samples.tsv", sep="\t", dtype=str).set_index("sample")
    pairs = pd.read_csv(
        out / "clusters_neo/pairs.tsv", sep="\t", usecols=[0, 1, 5]
    ).set_axis(["a", "b", "d"], axis=1)
    neo = lit[lit.species == "C.neoformans"]
    holdouts = {}  # tag -> (members, how)

    def by_label(level, value, parent_level=None, parent_value=None):
        lab = prop.loc[neo["sample"], level]
        lit_lab = neo.set_index("sample")[level]
        how = {}
        for s in neo["sample"]:
            if lit_lab.get(s) == value:
                how[s] = "literature"
            elif lab.get(s) == value:
                how[s] = "propagated"
        scope = neo["sample"]
        if parent_level:
            scope = scope[prop.loc[scope, parent_level].values == parent_value]
        unl = [s for s in scope if pd.isna(lab.get(s))]
        known = lab.dropna().index
        nn = nearest(pairs, unl, known)
        for s, r in nn.iterrows():
            if r.d <= NN_MAX and lab[r.n] == value:
                how[s] = f"nearest {r.n} {r.d:.0f}"
        return how

    for v in ["VNII", "VNB"]:
        holdouts[v] = by_label("lineage", v)
    for parent, subs in [("VNI", ["VNIa", "VNIb", "VNIc"]), ("VNB", ["VNBI", "VNBII"])]:
        for v in subs:
            holdouts[v] = by_label("sublineage", v, "lineage", parent)

    # VNIa complexes
    lit_vnia = list(neo.loc[neo.sublineage == "VNIa", "sample"])
    comps = components(pairs, lit_vnia, VNIA_LINK)[:3]
    other_vnia = [
        s
        for s in neo["sample"]
        if s not in set(lit_vnia)
        and prop.loc[s, "lineage"] in ("VNI", None)
        and prop.loc[s, "sublineage"] in ("VNIa",)
        or (s not in set(lit_vnia) and pd.isna(prop.loc[s, "sublineage"]))
    ]
    nn = nearest(pairs, other_vnia, lit_vnia)
    lit_sub = neo.set_index("sample").sublineage.dropna()
    nn_any = nearest(pairs, other_vnia, list(lit_sub.index))
    for i, comp in enumerate(comps, 1):
        cs = set(comp)
        how = {s: "literature" for s in comp}
        for s, r in nn.iterrows():
            if r.n in cs and r.d <= VNIA_LINK and lit_sub.get(nn_any.loc[s, "n"]) == "VNIa":
                how[s] = f"nearest {r.n} {r.d:.0f}"
        holdouts[f"VNIa_C{i}"] = how

    # other species: deep clades at 1,000 SNPs/Mb
    co = pd.read_csv(out / "clusters_other/clusters.tsv", sep="\t", dtype=str)
    co = co.merge(lit[["sample", "species"]], on="sample")
    for rank in (1, 2):
        how = {}
        for sp, d in co.groupby("species"):
            if sp == "C.decagattii":
                continue
            top = d[OTHER_CUT].value_counts()
            if len(top) < rank:
                continue
            cl = top.index[rank - 1]
            for s in d.loc[d[OTHER_CUT] == cl, "sample"]:
                how[s] = f"{sp} {cl}"
        holdouts[f"other{rank}"] = how

    r2 = out / "r2"
    summary = []
    for tag, how in holdouts.items():
        for near_miss in (True, False):
            name = tag if near_miss else f"exact_{tag}"
            d = r2 / name
            d.mkdir(parents=True, exist_ok=True)
            sheet = lit.copy()
            sheet.loc[sheet["sample"].isin(how), "train"] = "no"
            sheet.to_csv(d / "samples.tsv", sep="\t", index=False)
            (d / "config.toml").write_text(
                CONFIG.format(near_miss=str(near_miss).lower(), threads=8)
            )
            held = lit[lit["sample"].isin(how)].copy()
            held["how"] = held["sample"].map(how)
            held.to_csv(d / "heldout.tsv", sep="\t", index=False)
        src = pd.Series(how).str.split().str[0].value_counts().to_dict()
        summary.append((tag, len(how), src))
    d = r2 / "exact_full"
    d.mkdir(parents=True, exist_ok=True)
    lit.to_csv(d / "samples.tsv", sep="\t", index=False)
    (d / "config.toml").write_text(CONFIG.format(near_miss="false", threads=8))
    for tag, n, src in summary:
        print(f"{tag:10s} {n:5d} {src}")


if __name__ == "__main__":
    main(sys.argv[1])
