"""Second-round lineage labels for runs the literature does not label, from a first-round bundle.

The literature labels only 37% of the census C. neoformans runs (VNII: 19 training runs), which
leaves the lineage benchmark thin and most runs usable only at species level. A run with no label
is ignored for every lineage target, so it contributes nothing there.

This labels such a run with the lineage whose markers it carries, when the call is unambiguous:
top score >= MIN_TOP and the next lineage <= MAX_SECOND. Sublineages are labelled the same way,
within the lineage just assigned. A run that fails either rule keeps an empty label and stays
ignored, so uncertainty costs nothing.

This is self-training: the labels come from markers built on the literature-labelled runs, so a
bundle built from them cannot be validated by the same scores. Two checks that do not use the
markers are printed (and must be read before adopting the output):
  - recovery: every literature-labelled training run's top score must be its own label. These are
    leave-one-out scores for the run's own groups, so a run is not scored against itself.
  - genotype clusters (census/out/clusters/clusters.tsv, split-k-mer SNP distances - independent of
    the markers): a propagated label must agree with the literature labels of its cluster-mates.
Report both in the paper, keep the literature-only bundle as the primary result, and validate the
propagated one by leave-one-clade-out.

Usage: micromamba run -n cladeprint python census/propagate_labels.py [OUT] [SCORES_DIR]
  OUT         census output directory (default out): reads OUT/samples.tsv and OUT/clusters/
  SCORES_DIR  directory holding scores.tsv of the first-round run (default OUT)
Writes OUT/samples_propagated.tsv (+ a `label_source` column) and prints the two checks.
"""
import os
import sys

import pandas as pd

os.chdir(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1] if len(sys.argv) > 1 else "out"
SCORES = sys.argv[2] if len(sys.argv) > 2 else OUT
MIN_TOP, MAX_SECOND = 0.85, 0.10
SPECIES = "C.neoformans"
LINEAGES = ["VNI", "VNII", "VNB"]
SUBLINEAGES = {"VNI": ["VNIa", "VNIb", "VNIc"], "VNB": ["VNBI", "VNBII"]}

sheet = pd.read_csv(f"{OUT}/samples.tsv", sep="\t", dtype=str).fillna("")
scores = pd.read_csv(f"{SCORES}/scores.tsv", sep="\t", dtype={"sample": str})
clusters = pd.read_csv(f"{OUT}/clusters/clusters.tsv", sep="\t", dtype=str).set_index("sample")
d = sheet.merge(scores, on="sample", suffixes=("", "_sc"))
# both tables have a `train` column: the sheet's override (auto/yes/no) and the run's QC decision
d["is_train"] = d["train_sc"].astype(str).str.lower() == "true"


def call(row, groups, prefix):
    """(label, top, second) for one run over `groups`, or ('', ...) when no column exists."""
    cols = [f"{prefix}/{g}" for g in groups]
    vals = sorted(
        ((pd.to_numeric(row.get(c), errors="coerce"), g) for c, g in zip(cols, groups)),
        key=lambda t: (-(t[0] if pd.notna(t[0]) else -1), t[1]),
    )
    if not vals or pd.isna(vals[0][0]):
        return "", float("nan"), float("nan")
    top, second = vals[0][0], (vals[1][0] if len(vals) > 1 and pd.notna(vals[1][0]) else 0.0)
    return (vals[0][1] if top >= MIN_TOP and second <= MAX_SECOND else ""), top, second


rows = []
for _, r in d.iterrows():
    lin, lin_top, lin_2 = call(r, LINEAGES, SPECIES)
    sub = ""
    if r.species == SPECIES:
        parent = r.lineage or lin
        if parent in SUBLINEAGES:
            sub, _, _ = call(r, SUBLINEAGES[parent], f"{SPECIES}/{parent}")
    # only runs the first round trained on are labelled: the rest failed QC (low coverage,
    # mixture, foreign), where a confident-looking call would be misleading
    keep = r.is_train and r.species == SPECIES
    rows.append((lin if keep else "", lin_top, lin_2, sub if keep else ""))
d["p_lineage"], d["p_top"], d["p_second"], d["p_sublineage"] = zip(*rows)

lit = d[d.is_train & (d.species == SPECIES) & (d.lineage != "")]
best = [
    max(((pd.to_numeric(r.get(f"{SPECIES}/{g}"), errors="coerce"), g) for g in LINEAGES))[1]
    for _, r in lit.iterrows()
]
agree = sum(a == b for a, b in zip(best, lit.lineage))
print(f"recovery: top lineage == literature label for {agree}/{len(lit)} labelled training runs")
print(f"  of those the rule calls confidently: {int((lit.p_lineage == lit.lineage).sum())}")
if agree < len(lit):
    for b, (_, r) in zip(best, lit.iterrows()):
        if b != r.lineage:
            print(f"  MISMATCH {r['sample']}: literature {r.lineage}, top {b}")

# cluster cross-check (independent of the markers)
d["cluster"] = d["sample"].map(clusters["cluster_5"]) if "cluster_5" in clusters.columns else d.get("cluster")
tr = d[d.is_train]
ok = bad = lonely = 0
conflicts = []
for cl, grp in tr.groupby("cluster"):
    labs = set(grp.loc[grp.lineage != "", "lineage"])
    new = grp[(grp.lineage == "") & (grp.p_lineage != "")]
    if not len(new):
        continue
    if not labs:
        lonely += len(new)
        continue
    for _, r in new.iterrows():
        if labs == {r.p_lineage}:
            ok += 1
        else:
            bad += 1
            conflicts.append((cl, r["sample"], r.p_lineage, sorted(labs)))
print(f"cluster cross-check: {ok} agree, {bad} disagree, {lonely} in clusters with no labelled member")
for c in conflicts[:10]:
    print("  conflict:", c)

out = sheet.copy()
src = []
for i, r in d.iterrows():
    if not r.is_train:
        src.append("")
    elif r.lineage:
        src.append("literature")
    elif r.p_lineage:
        out.loc[i, "lineage"] = r.p_lineage
        if r.p_sublineage:
            out.loc[i, "sublineage"] = r.p_sublineage
        src.append("propagated")
    else:
        src.append("")
    if r.lineage and not r.sublineage and r.p_sublineage:
        out.loc[i, "sublineage"] = r.p_sublineage
out["label_source"] = src
out.to_csv(f"{OUT}/samples_propagated.tsv", sep="\t", index=False)
before = sheet[sheet.species == SPECIES].groupby(["lineage", "sublineage"]).size()
after = out[out.species == SPECIES].groupby(["lineage", "sublineage"]).size()
print("\nC. neoformans labels, before -> after:")
print(pd.DataFrame({"before": before, "after": after}).fillna(0).astype(int).to_string())
print(f"\nwrote {OUT}/samples_propagated.tsv ({int((out.label_source == 'propagated').sum())} propagated)")
