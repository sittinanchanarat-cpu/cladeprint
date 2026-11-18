"""Summarise the R5 published-assay audit.

Input: census/out/r5_sites.tsv (examples/primer_sites output: best binding site of each oligo in
each census genome), labels/assays/r5_oligos_map.tsv, labels/assays/assay_list_for_review_audited.tsv,
census/out/prop/samples.tsv (labels) and census/out/qc.tsv.

An oligo is taken to bind when its best site has at most MAX_MM mismatches and none of them in the
3' tail primer_sites reports (last 5 bases): internal mismatches are usually tolerated, a 3' one
is not. The median mismatch count is reported alongside, since many older assays sit at 2-3
mismatches from every modern genome. An assay is taken to
work on a genome when every one of its oligos binds. Allele-discriminating assays (MAMA, allele-
specific) are reported per oligo instead, since their call compares two primers rather than asking
whether one binds.

Usage: python census/r5_summary.py [out_dir]
Writes out_dir/r5_assay_by_group.tsv and prints the headline table.
"""

import sys
from pathlib import Path

import pandas as pd

MAX_MM = 3  # primer_sites searches to 3; 'none' means >= 4 and is treated as no binding
TAIL_STRICT = 2
LEVELS = ["species", "lineage", "sublineage"]


def main(out="census/out"):
    out = Path(out)
    sites = pd.read_csv(out / "r5_sites.tsv", sep="\t", dtype={"mismatches": str, "mm_3prime": str})
    amap = pd.read_csv("labels/assays/r5_oligos_map.tsv", sep="\t")
    alist = pd.read_csv("labels/assays/assay_list_for_review_audited.tsv", sep="\t", dtype=str).fillna("")
    s = pd.read_csv(out / "prop/samples.tsv", sep="\t", dtype=str)
    qc = pd.read_csv(out / "qc.tsv", sep="\t", dtype=str).set_index("sample")
    s = s[qc.reindex(s["sample"]).train.values == "true"].set_index("sample")

    sites["mm"] = pd.to_numeric(sites.mismatches, errors="coerce")
    sites["tailmm"] = pd.to_numeric(sites.mm_3prime, errors="coerce")
    # PCR tolerates internal mismatches but not one at the 3' end, which blocks extension:
    # bind = a site within MAX_MM mismatches with none in the last 5 bases (primer_sites' tail).
    sites["binds"] = sites.mm.le(MAX_MM) & sites.tailmm.eq(0)
    sites["perfect"] = sites.mm.eq(0)

    # map the deduplicated oligo rows back to every assay oligo that shares the sequence
    seq_of = dict(zip(amap.oligo, amap.sequence))
    first_for_seq = {}
    for o, q in seq_of.items():
        first_for_seq.setdefault(q, o)
    rep = {o: first_for_seq[q] for o, q in seq_of.items()}
    by = sites.set_index(["oligo", "sample"]).binds

    meta = alist.set_index("id")
    amap["id"] = amap.oligo.str.split("|").str[0]
    amap["assay"] = amap.id.map(meta.assay)
    amap["role"] = amap.id.map(meta.role).fillna("")
    amap["notes"] = amap.id.map(meta.notes).fillna("")
    amap["discriminating"] = (
        amap.assay.str.contains("MAMA", case=False, na=False)
        | amap.role.str.contains("allele-specific", case=False, na=False)
    )

    def path(row, level):
        cols = LEVELS[: LEVELS.index(level) + 1]
        v = [row.get(c) for c in cols]
        return None if any(pd.isna(x) or x == "" for x in v) else "/".join(v)

    groups = {}
    for lv in LEVELS:
        p = s.apply(lambda r: path(r, lv), axis=1)
        for g in sorted(p.dropna().unique()):
            m = set(p.index[p == g])
            if len(m) >= 5:
                groups[g] = m
    samples = list(s.index)

    rows = []
    for assay, grp in amap.groupby("assay"):
        if grp.discriminating.any():
            continue  # reported per oligo elsewhere
        oligos = [rep[o] for o in grp.oligo if rep[o] in by.index.get_level_values(0)]
        if not oligos:
            continue
        works = None
        for o in oligos:
            v = by.loc[o].reindex(samples).fillna(False)
            works = v if works is None else (works & v)
        for g, mem in groups.items():
            inm = works.reindex(list(mem)).mean()
            outm = works.reindex([x for x in samples if x not in mem]).mean()
            mm_med = float(pd.concat([sites[sites.oligo == o].set_index("sample").mm.reindex(list(mem))
                                      for o in oligos]).median())
            rows.append(dict(assay=assay, n_oligos=len(oligos), group=g, n_members=len(mem),
                             member_binds=round(float(inm), 4), nonmember_binds=round(float(outm), 4),
                             member_median_mm=mm_med))
    d = pd.DataFrame(rows)
    d.to_csv(out / "r5_assay_by_group.tsv", sep="\t", index=False)
    print(f"{d.assay.nunique()} non-discriminating assays x {d.group.nunique()} groups "
          f"-> {out}/r5_assay_by_group.tsv")
    top = d[d.member_binds >= 0.95].sort_values(["assay", "nonmember_binds"])
    print("\nassay x group with >= 95% of members bound (lowest off-target first):")
    print(top.head(40).to_string(index=False))


if __name__ == "__main__":
    main(*sys.argv[1:])
