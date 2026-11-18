"""Summarise the leave-one-clade-out runs (R2; census/holdouts.py).

For each holdout, over the held-out runs that pass QC in the v2 run (so they would have trained):
  - ancestor levels (a group above the clade, still trained): score for that group and whether
    classify calls it; compared with the same runs' leave-one-out score in v2 (clade present).
  - the clade's own level: its group is gone (or, for a VNIa complex / deep species clade, trained
    without it). Right answers: the trained group for a complex or deep clade; NO_CALL for a whole
    held-out lineage or sublineage. A sibling call is an error.
Literature-labelled members are reported separately from those added by propagation or nearest
neighbour.

Usage: python census/r2_summary.py census/out [prefix]   (prefix "" = cladeprint, "exact_" = baseline)
"""

import sys
from pathlib import Path

import pandas as pd

LEVELS = ["species", "lineage", "sublineage"]


def label_path(r, level):
    cols = LEVELS[: LEVELS.index(level) + 1]
    v = [r[c] for c in cols]
    return None if any(pd.isna(x) for x in v) else "/".join(v)


def clade_level(tag):
    if tag in ("VNII", "VNB"):
        return "lineage"
    if tag.startswith("VNIa_C") or tag.startswith("other"):
        return None  # within a still-trained group
    return "sublineage"


def main(out, prefix=""):
    out = Path(out)
    qc = pd.read_csv(out / "qc.tsv", sep="\t", dtype=str).set_index("sample")
    loo = pd.read_csv(out / "loo.tsv", sep="\t").set_index(["sample", "group"]).sensitivity
    rows = []
    for d in sorted((out / "r2").iterdir()):
        tag = d.name
        if prefix and not tag.startswith(prefix) or not prefix and tag.startswith("exact_"):
            continue
        tag0 = tag[len(prefix):]
        calls = d / "heldout_calls.tsv"
        if tag0 == "full" or not calls.exists():
            continue
        h = pd.read_csv(d / "heldout.tsv", sep="\t", dtype=str).set_index("sample")
        c = pd.read_csv(calls, sep="\t", dtype=str).set_index("sample")
        h = h[qc.reindex(h.index).train == "true"]
        h["src"] = "literature" if tag0.startswith("other") else h.how.str.split().str[0]
        cl = clade_level(tag0)
        for src, hs in h.groupby(h.src.where(h.src == "literature", "added")):
            for level in LEVELS:
                truth = hs.apply(lambda r: label_path(r, level), axis=1)
                if truth.isna().all():
                    continue
                got = c.loc[hs.index, level]
                ok = truth.notna()
                if cl is not None and LEVELS.index(level) >= LEVELS.index(cl):
                    # the clade's own level (or below): its group is not trained
                    nocall = (got == "NA") | got.isna()
                    rows.append(dict(holdout=tag0, src=src, level=level, role="held-out level",
                                     n=int(ok.sum()), no_call=int(nocall[ok].sum()),
                                     wrong=int((~nocall)[ok].sum()), right=None,
                                     score_med=None, loo_med=None, score_p05=None))
                    continue
                # trained ancestor (or, for a complex/deep clade, the group itself)
                sc = pd.Series({s: float(c.at[s, "score:" + t]) if "score:" + t in c.columns else float("nan")
                                for s, t in truth[ok].items()})
                lo = pd.Series({s: loo.get((s, t), float("nan")) for s, t in truth[ok].items()})
                rows.append(dict(holdout=tag0, src=src, level=level, role="trained",
                                 n=int(ok.sum()), right=int((got[ok] == truth[ok]).sum()),
                                 wrong=int(((got[ok] != truth[ok]) & (got[ok] != "NA")).sum()),
                                 no_call=int((got[ok] == "NA").sum()),
                                 score_med=sc.median(), score_p05=sc.quantile(0.05),
                                 loo_med=lo.median()))
    r = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    print(r.round(3).to_string(index=False))
    r.to_csv(out / f"r2/summary{('_' + prefix.rstrip('_')) if prefix else ''}.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main(*sys.argv[1:])
