"""Census-bundle classify of isolates no census group was trained on: 30 hybrid/diploid isolates and
22 untrained pure isolates (21 low coverage at 1M read pairs) from the pilot external set.

census/ext/sheet.tsv + db/ (symlinks to external/db, stats from external/stats.tsv) ->
`cladeprint classify --kmc-db` -> calls.tsv. Labels: external/samples.tsv, except 102-14, whose
VNI/VNB label is corrected to VNII/VNB in the census sheet. Writes census/ext/summary.tsv.
Usage: python census/ext/summarize.py
"""
import os

import pandas as pd

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
c = pd.read_csv("census/ext/calls.tsv", sep="\t", keep_default_na=False)
sh = pd.read_csv("census/ext/sheet.tsv", sep="\t").set_index("sample")
c["isolate"] = c["sample"].map(sh.isolate)
c["label"] = c["sample"].map(sh.lineage)
c.loc[c.isolate == "102-14", "label"] = "VNII/VNB"
c["kind"] = c.label.map(lambda l: "pure" if "/" not in l else
                        ("AD hybrid" if l.endswith("VNIV") else "A x gattii" if l.endswith("VG") else
                         "diploid, one lineage" if len(set(l.split("/"))) == 1 else "lineage hybrid"))
c["mixed"] = c["flags"].str.contains("MIXED")
c["excess"] = c["flags"].str.contains("EXCESS_KMERS")


def parents_found(r):
    a, b = r.label.split("/")
    if r.kind == "AD hybrid":
        return {r.species, *[f.split(":")[-1].split("=")[0] for f in r.flags.split(",") if "species:MIXED" in f]} \
            >= {"C.neoformans", "C.deneoformans"}
    if r.kind == "A x gattii":
        return "MIXED:C.gattii" in r.flags
    called = {r.lineage.split("/")[-1]} | {f.split("/")[-1].split("=")[0] for f in r.flags.split(",") if "MIXED" in f}
    return {a, b} <= called


c["both_parents"] = [parents_found(r) if r.kind not in ("pure", "diploid, one lineage") else None
                     for r in c.itertuples()]


PARENT = {"VNBI": "VNB", "VNBII": "VNB", "VNIa": "VNI", "VNIb": "VNI", "VNIc": "VNI"}


def right(r):
    """Label is a sublineage: right when the deepest call is it or its lineage (sublineage not called)."""
    if r.species == "NA":
        return "no call"
    if r.species != "C.neoformans":
        return "WRONG"
    if r.lineage == "NA":
        return "species only"
    lin, sub = r.lineage.split("/")[-1], r.sublineage.split("/")[-1]
    if lin != PARENT[r.label]:
        return "WRONG"
    return "right" if sub in (r.label, "NA") else "WRONG"


c["outcome"] = [right(r) if r.kind == "pure" else "" for r in c.itertuples()]
c.to_csv("census/ext/summary.tsv", sep="\t", index=False)
for k, g in c.groupby("kind"):
    print(f"{k:22s} n={len(g):2d}  MIXED {g.mixed.sum():2d}  EXCESS {g.excess.sum():2d}  both parents "
          f"{int(g.both_parents.fillna(False).sum())}  low cov {g.low_coverage.astype(str).str.lower().eq('true').sum()}")
p = c[c.kind == "pure"]
print("pure:", p.outcome.value_counts().to_dict(), "| confidence:", p.lineage_confidence.value_counts().to_dict())
print(c[c.kind == "lineage hybrid"][["isolate", "label", "lineage", "flags"]].to_string(index=False))
