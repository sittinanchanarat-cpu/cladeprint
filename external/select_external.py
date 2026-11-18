"""External test set: published-labelled isolates never used in training.

From labels/literature_labels.tsv, drop isolates with any run already in kmers/samples.tsv, then take
  - 12 each of VNIa, VNIb, VNIc, VNBI, VNBII (pure lineages: specificity + external sensitivity);
  - every hybrid or diploid label (VNI/VNB, VNII/VNB, VNB/VNIV, VNII/VNII, VNI/VNI, VNB/VNB, VNI/VG):
    tests whether the scores name both parents.
One run per isolate: paired Illumina with >= 1.2M spots (enough for the 1M-pair subsample used in
training), smallest download first; otherwise the largest run. Pure-lineage isolates are taken
smallest download first within each country (round-robin over countries). Order of the sheet: pure
lineages interleaved (13 GB, the core specificity test), then hybrids smallest first, so a partial
overnight run loses only the largest hybrid runs. The next 12 per pure lineage are listed as
reserves (set = reserve), used by topup.py when a run fails or turns out low-coverage.

Writes external/samples.tsv (isolate, lineage, source, country, run, spots, size_MB, set).
Usage: micromamba run -n cladeprint python external/select_external.py
"""
import os
import numpy as np
import pandas as pd

os.chdir(os.path.dirname(os.path.abspath(__file__)))
PURE = ["VNIa", "VNIb", "VNIc", "VNBI", "VNBII"]
MIXED = ["VNII/VNII", "VNII/VNB", "VNI/VG", "VNB/VB", "VNB/VNB", "VNI/VNI", "VNI/VNB", "VNB/VNIV"]
N_PURE = 12
SEED = 20261013
MIN_SPOTS = 1_200_000

lab = pd.read_csv("../labels/literature_labels.tsv", sep="\t")
used = set(pd.read_csv("../kmers/samples.tsv", sep="\t").run)
lab = lab[lab.runs.notna()]
lab = lab[~lab.runs.str.split(";").map(lambda r: any(x in used for x in r))]
lab["country"] = lab.country.fillna("").str.strip()
lab = lab.drop_duplicates("isolate")

def read_runinfo(path):
    """The runinfo files are concatenated pages, each with its own header and column count."""
    import csv
    recs, head = [], None
    with open(path, newline="") as fh:
        for row in csv.reader(fh):
            if row and row[0] == "Run":
                head = row
            elif row and head and len(row) == len(head):
                recs.append(dict(zip(head, row)))
    d = pd.DataFrame(recs)[["Run", "spots", "size_MB", "LibraryLayout", "Platform"]]
    return d.assign(spots=pd.to_numeric(d.spots, "coerce").fillna(0), size_MB=pd.to_numeric(d.size_MB, "coerce"))


ri = pd.concat([read_runinfo(f"../census/{f}") for f in os.listdir("../census")
                if f.startswith("runinfo_")]).drop_duplicates("Run").set_index("Run")


def pick_run(runs):
    r = ri.reindex([x for x in runs.split(";") if x in ri.index])
    r = r[(r.Platform == "ILLUMINA") & (r.LibraryLayout == "PAIRED")]
    if r.empty:
        return None
    big = r[r.spots >= MIN_SPOTS]
    best = big.sort_values("size_MB").iloc[0] if len(big) else r.sort_values("spots").iloc[-1]
    return best.name, int(best.spots), float(best.size_MB)


rows = []
seen_runs = set()
rng = np.random.default_rng(SEED)
for cls in MIXED + PURE:
    g = lab[lab.lineage == cls]
    g = g.iloc[rng.permutation(len(g))]
    picks = [pick_run(r) for r in g.runs]
    g = g.assign(run=[p and p[0] for p in picks], spots=[p and p[1] for p in picks],
                 size_MB=[p and p[2] for p in picks]).dropna(subset=["run"])
    # the label table can list one isolate under two names sharing a run: one row per run, overall
    g = g[~g.run.isin(seen_runs)].drop_duplicates("run")
    if cls in PURE:  # smallest downloads first, round-robin over countries so one can't fill the quota
        g = g.sort_values("size_MB", kind="stable")
        g = g.assign(k=g.groupby("country").cumcount()).sort_values(["k", "size_MB"], kind="stable")
        g = g.head(2 * N_PURE).assign(set=["pure"] * min(N_PURE, len(g)) + ["reserve"] * max(0, min(N_PURE, len(g) - N_PURE)))
    else:
        g = g.assign(set="hybrid")
    for _, r in g.iterrows():
        rows.append(dict(isolate=r.isolate, lineage=cls, source=r.source, country=r.country,
                         run=r.run, spots=int(r.spots), size_MB=r.size_MB, set=r.set))
        seen_runs.add(r.run)

out = pd.DataFrame(rows)
mixed = out[out.set == "hybrid"].sort_values("size_MB")
pure = out[out.set != "hybrid"]
pure = pure.assign(k=pure.groupby(["set", "lineage"]).cumcount()).sort_values(["set", "k", "lineage"]).drop(columns="k")
out = pd.concat([pure[pure.set == "pure"], mixed, pure[pure.set == "reserve"]])
out.to_csv("samples.tsv", sep="\t", index=False)
print(out.groupby(["set", "lineage"]).size().to_string())
main = out[out.set != "reserve"]
print(f"{len(main)} isolates (+{(out.set == 'reserve').sum()} reserves), {main.size_MB.sum() / 1024:.1f} GB to download; "
      f"{(main.spots < MIN_SPOTS).sum()} below {MIN_SPOTS:,} spots")
print(main.groupby("lineage").country.agg(lambda c: dict(c.value_counts())).to_string())
