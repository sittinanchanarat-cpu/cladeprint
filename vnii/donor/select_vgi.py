"""Donor search: a country-stratified sample of public C. gattii (VGI) read sets.

From the census manifest, independent natural isolates filed as C. gattii, one run each (Illumina,
>= 1.2M spots; the smallest such run, as in external/select_external.py), at most MAX_MB to
download. Paired runs preferred; single-end allowed (all Malawi runs, PRJNA1103328, are
single-end). Quota per country, weighted to Africa, where VNII is common; random within a country. Species labels in NCBI are unreliable (pilot/README.md), so fetch_vgi.sh checks each
isolate's genome-wide distance from WM276 before it counts as VGI.

Writes vgi_samples.tsv (biosample, strain, country, run, spots, size_MB, layout).
Usage: micromamba run -n cladeprint python vnii/donor/select_vgi.py
"""
import csv, os
import numpy as np
import pandas as pd

os.chdir(os.path.dirname(os.path.abspath(__file__)))
SEED = 20261014
MIN_SPOTS = 1_200_000
MAX_MB = 1500
QUOTA = {"Malawi": 6, "Zambia": 3, "South Africa": 3, "USA": 3, "China": 3, "Colombia": 2, "Mexico": 2,
         "Brazil": 2, "Australia": 2, "India": 2, "Canada": 1, "New Zealand": 1, "Papua New Guinea": 1,
         "Venezuela": 1, "Paraguay": 1, "Urugay": 1}


def read_runinfo(path):
    """The runinfo files are concatenated pages, each with its own header and column count."""
    recs, head = [], None
    with open(path, newline="") as fh:
        for row in csv.reader(fh):
            if row and row[0] == "Run":
                head = row
            elif row and head and len(row) == len(head):
                recs.append(dict(zip(head, row)))
    d = pd.DataFrame(recs)[["Run", "spots", "size_MB", "LibraryLayout", "Platform"]]
    return d.assign(spots=pd.to_numeric(d.spots, "coerce").fillna(0), size_MB=pd.to_numeric(d.size_MB, "coerce"))


ri = pd.concat([read_runinfo(f"../../census/{f}") for f in os.listdir("../../census")
                if f.startswith("runinfo_")]).drop_duplicates("Run").set_index("Run")


def pick_run(runs):
    r = ri.reindex([x for x in runs.split(";") if x in ri.index])
    r = r[(r.Platform == "ILLUMINA") & (r.spots >= MIN_SPOTS) & (r.size_MB <= MAX_MB)]
    if r.empty:
        return None
    best = r.assign(single=r.LibraryLayout != "PAIRED").sort_values(["single", "size_MB"]).iloc[0]
    return best.name, int(best.spots), float(best.size_MB), best.LibraryLayout


m = pd.read_csv("../../census/cryptococcus_sporothrix_manifest.tsv", sep="\t")
m = m[(m.species == "gattii") & (m.independent == 1) & m.runs.notna()]
m["country"] = m.country.fillna("").str.strip()
rng = np.random.default_rng(SEED)
rows = []
for country, n in QUOTA.items():
    g = m[m.country == country]
    g = g.iloc[rng.permutation(len(g))]
    for _, r in g.iterrows():
        p = pick_run(r.runs)
        if p:
            rows.append(dict(biosample=r.biosample, strain=r.strain, country=country, run=p[0], spots=p[1], size_MB=p[2], layout=p[3]))
        if sum(x["country"] == country for x in rows) == n:
            break
out = pd.DataFrame(rows)
out.to_csv("vgi_samples.tsv", sep="\t", index=False)
print(out.groupby("country").size().to_string())
print(f"{len(out)} isolates, {out.size_MB.sum() / 1000:.1f} GB to download")
