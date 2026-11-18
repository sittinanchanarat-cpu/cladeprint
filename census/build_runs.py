"""Census run list for the AWS run: one sequencing run per independent isolate.

From the census manifest: biosamples marked independent, of natural or serial origin (serial =
the first isolate of a patient series), of a species of the Cryptococcus neoformans/gattii complex
(set cryptococcus) or of Sporothrix (set sporothrix, for its own bundle later). The NCBI species is
only recorded: species are assigned from the data (species screen in fetch_count.sh), since NCBI
labels are often wrong (labels/corrections.tsv). Excluded: class "unverified" (260 C. neoformans,
see census/classify.py), lab and hybrid-study material, and long-read-only isolates.

One run per biosample: short reads (Illumina, DNBSEQ, BGISEQ), paired preferred, >= 1.2M spots
(enough for the 1M-spot subsample), the smallest download among those; else the run with most spots.
Writes census/census_runs.tsv. Usage: micromamba run -n cladeprint python census/build_runs.py
"""
import csv, os
import pandas as pd

os.chdir(os.path.dirname(os.path.abspath(__file__)))
CRYPTO = ["neoformans", "deneoformans", "gattii", "deuterogattii", "bacillisporus", "tetragattii", "decagattii", "hyracis"]
SPORO = ["brasiliensis", "schenckii", "globosa", "other_Sporothrix"]
SHORT = {"ILLUMINA", "DNBSEQ", "BGISEQ"}
MIN_SPOTS = 1_200_000


def read_runinfo(path):
    """The runinfo files are concatenated pages, each with its own header and column count."""
    recs, head = [], None
    with open(path, newline="") as fh:
        for row in csv.reader(fh):
            if row and row[0] == "Run":
                head = row
            elif row and head and len(row) == len(head):
                recs.append(dict(zip(head, row)))
    d = pd.DataFrame(recs)[["Run", "spots", "size_MB", "LibraryLayout", "Platform", "avgLength"]]
    for c in ("spots", "size_MB", "avgLength"):
        d[c] = pd.to_numeric(d[c], errors="coerce").fillna(0)
    return d


ri = pd.concat([read_runinfo(f) for f in os.listdir(".") if f.startswith("runinfo_")]).drop_duplicates("Run").set_index("Run")
m = pd.read_csv("cryptococcus_sporothrix_manifest.tsv", sep="\t", dtype=str).fillna("")
m = m[(m.independent == "1") & m["class"].isin(["natural", "serial"]) & m.species.isin(CRYPTO + SPORO) & (m.runs != "")]

rows, skipped = [], {}
for _, b in m.iterrows():
    r = ri.reindex([x for x in b.runs.split(";") if x in ri.index])
    r = r[r.Platform.isin(SHORT) & (r.spots > 0)]
    if r.empty:
        skipped["no short-read run"] = skipped.get("no short-read run", 0) + 1
        continue
    ok = r[r.spots >= MIN_SPOTS].assign(single=lambda x: x.LibraryLayout != "PAIRED")
    best = ok.sort_values(["single", "size_MB"]).iloc[0] if len(ok) else r.sort_values("spots").iloc[-1]
    rows.append(dict(
        biosample=b.biosample, run=best.name, set="cryptococcus" if b.species in CRYPTO else "sporothrix",
        ncbi_species=b.species, origin=b["class"], strain=b.strain, country=b.country, bioproject=b.bioproject,
        platform=best.Platform, layout=best.LibraryLayout, spots=int(best.spots), size_MB=float(best.size_MB),
        read_length=int(best.avgLength), enough_spots=bool(best.spots >= MIN_SPOTS)))
out = pd.DataFrame(rows).drop_duplicates("run")
out.to_csv("census_runs.tsv", sep="\t", index=False)
print(out.groupby(["set", "ncbi_species"]).size().to_string())
print(f"{len(out)} runs, {out.size_MB.sum() / 1e6:.2f} TB to download; fewer than {MIN_SPOTS:,} spots: "
      f"{(~out.enough_spots).sum()}; paired {(out.layout == 'PAIRED').mean():.0%}; skipped: {skipped}")
