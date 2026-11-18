"""Sample sheet and config for `cladeprint run` on the census (after fetch_count.sh).

Set cryptococcus (default), levels: species -> lineage (VNI, VNII, VNB) -> sublineage (VNIa, VNIb,
VNIc, VNBI, VNBII). Set sporothrix: species only, from the Sporothrix panel (no hybrid rule).
- species from the data (OUT/species/RUN.tsv), as shares of the reads assigned to a species: the top
  species >= 80% with the next below 10%; C. neoformans and C. deneoformans both >= 20% is an AD
  hybrid (species "AD-hybrid", does not train); fewer than 5% of reads assigned, or anything else,
  is left blank and does not train. NCBI labels are not
  used (labels/corrections.tsv).
- lineage and sublineage from the literature labels (labels/literature_labels.tsv, by biosample),
  only for C. neoformans by the screen; corrections (labels/corrections.tsv, level lineage) first,
  then Desjardins 2017, Ashton 2019, Rhodes 2017. Hybrid or diploid labels ("/") leave the lineage
  blank and the sample does not train; disagreeing sources at the lineage level leave it blank.
- train: auto (cladeprint's QC: coverage, k-mer excess) unless set to no above.
- unmapped: share of the screened reads that map to no panel species. Over 50% is flagged in the note
  (old error-rich runs, contamination or a species outside the panel) for review; it does not by
  itself stop training.
Only runs with a finished database are included.
Config: markers.max_out = 1 (a marker may be in one out-group run: survives a single mislabel or
contaminated run; README.md). Writes OUT/samples.tsv and OUT/config.toml; for sporothrix OUT/sporothrix/ with a db link to the shared
OUT/db.
Clade weighting: with a third argument, the `cluster` column comes from that column of
OUT/clusters/clusters.tsv (`cladeprint clusters` over all runs, e.g. cluster_10); its mixed_share is
copied too, for review. Without it, no cluster column (unweighted).
Usage: micromamba run -n cladeprint python census/make_sheet.py [OUT] [cryptococcus|sporothrix] [CLUSTER_COLUMN]
"""
import os, sys
import pandas as pd

os.chdir(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1] if len(sys.argv) > 1 else "out"
SET = sys.argv[2] if len(sys.argv) > 2 else "cryptococcus"
assert SET in ("cryptococcus", "sporothrix"), SET
CLUSTER_COL = sys.argv[3] if len(sys.argv) > 3 else None
PANEL, GENUS = ("complex", "C.") if SET == "cryptococcus" else ("sporothrix", "S.")
SPECIES = open(f"../pilot/ref/panel/{PANEL}.species.txt").read().split()  # column order of species/RUN.tsv
HIGH_UNMAPPED = 0.5
SUB = {"VNIa": "VNI", "VNIb": "VNI", "VNIc": "VNI", "VNBI": "VNB", "VNBII": "VNB"}
PRECEDENCE = ["Desjardins2017", "Desjardins2017(via Ashton)", "Ashton2019", "Rhodes2017"]

runs = pd.read_csv("census_runs.tsv", sep="\t", dtype={"strain": str}).fillna("")
runs = runs[runs.set == SET]
done = {f[: -len(".stats.toml")] for f in os.listdir(f"{OUT}/db") if f.endswith(".stats.toml")}
runs = runs[runs.run.isin(done)]


def species_call(run):
    """Species from the screen. Shares are taken among the reads assigned to a species (MAPQ >= 20,
    nuclear contigs): old error-rich runs leave most reads unmapped, but their assigned reads are
    clear-cut. At least 5% of the screened reads must be assigned."""
    p = f"{OUT}/species/{run}.tsv"
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return "", None, "no species screen"
    f = open(p).read().split("\t")
    unmapped = float(f[2])
    frac = dict(zip(SPECIES, map(float, f[4:4 + len(SPECIES)])))
    assigned = sum(frac.values())
    if assigned < 0.05:
        return "", unmapped, f"species screen: only {assigned:.1%} of reads assigned"
    share = {k: v / assigned for k, v in frac.items()}
    if SET == "cryptococcus" and share["neoformans"] >= 0.2 and share["deneoformans"] >= 0.2:
        return "AD-hybrid", unmapped, "AD hybrid by species screen"
    top = sorted(share.items(), key=lambda kv: -kv[1])
    if top[0][1] >= 0.8 and top[1][1] < 0.1:
        return GENUS + top[0][0], unmapped, ""
    return "", unmapped, f"species unclear ({top[0][0]} {top[0][1]:.2f}, {top[1][0]} {top[1][1]:.2f} of assigned reads)"


lab = pd.read_csv("../labels/literature_labels.tsv", sep="\t", dtype=str).fillna("")
cor = pd.read_csv("../labels/corrections.tsv", sep="\t", dtype=str).fillna("")
cor = cor[cor.level == "lineage"].set_index("run").corrected_label.to_dict()
by_bs = {}
for _, r in lab.iterrows():
    for bs in r.biosamples.split(";"):
        if bs:
            by_bs.setdefault(bs, []).append((PRECEDENCE.index(r.source) if r.source in PRECEDENCE else 9, r.lineage, r.source))


def lineage_call(biosample, run):
    if run in cor:
        raw = cor[run].replace("Diploid ", "").replace("Haploid ", "")
        return ("", "", f"hybrid/diploid label {raw} (corrected)") if "/" in raw else (SUB.get(raw, raw), raw if raw in SUB else "", "")
    labs = sorted(by_bs.get(biosample, []))
    if not labs:
        return "", "", ""
    if any("/" in lin for _, lin, _ in labs):
        hyb = next(lin for _, lin, _ in labs if "/" in lin)
        return "", "", f"hybrid/diploid label {hyb}"
    majors = {SUB.get(lin, lin) for _, lin, _ in labs}
    if len(majors) > 1:
        return "", "", f"conflicting lineage labels {sorted(majors)}"
    best = labs[0][1]
    return SUB.get(best, best), best if best in SUB else "", ""


rows = []
for _, r in runs.iterrows():
    sp, unmapped, why_sp = species_call(r.run)
    why_un = f"high unmapped ({unmapped:.0%} of screened reads)" if unmapped is not None and unmapped > HIGH_UNMAPPED else ""
    lin, sub, why_lin = lineage_call(r.biosample, r.run) if sp == "C.neoformans" else ("", "", "")
    train = "no" if (sp in ("", "AD-hybrid") or why_lin.startswith("hybrid")) else "auto"
    rows.append(dict(sample=r.run, reads="fetched", species=sp, lineage=lin, sublineage=sub, train=train,
                     biosample=r.biosample, ncbi_species=r.ncbi_species, country=r.country,
                     unmapped="" if unmapped is None else f"{unmapped:.4f}",
                     note="; ".join(x for x in (why_sp, why_lin, why_un) if x)))
sheet = pd.DataFrame(rows)
if CLUSTER_COL:
    clu = pd.read_csv(f"{OUT}/clusters/clusters.tsv", sep="\t", dtype=str).set_index("sample")
    assert CLUSTER_COL in clu.columns, f"{CLUSTER_COL} not in clusters.tsv: {list(clu.columns)}"
    missing = sorted(set(sheet["sample"]) - set(clu.index))
    assert not missing, f"{len(missing)} runs missing from clusters.tsv, e.g. {missing[:3]}"
    sheet["cluster"] = sheet["sample"].map(clu[CLUSTER_COL])
    sheet["mixed_share"] = sheet["sample"].map(clu["mixed_share"])
levels = ["species", "lineage", "sublineage"] if SET == "cryptococcus" else ["species"]
if SET == "sporothrix":
    sheet = sheet.drop(columns=["lineage", "sublineage"])
dest = OUT if SET == "cryptococcus" else f"{OUT}/sporothrix"
os.makedirs(dest, exist_ok=True)
if SET == "sporothrix" and not os.path.lexists(f"{dest}/db"):
    os.symlink("../db", f"{dest}/db")  # databases are shared with the cryptococcus set
sheet.to_csv(f"{dest}/samples.tsv", sep="\t", index=False)
open(f"{dest}/config.toml", "w").write(
    f'samples = "samples.tsv"\noutput = "."\nlevels = {levels}\n\n'.replace("'", '"')
    + "[markers]\nmax_out = 1\n\n"  # tolerate one mislabelled or contaminated out-group run
    + "[resources]\nthreads = 0\nmemory_gb = 8\n")
print(sheet.groupby(levels, dropna=False).size().to_string())
print(f"{len(sheet)} samples, train=auto {int((sheet.train == 'auto').sum())}; species vs NCBI disagree: "
      f"{int((sheet.species.str.replace(GENUS, '', regex=False) != sheet.ncbi_species).sum())}; "
      f"high unmapped: {int(sheet.note.str.contains('high unmapped').sum())}")
