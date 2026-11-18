"""Pick literature-labelled isolates to genotype: every VNII, 8 anchors each for VNBI/VNBII/VNIb/VNIc,
and the isolates whose labels conflict (between papers, or with the census species).

Per isolate, use the smallest paired Illumina run with >= 1.1M spots (less to download).
Writes labels/validation_samples.tsv in the pilot_samples.tsv format (role = "labelled").
Usage: micromamba run -n cladeprint python labels/select_validation.py
"""
import csv
import random
import pandas as pd

random.seed(20261012)
N_ANCHORS = 8
pilot = set(pd.read_csv("pilot/pilot_samples.tsv", sep="\t").biosample)

runs = {}
for f in ["census/runinfo_5207.csv", "census/runinfo_29907.csv", "census/runinfo_1884637.csv"]:
    for r in csv.DictReader(open(f)):
        if r["LibraryLayout"] == "PAIRED" and r["Platform"] in ("ILLUMINA", "BGISEQ", "DNBSEQ") \
                and int(r["spots"] or 0) >= 1_100_000:
            runs.setdefault(r["BioSample"], []).append(r)

lab = pd.read_csv("labels/literature_labels.tsv", sep="\t").fillna("")
lab = lab[lab.biosamples != ""]
lab["major"] = lab.lineage.map(lambda s: "VNI" if s.startswith("VNI") and not s.startswith("VNII") else s)

# isolates whose labels disagree: between sources, or label vs census species
conflict = set()
for bs, d in lab.groupby("biosamples"):
    if d.ploidy.ne("haploid").any() and d.ploidy.eq("haploid").any():
        conflict.add(bs)
conflict |= set(lab[(lab.census_species != "neoformans")].biosamples)

clean = lab[(lab.ploidy == "haploid") & ~lab.species_check.isin(["gattii", "AD hybrid"])
            & (lab.census_species == "neoformans") & ~lab.biosamples.isin(conflict)]
clean = clean.sort_values("source").drop_duplicates("biosamples")

pick = [(bs, "VNII") for bs in clean[clean.major == "VNII"].biosamples]
for lin in ["VNBI", "VNBII", "VNIb", "VNIc"]:
    pool = sorted(set(clean[clean.lineage == lin].biosamples) & set(runs) - pilot)
    pick += [(bs, lin) for bs in random.sample(pool, min(N_ANCHORS, len(pool)))]
pick += [(bs, "conflict") for bs in sorted(conflict)]

info = lab.drop_duplicates("biosamples").set_index("biosamples")
out, missing = [], []
for bs, lin in pick:
    if bs in pilot:
        continue
    if bs not in runs:
        missing.append((bs, lin))
        continue
    r = min(runs[bs], key=lambda r: float(r.get("size_MB") or 0) or int(r["spots"]))
    x = info.loc[bs]
    label = x.label_raw if lin == "conflict" else lin
    out.append(dict(biosample=bs, run=r["Run"], role="conflict" if lin == "conflict" else "labelled", label=label, bioproject=r["BioProject"],
                    strain=x.isolate, country=x.country, spots=r["spots"], avgLength=r["avgLength"]))
o = pd.DataFrame(out)
o.to_csv("labels/validation_samples.tsv", sep="\t", index=False)
print(f"{len(o)} runs to genotype; skipped {len(missing)} with no paired run >= 1.1M spots")
print(o.groupby(["role", "label"]).size().to_string())
print(f"total download ~{sum(float(r.get('size_MB') or 0) for bs in o.biosample for r in runs[bs] if r['Run'] in set(o.run)) / 1e3:.1f} GB")
