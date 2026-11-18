"""Collect published lineage labels and link them to census read sets.

Sources (labels/raw/):
  Rhodes et al. 2017 Genetics, Table S1 (188 isolates; VNI/VNII/VNB, diploids and hybrids)
  Desjardins et al. 2017 Genome Res, Supplemental Table S1 (387 isolates; VNI/VNII/VNBI/VNBII;
    saved by hand: PMC blocks scripted downloads)
  Ashton et al. 2019 Nat Commun, Supplementary Data 1 (884 isolates incl. 185 from Desjardins 2017;
    VNIa/b/c sub-clades, VNII; mash species ID)

Output: labels/literature_labels.tsv, one row per published isolate.
Usage: micromamba run -n cladeprint python labels/build_labels.py
"""
import csv
import json
import re
import pandas as pd

ACC = re.compile(r"\b(SRX\d+|ERX\d+|SRS\d+|ERS\d+|SAM[NEA]+\d+|SRR\d+|ERR\d+)\b")

# accession -> (biosample, [runs]) from all census runinfo tables
by_key, runs_of = {}, {}
for f in ["census/runinfo_5207.csv", "census/runinfo_29907.csv", "census/runinfo_1884637.csv"]:
    for r in csv.DictReader(open(f)):
        bs = r.get("BioSample") or ""
        if not bs:
            continue
        runs_of.setdefault(bs, set()).add(r["Run"])
        for k in ("Run", "Experiment", "Sample", "BioSample"):
            if r.get(k):
                by_key[r[k]] = bs
census = json.load(open("census/classified.json"))


def link(acc_text):
    accs = ACC.findall(str(acc_text))
    bss = sorted({by_key[a] for a in accs if a in by_key})
    return accs, bss


def norm(lin):
    """Map a published label to (lineage, ploidy_note): haploid lineages stay; others are flagged."""
    lin = str(lin).strip()
    m = re.match(r"(Diploid|Hybrid)\s+(\S+)", lin)
    if m:
        return m.group(2), m.group(1).lower()
    return lin, "haploid"


rows = []
r = pd.read_excel("labels/raw/rhodes2017/327TableS1.xlsx")
for _, x in r.iterrows():
    lin, note = norm(x["Lineage"])
    accs, bss = link(x["SRA accession"])
    rows.append(dict(source="Rhodes2017", isolate=x["Isolate"], label_raw=x["Lineage"], lineage=lin,
                     ploidy=note, species_check="", country=x["Country"], accessions=";".join(accs),
                     biosamples=";".join(bss)))

d = pd.read_excel("labels/raw/desjardins2017/Supplemental_Table_S1.xlsx", header=1)
for _, x in d.iterrows():
    accs, bss = link(x["SRA Accession"])
    sub = x["VNI subdivision"]
    rows.append(dict(source="Desjardins2017", isolate=x["Strain"],
                     label_raw=f'{x["Group"]} {sub}' if isinstance(sub, str) else x["Group"],
                     lineage=sub if isinstance(sub, str) and sub.startswith("VNI") else x["Group"],
                     ploidy="haploid", species_check="", country=x["Country of origin"],
                     accessions=";".join(accs), biosamples=";".join(bss)))

a = pd.read_excel("labels/raw/ashton2019/41467_2019_10092_MOESM5_ESM.xlsx", header=2)
for _, x in a.iterrows():
    sp = str(x.iloc[2])
    accs, bss = link(f'{x["ENA accession"]} {x["Sequencing ID"]}')  # Desjardins rows: SRR in Sequencing ID
    sub = str(x["Sub-clade"])
    rows.append(dict(source="Ashton2019" if x["Study"] == "Ashton" else "Desjardins2017(via Ashton)",
                     isolate=x["Sequencing ID"], label_raw=sub,
                     lineage=re.match(r"(VN[A-Za-z]+)", sub).group(1) if sub.startswith("VN") else sub,
                     ploidy="haploid", species_check="gattii" if "gattii" in sp else
                     ("AD hybrid" if "hybrid" in sp else "neoformans"),
                     country=x["Country of origin"], accessions=";".join(accs), biosamples=";".join(bss)))

lab = pd.DataFrame(rows)
lab["n_biosamples"] = lab.biosamples.str.count("SAM")
lab["runs"] = lab.biosamples.map(lambda s: ";".join(sorted(set().union(*[runs_of.get(b, set()) for b in s.split(";") if b]))) if s else "")
lab["census_species"] = lab.biosamples.map(lambda s: ";".join(sorted({census[b]["sp"] for b in s.split(";") if b in census})))
lab["census_indep"] = lab.biosamples.map(lambda s: any(census.get(b, {}).get("indep") for b in s.split(";") if b))
lab.to_csv("labels/literature_labels.tsv", sep="\t", index=False)

print(f"{len(lab)} published isolates; linked to a census BioSample: {(lab.n_biosamples > 0).sum()}; "
      f"more than one BioSample: {(lab.n_biosamples > 1).sum()}")
print("\nUnlinked by source:")
print(lab[lab.n_biosamples == 0].groupby("source").size().to_string())
print("\nHaploid, species-confirmed labels by lineage (linked / total):")
ok = lab[(lab.ploidy == "haploid") & (lab.species_check != "gattii") & (lab.species_check != "AD hybrid")]
print(ok.groupby("lineage").agg(total=("isolate", "size"), linked=("n_biosamples", lambda s: (s > 0).sum())).to_string())
print("\nCensus species of linked isolates:")
print(lab[lab.n_biosamples > 0].census_species.value_counts().to_string())
