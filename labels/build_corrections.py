"""Corrected-labels table (labels/corrections.tsv): published or database labels that are wrong,
with the evidence. `literature_labels.tsv` stays as published, so every correction is traceable.

Levels:
  lineage            a published lineage label contradicted by the same paper and by our data
  species            NCBI species wrong (a different species of the complex)
  species_broad_name filed under "C. gattii", the pre-2015 name of the whole species complex;
                     the run belongs to another species of the complex (not wrong, but ambiguous)
  conflict_resolved  sources disagree; resolved by a source's own finer data
verified = yes when our own data support the correction; partial when they support it only in part;
no when it rests on the source only.

Sources: manual entries below (pilot and VNII analyses); Ashton 2019's species check
(literature_labels.tsv, species_check), verified with our species screen (pilot/screen/ashton_species.tsv);
the donor search (vnii/donor/vgi_dist.tsv).
Usage: micromamba run -n cladeprint python labels/build_corrections.py
"""
import os
import pandas as pd

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
ADDED = "2026-10-15"
rows = [
    dict(isolate="102-14", run="SRR1063248", source="Rhodes2017 Table S1", published_label="Diploid VNI/VNB",
         corrected_label="Diploid VNII/VNB", level="lineage", verified="yes",
         evidence="Rhodes 2017 Table S3: heterozygous at 70,596/70,950 VNII-predictive, 12,785/12,875 VNB-predictive, "
                  "59/32,631 VNI-predictive sites; text: 'Cng9, PMHc1045, and 102-14 have similar frequencies of VNII and "
                  "VNB alleles'; cladeprint scores VNII 0.83 + VNBII 0.77; carries 84% of the VNII gattii-type block"),
    dict(isolate="ZC0022", run="ERR1520679", source="NCBI BioSample SAMEA4067825 (PRJEB14735)",
         published_label="C. neoformans", corrected_label="C. gattii", level="species", verified="yes",
         evidence="species screen (pilot/screen/species.tsv): 81.3% of reads map best to WM276 (C. gattii), 0.1% to H99"),
    dict(isolate="NRHc5024A.ENR.CLIN.ISO", run="SRR1944711", source="NCBI census (organism C. bacillisporus)",
         published_label="C. bacillisporus", corrected_label="C. neoformans VNBI", level="species", verified="yes",
         evidence="Desjardins 2017 labels it VNBI; our pilot genotyping calls VNBI (pilot/pilot_lineages.tsv)"),
    dict(isolate="MW-RSA852", run="SRR652313", source="Rhodes2017 (Hybrid VNII/VNB) vs Desjardins2017 (VNII)",
         published_label="VNII/VNB hybrid | VNII", corrected_label="Haploid VNII/VNB hybrid, 86% VNII",
         level="conflict_resolved", verified="no",
         evidence="Rhodes 2017 Table S2 ancestry VNII 0.860 / VNB 0.140, 200 SNPs from CCTP51 (VNII 0.857 / VNB 0.143), "
                  "so both labels are partly right. Our run SRR652313 is too poor to check (24% of reads map); CCTP51 scores "
                  "VNII 0.95 with no VNB signal (minor untrained ancestry)"),
]

lab = pd.read_csv("labels/literature_labels.tsv", sep="\t", dtype=str).fillna("")
ash = lab[(lab.source == "Ashton2019") & lab.species_check.isin(["gattii", "AD hybrid"]) & (lab.census_species == "neoformans")]
# our species screen of these runs (pilot/species_screen.sh, OUT=screen/ashton_species.tsv)
scr = pd.read_csv("pilot/screen/ashton_species.tsv", sep="\t").set_index("run")
GATTII = ["gattii", "deuterogattii", "bacillisporus", "tetragattii", "decagattii"]
for _, r in ash.iterrows():
    s = scr.loc[r.runs] if r.runs in scr.index else None
    if s is None:
        label, verified, ours = ("C. gattii species complex" if r.species_check == "gattii" else "AD hybrid"), "no", "not screened"
    elif r.species_check == "AD hybrid":
        ok = s.neoformans > 0.2 and s.deneoformans > 0.2
        label, verified = "AD hybrid (C. neoformans x C. deneoformans)", "yes" if ok else "no"
        ours = f"our screen: {s.neoformans:.0%} of reads C. neoformans, {s.deneoformans:.0%} C. deneoformans"
    else:
        best = max(GATTII, key=lambda g: s[g])
        ok = s[best] > 0.4 and s.neoformans < 0.05
        label, verified = f"C. {best}", "yes" if ok else "no"
        ours = f"our screen: {s[best]:.0%} of reads C. {best}, {s.neoformans:.1%} C. neoformans"
    rows.append(dict(isolate=r.isolate, run=r.runs, source="NCBI census (organism C. neoformans)",
                     published_label="C. neoformans", corrected_label=label, level="species", verified=verified,
                     evidence=f"Ashton 2019 species check: {r.species_check} ({r.country}); {ours} (pilot/screen/ashton_species.tsv)"))

vgi = pd.read_csv("vnii/donor/vgi_dist.tsv", sep="\t", dtype={"strain": str}).fillna("")
for _, r in vgi[vgi.species != "gattii"].iterrows():
    far = r.species_dist > 10  # same-species isolates sit at ~1-8/kb; far from every panel genome
    rows.append(dict(isolate=r.strain if r.strain not in ("", "missing") else r.run, run=r.run,
                     source="NCBI census (organism C. gattii)", published_label="C. gattii",
                     corrected_label=f"C. {r.species}" + (" (tentative)" if far else ""), level="species_broad_name",
                     verified="partial" if far else "yes",
                     evidence=f"VNII-block-region consensus is {r.species_dist}/kb from the C. {r.species} panel genome, "
                              f"next species {r.next_species_dist}/kb ({r.country}; vnii/donor/vgi_dist.tsv)"
                              + ("; far from every panel genome: a divergent lineage or a hybrid, check genome-wide" if far else "")))

out = pd.DataFrame(rows).assign(added=ADDED)
out = out[["isolate", "run", "source", "published_label", "corrected_label", "level", "verified", "evidence", "added"]]
out.to_csv("labels/corrections.tsv", sep="\t", index=False)
print(out.groupby(["level", "verified"]).size().to_string())
print(f"{len(out)} corrections")
