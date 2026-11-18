"""One row per assay for Sittinan's include decision (main table / supplement / exclude).
Proposals and one-line in-silico results from the R5 audits (census unitigs + NCBI assemblies,
19 Oct). Writes include_review.tsv and .xlsx; the decision column is filled by hand, then copied
to the oligo rows of assay_list_for_review_audited.tsv (include column)."""
import csv
from collections import OrderedDict
import pandas as pd

M, S, X = "main", "supplement", "exclude"
P = {  # assay -> (proposal, in-silico result, reason)
 "Mitchell 1994 CN rDNA species PCR": (S, "complex 0.88 (=all complete assemblies; 25 scaffolds lack rDNA)", "pan-complex control"),
 "Rappelli 1998 nested CSF PCR": (S, "complex 0.88 (rDNA artifact as above)", "pan-complex control"),
 "Bialek 2002 nested 18S PCR + LightCycler FRET": (S, "complex 0.88 (rDNA artifact)", "pan-complex control"),
 "Aoki 1999 serotype A RAPD-derived PCR (CNa-70)": (M, "C. neoformans 0.99, non-target 0.10", "species assay, works"),
 "Aoki 1999 serotype A RAPD-derived PCR (CNa-29)": (M, "C. neoformans 1.00, non-target 0.10 strict / 0.45 soft", "species assay; cross-reaction worth showing"),
 "Aoki 1999 serotype B PCR (CNb-49)": (M, "C. deuterogattii 14/14, C. gattii 0/6: VGII-only, not 'serotype B'", "finding: narrower than published"),
 "Aoki 1999 preliminary serotype AD PCR (CNad-4)": (S, "C. deneoformans 1.00, AD hybrids 6/9, non-target 0", "preliminary in source"),
 "Aoki 1999 preliminary serotype D PCR (CNd-2)": (S, "no product: both primers on the same strand in all 18 D genomes", "not reproducible in silico; preliminary"),
 "Enache-Angoulvant 2007 CAP59 PCR-REA / fragment-length serotyping": (M, "sizes A 390 / D 380 / B 370 / C 360 (published 387/375/373/359)", "serotyping, reproduces"),
 "Feng 2013 STR1 singleplex species PCR": (M, "270 / 220 / 170 bp (published 274/224/170)", "species assay, reproduces"),
 "Halliday 1999 MFalpha/660 mating-type + variety coamplification": (M, "MFalpha 110 (109); 660-locus 220 only in C. gattii (216)", "reproduces"),
 "Halliday 1999 STE12alpha mating-type PCR": (S, "complex 0.88", "MAT only: no MAT labels to judge"),
 "Lengeler 2001 serotype/mating-type allele panel": (M, "serotype pairs work (CLA4/CNA1/GPA1-D 0.77-1.00); MFalpha2 D-specific (terminal mismatch); STE11alpha-D works after tail trim; JOHE2189 printed 3' end wrong; STE20alpha-D cross-reacts 0.87", "core serotype x MAT panel; several method points"),
 "Yan 2002 STE20/STE12 serotype + mating-type PCR": (M, "STE20Aa 0.084 + STE20Aalpha 0.916 = 1.000; STE20Dalpha cross-reacts 0.87", "widely used; internal check"),
 "Carvalho 2007 multiplex serotype+MAT (alphaA/aD and aA/alphaD)": (S, "MATalpha-A 0.92 (lenient); MATa-A JOHE7270 has an extra GAG vs all 19 genomes", "built from Yan/Lengeler primers"),
 "Tay 2022 Cryptococcus cyt b real-time PCR + HRM": (S, "primers 82 bp in the 120 assemblies with an mt genome; published probe fails", "mitochondrial; probe discrepancy"),
 "Gago 2014 multiplex real-time PCR (molecular beacon)": (S, "complex 0.89", "beacon scored on its loop only (approximate)"),
 "Mbangiwa 2024 QSP1 qPCR (C. neoformans)": (M, "C. neoformans 0.97, non-target 0.10", "modern species qPCR"),
 "Mbangiwa 2024 QSP1 qPCR (C. gattii)": (M, "gattii complex 0.92, non-target 0", "modern species qPCR"),
 "Mbangiwa 2024 QSP1 qPCR (C. deneoformans)": (M, "C. deneoformans 0.92, non-target 0", "modern species qPCR"),
 "Mbangiwa 2024 pan-Cryptococcus 28S qPCR": (S, "complex 0.88 (rDNA artifact)", "pan-complex control"),
 "URA5 PCR for URA5-RFLP molecular typing": (X, "", "duplicate of URA5-RFLP (Meyer et al. 2003)"),
 "Gao 2026 duplex real-time DNA RPA (C. neoformans channel)": (S, "census: C. neoformans binds well", "RPA scored as PCR (approximate)"),
 "Lucas 2010 CAP59 allele-specific LAMP": (S, "not audited", "LAMP: no simple product model"),
 "Ito-Kuwa 2007 LAC1/CAP64 multiplex serotype PCR": (S, "CAP64 works (A 421 / D 397); LAC1 D 763+1051, A bands not reproducible", "partly not reproducible; source not accessible"),
 "Cogliati 2015 C. gattii sero-mating-type multiplex PCR": (S, "band patterns per gattii species; none in C. neoformans", "MAT part cannot be judged"),
 "Leal 2008 Cryptococcus species multiplex PCR": (S, "CNa 700 in C. neoformans; CNb 450 only in C. deuterogattii", "multiplex of the Aoki primers (same result)"),
 "Veron 2009 Cryptococcus TaqMan real-time PCR": (S, "complex 0.87-0.89", "pan-complex control"),
 "Meyer 2003 URA5-RFLP original description": (X, "", "not retrieved"),
 "PCR-fingerprinting (M13)": (S, "not audited", "multi-band fingerprint (list as not auditable)"),
 "PCR-fingerprinting (microsatellite)": (S, "not audited", "multi-band fingerprint (list as not auditable)"),
 "URA5-RFLP (Meyer et al. 2003)": (M, "census: primers fail on ~5% of C. neoformans; assemblies 0.93; digest patterns (r5_rflp)", "standard typing method"),
 "PLB1-RFLP (Latouche et al. 2003)": (M, "assemblies 0.98; digest patterns (r5_rflp)", "standard typing method"),
 "QSP1 species-specific qPCR (Mbangiwa et al. 2024)": (X, "", "duplicate listing of the four Mbangiwa assays"),
 "Gago 2011 ITS1 high-resolution melting PCR": (S, "complex 0.90", "pan-complex (HRM discriminates downstream)"),
 "Lucas 2010 CAP59 allele-specific LAMP (duplicate of A103)": (X, "", "duplicate"),
}
for v in ("VGI", "VGII", "VGIIa", "VGIIb", "VGIIc", "VGIII", "VGIV"):
    P[f"MAMA SYBR real-time, {v} (Kelley et al. 2014)"] = (M, "head-to-head with our designed primers (census)", "lineage assay; comparison table")

rows = list(csv.DictReader(open("labels/assays/assay_list_for_review_audited.tsv"), delimiter="\t"))
by = OrderedDict()
for r in rows:
    by.setdefault(r["assay"], []).append(r)
out = []
for a, rs in by.items():
    prop, res, why = P[a]
    tg = "; ".join(dict.fromkeys(r["intended_target"] for r in rs))
    out.append({"assay": a, "n_oligos": len(rs), "family": rs[0]["family"], "intended_target": tg[:120],
                "in_silico_result": res, "proposed": prop, "reason": why, "your_decision": ""})
d = pd.DataFrame(out)
d.to_csv("labels/assays/include_review.tsv", sep="\t", index=False)
with pd.ExcelWriter("labels/assays/include_review.xlsx") as w:
    d.to_excel(w, index=False, sheet_name="include")
    ws = w.sheets["include"]
    for col, width in zip("ABCDEFGH", (55, 9, 22, 45, 70, 12, 45, 14)):
        ws.column_dimensions[col].width = width
    ws.freeze_panes = "B2"
print(d.proposed.value_counts().to_dict())
