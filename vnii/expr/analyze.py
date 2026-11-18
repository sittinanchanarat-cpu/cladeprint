"""Expression of the C. gattii-type block genes in HC1 (VNII), Chen et al. 2014 RNA-seq (GSE51573).

Counts: featureCounts on CDS exons of miniprot models (genes.saf), MAPQ >= 10, unstranded.
- YPD at full depth (align.sh -> ypd_full_s{0,1,2}.tsv; s1/s2 = sense/antisense, strandedness check)
- YPD, ex vivo CSF and in vivo human CSF as 10M-spot subsamples (align_fast.sh -> sub/<run>.counts.tsv)
TPM per gene; each block gene's TPM and percentile among all 6,756 genes with a model; CSF/YPD
log2 fold change against the genome-wide spread (one library per condition: no replicates, no test).
Writes expression.tsv. Usage: micromamba run -n cladeprint python analyze.py
"""
import os
import numpy as np
import pandas as pd

os.chdir(os.path.dirname(os.path.abspath(__file__)))
NAMES = {"XP_003194833.1": "F4600C flank (CNAG_05835)", "XP_003194826.1": "F4610W HOC1 a-1,6-mannosyltransferase",
         "XP_003194832.1": "F4620C uncharacterized", "XP_003194831.1": "F4630C hypothetical 807 aa",
         "XP_003194830.1": "F4640C amino-acid transporter", "XP_003194827.1": "F4650W uncharacterized",
         "XP_003194829.1": "F4660C dihydropyrimidinase", "XP_003194828.1": "F4670W 2OG-oxygenase-like",
         "XP_003195013.1": "F4680W uncharacterized", "XP_003195014.1": "F4690W RGD1 Rho-GAP",
         "XP_003195015.1": "F4700C flank (partial model)"}


def counts(path):
    d = pd.read_csv(path, sep="\t", comment="#").set_index("Geneid")
    return d.iloc[:, -1], d.Length


ypd_full, length = counts("ypd_full_s0.tsv")
sense, anti = counts("ypd_full_s1.tsv")[0].sum(), counts("ypd_full_s2.tsv")[0].sum()
print(f"YPD strandedness: {sense:,} reads on the gene strand vs {anti:,} opposite -> unstranded library")
c = pd.DataFrame({"YPD_full": ypd_full, "YPD": counts("sub/SRR1015479.counts.tsv")[0],
                  "exvivo_CSF": counts("sub/SRR1015480.counts.tsv")[0],
                  "invivo_CSF": counts("sub/SRR1015481.counts.tsv")[0]})
rpk = c.div(length / 1000, axis=0)
tpm = rpk / rpk.sum() * 1e6
pct = tpm.rank(pct=True) * 100
print("assigned reads/pairs:", c.sum().astype(int).to_dict(), "| median TPM:", tpm.median().round(1).to_dict())

r = np.corrcoef(np.log2(tpm.YPD_full + 1), np.log2(tpm.YPD + 1))[0, 1]
print(f"depth check, YPD full vs 10M subsample: Pearson r (log2 TPM+1) = {r:.4f}")

lfc = {k: np.log2((tpm[k] + 1) / (tpm.YPD + 1)) for k in ("exvivo_CSF", "invivo_CSF")}
for k, v in lfc.items():
    print(f"genome-wide log2FC {k}/YPD: 5th-95th percentile {np.percentile(v, 5):.2f} to {np.percentile(v, 95):.2f}")

block = [g for g in tpm.index if g.startswith("BLOCK:")]
rows = []
for g in block:
    row = {"gene": NAMES[g.split(":", 1)[1]], "cds_bp": int(length[g])}
    for k in c.columns:
        row[f"{k}_reads"] = int(c.at[g, k])
        row[f"{k}_TPM"] = round(tpm.at[g, k], 1)
        row[f"{k}_pct"] = int(round(pct.at[g, k]))
    for k, v in lfc.items():
        row[f"{k}_log2FC"] = round(v[g], 2)
        row[f"{k}_log2FC_pct"] = int(round((v <= v[g]).mean() * 100))
    rows.append(row)
out = pd.DataFrame(rows)
out.to_csv("expression.tsv", sep="\t", index=False)
show = ["gene"] + [f"{k}_{x}" for k in ("YPD_full", "exvivo_CSF", "invivo_CSF") for x in ("TPM", "pct")] + \
       [f"{k}_log2FC" for k in lfc]
pd.set_option("display.width", 250)
print(out[show].to_string(index=False))
