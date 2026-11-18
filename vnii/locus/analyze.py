"""Per-gene depth and read identity at the HOC1-RGD1 locus, from the competitive H99 + WM276 remaps.

For each run and gene: primary reads overlapping the gene (MAPQ >= 20), depth relative to the run's
genome-wide median (H99 chr1), and read identity = 1 - sum(NM) / sum(aligned bases). A VNII allele of
C. gattii type should recruit its reads to the WM276 gene at ~98% and leave the H99 gene empty.
Usage: python analyze.py RUN...   (table to locus_genes.tsv, or OUT=file)
"""
import os, subprocess, sys
import numpy as np, pandas as pd

H = "neoformans|NC_026751.1"
G = "gattii|NC_014943.1"
GENES = [  # name, contig, start, end (H99 coordinates from NCBI Gene; WM276 from wm276.gff)
    ("HOC1 (H99)", H, 1004386, 1006751),
    ("CNAG_05837 (H99)", H, 1006601, 1007693),
    ("RGD1 5' (H99, not in block)", H, 1007283, 1008830),
    ("RGD1 3' (H99)", H, 1008831, 1012078),
    ("H99 left flank 990-1000k", H, 990000, 1000000),
    ("H99 right flank 1015-1025k", H, 1015000, 1025000),
    ("H99 baseline chr1 0.2-1.2M", "neoformans|NC_026745.1", 200000, 1200000),
]
gff = pd.read_csv("../wm276.gff", sep="\t", comment="#", header=None)
gff = gff[(gff[0] == G.split("|")[1]) & (gff[2] == "gene") & (gff[3] > 915000) & (gff[4] < 955000)]
for _, r in gff.iterrows():
    loc = r[8].split("locus_tag=")[1].split(";")[0]
    GENES.append((f"{loc} (WM276)", G, int(r[3]), int(r[4])))


def view(bam, region, q=20):
    out = subprocess.run(["samtools", "view", "-F", "0x904", "-q", str(q), bam, region],
                         capture_output=True, text=True, check=True).stdout
    nm = al = n = 0
    for line in out.splitlines():
        f = line.split("\t")
        ops = __import__("re").findall(r"(\d+)([MIDNSHP=X])", f[5])
        al += sum(int(k) for k, o in ops if o in "M=XI")
        nm += next(int(t[5:]) for t in f[11:] if t.startswith("NM:i:"))
        n += 1
    return n, al, nm


def variants(bam, ctg, s, e):
    """Haploid calls against the reference the reads chose: SNPs + indel bases per kb of covered gene."""
    call = subprocess.run(f"bcftools mpileup -q 20 -Q 20 -a AD -f ref2.fna -r '{ctg}:{s}-{e}' {bam} 2>/dev/null"
                          f" | bcftools call -m --ploidy 1 -V indels 2>/dev/null | bcftools view -H -i 'QUAL>=30 && DP>=3' -", shell=True,
                          capture_output=True, text=True).stdout
    snps = len([l for l in call.splitlines() if l.split("\t")[4] != "."])
    cov = subprocess.run(["samtools", "depth", "-Q", "20", "-r", f"{ctg}:{s}-{e}", bam], capture_output=True, text=True).stdout
    covered = sum(1 for l in cov.splitlines() if int(l.split("\t")[2]) >= 3)
    return round(1000 * snps / covered, 1) if covered > 200 else np.nan


def median_depth(bam):
    out = subprocess.run(["samtools", "depth", "-a", "-Q", "20", "-r", "neoformans|NC_026745.1:200000-1200000", bam],
                         capture_output=True, text=True, check=True).stdout
    return float(np.median([int(l.split("\t")[2]) for l in out.splitlines()]))


rows = []
for run in sys.argv[1:]:
    bam = f"{run}.bam"
    med = median_depth(bam)
    for name, ctg, s, e in GENES:
        n, al, nm = view(bam, f"{ctg}:{s}-{e}")
        rows.append(dict(run=run, gene=name, reads=n, rel_depth=round(al / (e - s + 1) / med, 2) if med else np.nan,
                         identity=round(1 - nm / al, 4) if al else np.nan,
                         snps_per_kb=variants(bam, ctg, s, e)))
df = pd.DataFrame(rows)
order = [g[0] for g in GENES]
for col in ("reads", "rel_depth", "identity", "snps_per_kb"):
    t = df.pivot(index="gene", columns="run", values=col).loc[order, sys.argv[1:]]
    print(f"\n{col}\n{t.to_string()}")
df.to_csv(os.environ.get("OUT", "locus_genes.tsv"), sep="\t", index=False)
