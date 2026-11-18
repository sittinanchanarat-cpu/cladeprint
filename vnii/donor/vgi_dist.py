"""Donor search: how close is each VGI isolate's block to the VNII block?

Per isolate (bam/RUN.bam from fetch_vgi.sh, mapped to H99 + WM276):
  - species check: share of baseline reads on WM276 chr A vs H99 chr 1, and SNPs/kb against WM276
    over WM276 chr A 0.2-1.2 Mb (VGI ~ a few/kb; other gattii-complex species ~50+)
  - block: haploid consensus over WM276 chrF 923,314-947,300 (positions with depth < 3 masked),
    then differences per kb against WM276 and against the VNII consensus (vnii_SRR18672038.fa;
    the 5 VNII consensuses agree to within 1/kb), counted only where both sequences are called.
  - species: the block consensus is compared with the 5 gattii-complex panel genomes; closest wins
    (NCBI labels are unreliable: Malawi and India runs filed as C. gattii include C. tetragattii)
A donor-like isolate has a block much closer to VNII than the ~15/kb of the public assemblies.
Writes vgi_dist.tsv. Usage: micromamba run -n cladeprint python vgi_dist.py
"""
import os, re, subprocess
from pathlib import Path
import pandas as pd

os.chdir(os.path.dirname(os.path.abspath(__file__)))
REF = "../locus/ref2.fna"
BLOCK = "gattii|NC_014943.1:923314-947300"
BASE_G = "gattii|NC_014938.1:200000-1200000"
BASE_N = "neoformans|NC_026745.1:200000-1200000"


def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout


def consensus(bam, region):
    """Haploid consensus (str) over region, 'N' where depth < 3, plus SNPs/kb over called positions."""
    vcf = f"tmp/{Path(bam).stem}.vcf.gz"
    sh(f"bcftools mpileup -q 20 -Q 20 -f {REF} -r '{region}' {bam} 2>/dev/null | bcftools call -m --ploidy 1 2>/dev/null"
       f" | bcftools view -i 'ALT!=\".\" && QUAL>=30 && DP>=3' -Oz -o {vcf} 2>/dev/null; bcftools index -f {vcf}")
    low = f"tmp/{Path(bam).stem}.low.bed"
    sh(f"samtools depth -a -Q 20 -r '{region}' {bam} | awk '$3<3{{print $1\"\\t\"$2-1\"\\t\"$2}}' > {low}")
    seq = "".join(sh(f"samtools faidx {REF} '{region}' | bcftools consensus -m {low} {vcf} 2>/dev/null").split("\n")[1:])
    nvar = int(sh(f"bcftools view -H {vcf} | wc -l"))
    called = sum(c != "N" for c in seq)
    return seq, (1000 * nvar / called if called else float("nan")), called


SPECIES = ["gattii", "deuterogattii", "bacillisporus", "tetragattii", "decagattii"]


def closest_species(seq, run):
    """Closest gattii-complex panel genome to a block consensus: (species, diffs/kb, runner-up diffs/kb)."""
    fa = f"tmp/blk_{run}.fa"
    with open(fa, "w") as f:
        f.write(f">{run}\n{seq}\n")
    out = sh("python dist.py " + fa + " " + " ".join(f"../../pilot/ref/panel/{s}.fna" for s in SPECIES))
    d = sorted((float(l.split("\t")[3]), l.split("\t")[0]) for l in out.splitlines()[1:])
    return d[0][1], d[0][0], d[1][0]


def diffs_per_kb(a, b):
    """Differences per kb between two consensus strings: minimap2 alignment, substitutions and
    indel events from the cs tag, skipping positions where either side is N."""
    with open("tmp/a.fa", "w") as f:
        f.write(f">a\n{a}\n")
    with open("tmp/b.fa", "w") as f:
        f.write(f">b\n{b}\n")
    out = sh("minimap2 -c --cs -x asm20 tmp/b.fa tmp/a.fa 2>/dev/null")
    al = df = 0
    for line in out.splitlines():
        f = line.split("\t")
        if "tp:A:P" not in f:
            continue
        cs = next(t for t in f if t.startswith("cs:Z:"))[5:]
        for op, val in re.findall(r"([:*+\-])([0-9a-z]+)", cs):
            if op == ":":
                al += int(val)
            elif op == "*":
                if "n" in val:
                    continue
                al += 1
                df += 1
            elif "n" not in val:
                df += 1
    return (1000 * df / al if al else float("nan")), al


os.makedirs("tmp", exist_ok=True)
vnii = "".join(open("vnii_SRR18672038.fa").read().split("\n")[1:])
meta = pd.read_csv("vgi_samples.tsv", sep="\t").set_index("run")
rows = []
for bam in sorted(Path("bam").glob("*.bam")):
    run = bam.stem
    nG = int(sh(f"samtools view -c -F 0x904 -q 20 {bam} '{BASE_G}'"))
    nN = int(sh(f"samtools view -c -F 0x904 -q 20 {bam} '{BASE_N}'"))
    _, base_snps, base_called = consensus(str(bam), BASE_G)
    blk, blk_snps, blk_called = consensus(str(bam), BLOCK)
    d_vnii, al = diffs_per_kb(blk, vnii)
    sp, sp_d, sp_next = closest_species(blk, run)
    m = meta.loc[run] if run in meta.index else {}
    rows.append(dict(run=run, strain=m.get("strain", ""), country=m.get("country", ""),
                     species=sp, species_dist=sp_d, next_species_dist=sp_next,
                     gattii_read_share=round(nG / max(1, nG + nN), 3), base_snps_kb_vs_WM276=round(base_snps, 1),
                     block_called_bp=blk_called, block_vs_WM276=round(blk_snps, 1),
                     block_vs_VNII=round(d_vnii, 1), compared_bp=al))
    print(rows[-1], flush=True)
res = pd.DataFrame(rows).sort_values("block_vs_VNII")
res.to_csv("vgi_dist.tsv", sep="\t", index=False)
print(res.to_string(index=False))
