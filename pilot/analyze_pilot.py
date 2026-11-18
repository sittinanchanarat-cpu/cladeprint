"""Pilot analysis: SNP matrix -> pairwise distances -> anchor-based lineage calls,
clonal pairs, and QC (mapping rate, callable fraction, mixed alleles for hybrids/diploids).

Sparse samples (low mapping, or position-sorted uploads where the first 1M pairs cover
only chr1) must not veto sites: sites are chosen on well-covered samples only, and every
distance is a proportion over the sites both samples have called.

Usage: micromamba run -n cladeprint python analyze_pilot.py
"""
import os
import subprocess
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import average, fcluster
from scipy.spatial.distance import squareform

VCF, CACHE = "calls.raw.vcf.gz", "calls.snps.npz"
MITO = "NC_018792.1"
MIN_QUAL, MIN_DP, MAX_MISS = 30, 5, 0.10
GOOD_MAX_MISS = 0.50  # samples missing more than this are reported but don't choose sites
MIN_SHARED = 2000     # fewer shared sites than this -> distance flagged as unreliable
CLONAL_SNPS = 200     # genome-wide SNP distance below which two isolates count as one
# species screen (screen/species.tsv): not C. neoformans -> dropped entirely
EXCLUDE = {"ERR1520679"}  # ZC0022, C. gattii VGI
# kept as samples but not used as lineage anchors: mixed alleles suggest mixed infection/diploid
NOT_ANCHOR = {"SRR34291933", "SRR34291931"}  # MCM450-N, MCM531-C (Duke "VNII")


def load_matrix():
    if os.path.exists(CACHE) and os.path.getmtime(CACHE) > os.path.getmtime(VCF):
        z = np.load(CACHE, allow_pickle=True)
        return list(z["samples"]), z["G"], z["H"]
    samples = subprocess.run(["bcftools", "query", "-l", VCF], capture_output=True, text=True,
                             check=True).stdout.split()
    view = subprocess.Popen(["bcftools", "view", "-v", "snps", "-m2", "-M2", "-t", f"^{MITO}",
                             "-i", f"QUAL>={MIN_QUAL}", "-Ou", VCF], stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL)
    query = subprocess.Popen(["bcftools", "query", "-f", r"[%GT\t%DP\t%AD{0}\t%AD{1}\t]\n"],
                             stdin=view.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    t = pd.read_csv(query.stdout, sep="\t", header=None, na_values=".", dtype="float32",
                    usecols=range(4 * len(samples))).to_numpy()
    gt, dp, ref, alt = (t[:, k::4] for k in range(4))
    G = np.where(np.nan_to_num(dp) >= MIN_DP, gt, np.nan).astype("float32")
    ref, alt = np.nan_to_num(ref), np.nan_to_num(alt)
    tot = ref + alt
    H = (tot >= 10) & (np.minimum(ref, alt) >= 0.2 * np.maximum(tot, 1))
    np.savez(CACHE, samples=np.array(samples), G=G, H=H)
    return samples, G, H


SHEETS = ["pilot_samples.tsv", "../labels/validation_samples.tsv"]  # role: anchor / pilot / labelled / conflict
meta = pd.concat([pd.read_csv(f, sep="\t") for f in SHEETS if os.path.exists(f)]).drop_duplicates("run").set_index("run")
samples, G, H = load_matrix()
keep_s = [i for i, s in enumerate(samples) if s not in EXCLUDE]
samples, G, H = [samples[i] for i in keep_s], G[:, keep_s], H[:, keep_s]
meta.loc[meta.index.isin(NOT_ANCHOR), "role"] = "demoted"
print(f"{G.shape[0]} biallelic SNP sites (QUAL>={MIN_QUAL}, nuclear), {len(samples)} samples")

# per-sample QC
def map_rate_of(run):
    """Fraction of reads in the BAM that map to H99 (from the index, so it's quick)."""
    out = subprocess.run(["samtools", "idxstats", f"bam/{run}.bam"], capture_output=True, text=True).stdout
    rows = [l.split("\t") for l in out.splitlines()]
    mapped, unmapped = sum(int(r[2]) for r in rows), sum(int(r[3]) for r in rows)
    return mapped / max(1, mapped + unmapped)


map_rate = np.array([map_rate_of(s) for s in samples])
called = ~np.isnan(G)
miss_s = 1 - called.mean(0)
het_frac = H.sum(0) / np.maximum(1, called.sum(0))
good = miss_s <= GOOD_MAX_MISS
print(f"{good.sum()} samples with <= {GOOD_MAX_MISS:.0%} missing choose the sites; "
      f"{(~good).sum()} sparse samples are placed on whatever sites they have")

Gg = G[:, good]
keep = (np.isnan(Gg).mean(1) <= MAX_MISS) & (np.nanmax(Gg, 1) != np.nanmin(Gg, 1))
X = G[keep]
print(f"{X.shape[0]} polymorphic sites with <= {MAX_MISS:.0%} missing among good samples")

# pairwise distance = proportion of shared called sites that differ, scaled to all kept sites
C = (~np.isnan(X)).astype("float32")
A = np.nan_to_num(X).astype("float32")
shared = C.T @ C
diff = (A * C).T @ ((1 - A) * C) + ((1 - A) * C).T @ (A * C)
D = diff / np.maximum(1, shared) * X.shape[0]
np.fill_diagonal(D, 0)

# lineage = label of nearest good anchor; report the distance so weak calls are visible
# anchors with many mixed-allele sites (mixed infection/diploid) sit between lineages: never use them
MAX_ANCHOR_MIXED = 0.05
is_anchor = np.array([meta.loc[s, "role"] in ("anchor", "labelled") for s in samples])
for i in np.where(is_anchor & (het_frac > MAX_ANCHOR_MIXED))[0]:
    meta.loc[samples[i], "role"] = "demoted"
is_anchor &= het_frac <= MAX_ANCHOR_MIXED
anchors = np.where(is_anchor & good)[0]
rows = []
for i, s in enumerate(samples):
    # leave-one-out: an isolate is never called from itself, so labelled isolates test the method
    others = [a for a in anchors if a != i and shared[i, a] >= MIN_SHARED] or [a for a in anchors if a != i]
    a = min(others, key=lambda a: D[i, a])
    lbl = meta.loc[s, "label"]
    rows.append(dict(run=s, biosample=meta.loc[s, "biosample"], role=meta.loc[s, "role"],
                     given_label=lbl if isinstance(lbl, str) else "",
                     called=meta.loc[samples[a], "label"], nearest_anchor_snps=round(float(D[i, a])),
                     shared_sites=int(shared[i, a]),
                     map_rate=round(float(map_rate[i]), 3), missing=round(float(miss_s[i]), 3),
                     good=bool(good[i]), mixed_allele_frac=round(float(het_frac[i]), 4),
                     bioproject=meta.loc[s, "bioproject"], strain=meta.loc[s, "strain"],
                     country=meta.loc[s, "country"]))
res = pd.DataFrame(rows)
res["confidence"] = np.where(~res.good | (res.shared_sites < MIN_SHARED), "low", "ok")

# clonal groups among good samples only (sparse samples have too few shared sites)
gi = np.where(good)[0]
cl = fcluster(average(squareform(D[np.ix_(gi, gi)], checks=False)), CLONAL_SNPS, criterion="distance")
res["clonal_group"] = 0
res.loc[gi, "clonal_group"] = cl
sizes = res[res.clonal_group > 0].clonal_group.value_counts()
res["clonal_group_size"] = res.clonal_group.map(sizes).fillna(0).astype(int)

res.to_csv("pilot_lineages.tsv", sep="\t", index=False)
pd.DataFrame(D, index=samples, columns=samples).round(0).to_csv("pilot_snp_distances.tsv", sep="\t")

pd.set_option("display.width", 200)
print("\nSample QC (sorted by missing):")
print(res.sort_values("missing", ascending=False)[["run", "role", "map_rate", "missing", "good", "bioproject"]]
      .head(12).to_string(index=False))
print("\nAnchor check (anchor labelled by its nearest *other* good anchor):")
print(res[res.role.isin(["anchor", "demoted"])][["run", "role", "given_label", "called", "nearest_anchor_snps", "shared_sites",
                                  "confidence"]].to_string(index=False))
print("\nAnchor-to-anchor distances (good anchors, scaled SNPs):")
an = [samples[a] for a in anchors]
lab = [f"{meta.loc[s, 'label']}:{s[-4:]}" for s in an]
print(pd.DataFrame(D[np.ix_(anchors, anchors)], index=lab, columns=lab).round(0).astype(int).to_string())
print("\nLineage calls for pilot isolates (confidence ok / low):")
print(pd.crosstab(res[res.role == "pilot"].called, res[res.role == "pilot"].confidence).to_string())
print("\nPossible hybrids/diploids (mixed_allele_frac > 0.05):")
print(res[res.mixed_allele_frac > 0.05][["run", "strain", "country", "given_label", "called", "mixed_allele_frac", "missing"]]
      .to_string(index=False) or "none")
print(f"\nClonal groups (< {CLONAL_SNPS} SNPs) with >1 member:")
cg = res[res.clonal_group_size > 1].sort_values("clonal_group")
print(cg[["clonal_group", "run", "strain", "bioproject", "country"]].to_string(index=False) if len(cg) else "none")

lab = res[res.role.isin(["labelled", "demoted", "conflict"])].copy()
if len(lab):
    major = lambda l: "VNI" if str(l).startswith("VNI") and not str(l).startswith("VNII") else str(l)
    lab["agree_major"] = lab.given_label.map(major) == lab.called.map(major)
    lab["agree_exact"] = lab.given_label == lab.called
    print("\nValidation: published label vs leave-one-out call (labelled isolates):")
    v = lab[lab.role == "labelled"]
    print(v.groupby("given_label").agg(n=("run", "size"), major_ok=("agree_major", "sum"), exact_ok=("agree_exact", "sum"),
                                       mixed=("mixed_allele_frac", lambda x: (x > MAX_ANCHOR_MIXED).sum()),
                                       low_conf=("confidence", lambda x: (x == "low").sum())).to_string())
    bad = lab[~lab.agree_major | (lab.mixed_allele_frac > MAX_ANCHOR_MIXED) | (lab.role != "labelled")]
    print("\nDisagreements, mixed, demoted and conflict isolates:")
    print(bad[["run", "role", "strain", "country", "given_label", "called", "nearest_anchor_snps", "mixed_allele_frac",
               "missing", "map_rate"]].to_string(index=False) if len(bad) else "none")

