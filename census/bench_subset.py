"""Benchmark subset for the tool comparison (R4: cladeprint vs exact k-mers vs Neptune vs Kover).

Neptune and Kover take genomes, not presence databases, and are built for tens to hundreds of
bacterial genomes, so the comparison runs on a subset; every tool sees the same k-mers (each
selected run's presence database as unitigs, bench_genomes.sh).

Selection: literature-labelled runs that pass QC in the v2 census run, one per clone cluster
(cluster_5), at most PER_LEAF per leaf group (the six non-neoformans species; VNII; the five VNI/VNB
sublineages), drawn with a fixed seed; within each leaf, every third run (after shuffling) is test,
the rest train. Clone clusters never straddle the split.

Writes <out>/bench/samples.tsv (sample, leaf, species, lineage, sublineage, split).
Usage: python census/bench_subset.py census/out
"""

import sys
from pathlib import Path

import pandas as pd

PER_LEAF = 18
SEED = 20261018


def main(out):
    out = Path(out)
    s = pd.read_csv(out / "samples.tsv", sep="\t", dtype=str)
    qc = pd.read_csv(out / "qc.tsv", sep="\t", dtype=str)[["sample", "train"]]
    s = s.merge(qc.rename(columns={"train": "qc_train"}), on="sample")
    s = s[s.qc_train == "true"]
    neo = s.species == "C.neoformans"
    s = s[(~neo & s.species.notna()) | (neo & (s.sublineage.notna() | (s.lineage == "VNII")))]
    s["leaf"] = s.species.where(~(s.species == "C.neoformans"), s.lineage.where(s.lineage == "VNII", s.sublineage))
    s = s.sample(frac=1, random_state=SEED).drop_duplicates("cluster")
    picks = []
    for leaf, d in s.groupby("leaf"):
        d = d.head(PER_LEAF).copy()
        d["split"] = ["test" if i % 3 == 2 else "train" for i in range(len(d))]
        picks.append(d)
    b = pd.concat(picks)
    (out / "bench").mkdir(exist_ok=True)
    b[["sample", "leaf", "species", "lineage", "sublineage", "split"]].to_csv(
        out / "bench/samples.tsv", sep="\t", index=False
    )
    print(b.groupby(["leaf", "split"]).size().unstack().to_string())
    print("total", len(b))


if __name__ == "__main__":
    main(sys.argv[1])
