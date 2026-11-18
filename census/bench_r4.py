"""R4 tool comparison on the benchmark subset (census/bench_subset.py): inputs for each tool.

Targets: every group cladeprint builds on the subset (7 species, 3 lineages, 5 sublineages), each as
one-vs-rest. Training genomes only for discovery; test genomes only for classification. Every
tool gets the same genomes: the runs' presence databases (cladeprint) or the unitigs made from
them (Neptune, Kover), which hold exactly the same 31-mers.

  python census/bench_r4.py prepare census/out/bench
writes, under bench/:
  cp/samples.tsv, cp/config.toml        cladeprint (near-miss on)
  cpx/samples.tsv, cpx/config.toml      exact-k-mer baseline (near_miss = false)
  test.tsv                              test genomes (for classify)
  targets.tsv                           target, level, members (train), members (test)
  kover/genomes.tsv                     id <tab> unitig FASTA (path inside the container)
  kover/<target>.pheno, train.ids, test.ids
  neptune/<target>.inc, .exc            inclusion / exclusion FASTA lists (train only)
"""

import sys
from pathlib import Path

import pandas as pd

CONFIG = """samples = "samples.tsv"
output = "."
levels = ["species", "lineage", "sublineage"]

[markers]
max_out = 1
near_miss = {near_miss}

[resources]
threads = 8
memory_gb = 8
"""


def targets(b):
    """(name, level, mask over b) for every group with >= 3 training genomes."""
    out = []
    for level, cols in [
        ("species", ["species"]),
        ("lineage", ["species", "lineage"]),
        ("sublineage", ["species", "lineage", "sublineage"]),
    ]:
        path = b[cols].apply(lambda r: None if r.isna().any() else "/".join(r), axis=1)
        for g in sorted(path.dropna().unique()):
            m = path == g
            if (m & (b.split == "train")).sum() >= 3:
                out.append((g, level, m))
    return out


def prepare(bench):
    bench = Path(bench)
    b = pd.read_csv(bench / "samples.tsv", sep="\t", dtype=str)
    sheet = b[["sample", "species", "lineage", "sublineage"]].copy()
    sheet.insert(1, "reads", "fetched")
    sheet["train"] = b.split.map({"train": "yes", "test": "no"})
    for d, nm in [("cp", "true"), ("cpx", "false")]:
        (bench / d).mkdir(exist_ok=True)
        sheet.to_csv(bench / d / "samples.tsv", sep="\t", index=False)
        (bench / d / "config.toml").write_text(CONFIG.format(near_miss=nm))
    sheet[b.split == "test"].to_csv(bench / "test.tsv", sep="\t", index=False)

    (bench / "kover").mkdir(exist_ok=True)
    (bench / "neptune").mkdir(exist_ok=True)
    pd.DataFrame({"id": b["sample"], "path": "/data/genomes/" + b["sample"] + ".fa"}).to_csv(
        bench / "kover/genomes.tsv", sep="\t", index=False, header=False
    )
    b.loc[b.split == "train", "sample"].to_csv(bench / "kover/train.ids", index=False, header=False)
    b.loc[b.split == "test", "sample"].to_csv(bench / "kover/test.ids", index=False, header=False)
    rows = []
    train = b.split == "train"
    for name, level, m in targets(b):
        tag = name.replace("/", "_")
        pd.DataFrame({"id": b["sample"], "y": m.astype(int)}).to_csv(
            bench / f"kover/{tag}.pheno", sep="\t", index=False, header=False
        )
        fa = "genomes/" + b["sample"] + ".fa"
        fa[train & m].to_csv(bench / f"neptune/{tag}.inc", index=False, header=False)
        fa[train & ~m].to_csv(bench / f"neptune/{tag}.exc", index=False, header=False)
        rows.append((name, tag, level, int((train & m).sum()), int((~train & m).sum())))
    t = pd.DataFrame(rows, columns=["target", "tag", "level", "n_train", "n_test"])
    t.to_csv(bench / "targets.tsv", sep="\t", index=False)
    print(t.to_string(index=False))


if __name__ == "__main__":
    {"prepare": prepare}[sys.argv[1]](*sys.argv[2:])
