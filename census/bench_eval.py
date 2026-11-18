"""R4 evaluation: every tool's markers measured with the same yardstick.

  python census/bench_eval.py kmers <bench>      (on the instance: needs kmc_tools, the tool outputs)
      -> <bench>/eval/kmers.txt (sampled marker 31-mers), index.tsv (tool, target, start, end, n_total),
         genomes.txt (all benchmark runs)
  examples/marker_dist eval/kmers.txt eval/genomes.txt <db> eval/dist.u8
  python census/bench_eval.py evaluate <bench>
      -> <bench>/eval/summary.tsv

Per tool and target (markers sampled to at most N_SAMPLE k-mers, fixed seed):
  in_test_med / in_test_min   share of markers present in each test in-group genome
  out_test_max                largest share present in a test out-group genome
  hit0_train / hit0_test      share of markers present (exactly) in >= 1 out-group genome
  near1_* / near2_*           share within 1 / 2 (Hamming or edit) of >= 1 out-group genome
"train" out-group genomes are those the tool saw; "test" ones it never saw.
"""

import random
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

N_SAMPLE = 2000
SEED = 7
K = 31
RC = str.maketrans("ACGT", "TGCA")


def kmers_of(seq):
    seq = seq.upper()
    return [seq[i : i + K] for i in range(len(seq) - K + 1) if set(seq[i : i + K]) <= set("ACGT")]


def read_fasta(path):
    seqs, cur = [], []
    for line in open(path):
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            cur = []
        else:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


def kmc_dump(prefix, tmp):
    subprocess.run(["kmc_tools", "transform", str(prefix), "dump", str(tmp)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return [l.split()[0] for l in open(tmp)]


def kover_kmers(res):
    """Rule k-mers of a Kover model (model.fasta), with equivalent rules when written."""
    out = []
    for f in ["model.fasta"] + sorted(p.name for p in res.glob("model_rule_*_equiv.fa*")):
        p = res / f
        if p.exists():
            out += [s for s in read_fasta(p) if len(s) == K]
    return out


def kmers(bench):
    bench = Path(bench)
    ev = bench / "eval"
    ev.mkdir(exist_ok=True)
    t = pd.read_csv(bench / "targets.tsv", sep="\t")
    rng = random.Random(SEED)
    rows, allk = [], []
    for _, tr in t.iterrows():
        f = tr.target.split("/")
        fname = "__".join([f[0]] + f[1:]) if len(f) > 1 else f[0]
        for tool in ["cp", "cpx", "neptune", "kover"]:
            if tool in ("cp", "cpx"):
                pre = bench / tool / "bundle/markers" / tr.level / fname
                ks = kmc_dump(pre, ev / "dump.tmp") if Path(str(pre) + ".kmc_pre").exists() else []
            elif tool == "neptune":
                fa = bench / f"neptune/out_{tr.tag}/consolidated/consolidated.fasta"
                ks = [k for s in read_fasta(fa) for k in kmers_of(s)] if fa.exists() else []
            else:
                ks = kover_kmers(bench / f"kover/res_{tr.tag}")
            ks = sorted(set(ks))
            n_total = len(ks)
            # Kover: the model's own rule k-mers always come first (the rest are its equivalents)
            rule = []
            if tool == "kover":
                mf = bench / f"kover/res_{tr.tag}/model.fasta"
                rule = sorted({s for s in read_fasta(mf) if len(s) == K}) if mf.exists() else []
                ks = [k for k in ks if k not in set(rule)]
            if len(ks) > N_SAMPLE - len(rule):
                ks = rng.sample(ks, N_SAMPLE - len(rule))
            ks = rule + ks
            rows.append((tool, tr.target, len(allk), len(allk) + len(ks), n_total))
            allk += ks
    (ev / "dump.tmp").unlink(missing_ok=True)
    (ev / "kmers.txt").write_text("\n".join(allk) + "\n")
    pd.DataFrame(rows, columns=["tool", "target", "start", "end", "n_total"]).to_csv(
        ev / "index.tsv", sep="\t", index=False
    )
    b = pd.read_csv(bench / "samples.tsv", sep="\t", dtype=str)
    (ev / "genomes.txt").write_text("\n".join(b["sample"]) + "\n")
    print(pd.DataFrame(rows, columns=["tool", "target", "start", "end", "n_total"]).to_string(index=False))


def evaluate(bench):
    bench = Path(bench)
    ev = bench / "eval"
    idx = pd.read_csv(ev / "index.tsv", sep="\t")
    genomes = (ev / "genomes.txt").read_text().split()
    d = np.fromfile(ev / "dist.u8", dtype=np.uint8).reshape(len(genomes), -1)
    b = pd.read_csv(bench / "samples.tsv", sep="\t", dtype=str).set_index("sample").loc[genomes]
    levels = ["species", "lineage", "sublineage"]
    t = pd.read_csv(bench / "targets.tsv", sep="\t").set_index("target")

    def member(target):
        lv = t.loc[target, "level"]
        cols = levels[: levels.index(lv) + 1]
        p = b[cols].apply(lambda r: None if r.isna().any() else "/".join(r), axis=1)
        return (p == target).values

    train = (b.split == "train").values
    rows = []
    for _, r in idx.iterrows():
        m = d[:, r.start : r.end]
        n = m.shape[1]
        row = dict(tool=r.tool, target=r.target, n_total=r.n_total, n_eval=n)
        if n:
            ing = member(r.target)
            pres = (m == 0).mean(axis=1)
            row.update(
                in_test_med=np.median(pres[ing & ~train]),
                in_test_min=pres[ing & ~train].min(),
                out_test_max=pres[~ing & ~train].max(),
            )
            for name, sel in [("train", ~ing & train), ("test", ~ing & ~train)]:
                mo = m[sel].min(axis=0)
                row[f"hit0_{name}"] = (mo == 0).mean()
                row[f"near1_{name}"] = (mo <= 1).mean()
                row[f"near2_{name}"] = (mo <= 2).mean()
        rows.append(row)
    s = pd.DataFrame(rows)
    s.to_csv(ev / "summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(s.round(3).to_string(index=False))


if __name__ == "__main__":
    {"kmers": kmers, "evaluate": evaluate}[sys.argv[1]](*sys.argv[2:])
