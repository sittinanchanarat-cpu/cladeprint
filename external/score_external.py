"""Score the external test isolates against the lineage markers (groups/T.ed2).

Score(isolate, T) = share of T's markers present in the isolate's presence database. None of these
isolates was used in training, so off-target scores measure specificity directly.

Call: the major lineage (VNII, VNBI, VNBII, VNI) with the highest score if it reaches CALL; for VNI,
also the best of VNIa/b/c. A second major lineage reaching CALL is reported too (hybrid/mixed).
CALL = 0.3: training leave-one-out sensitivity never fell below 0.48, while off-target scores in the
unseen training-excluded runs stayed <= 0.012.

Writes external/scores.tsv and prints per-class summaries.
Usage: micromamba run -n cladeprint python external/score_external.py
"""
import os
import subprocess
import numpy as np
import pandas as pd

K = 31
MAJOR = ["VNII", "VNBI", "VNBII", "VNI"]
SUB = ["VNIa", "VNIb", "VNIc"]
CALL = 0.3
LUT = np.full(256, 255, np.uint8)
for i, b in enumerate(b"ACGT"):
    LUT[b] = i

os.chdir(os.path.dirname(os.path.abspath(__file__)))
MDB = "../kmers/groups/markers"  # every marker (and LOO extra) k-mer; built by score_isolates.py


def encode(seqs):
    codes = LUT[np.asarray(seqs, dtype=f"S{K}").view(np.uint8).reshape(-1, K)]
    v = np.zeros(len(codes), np.uint64)
    for i in range(K):
        v = (v << np.uint64(2)) | codes[:, i].astype(np.uint64)
    return v


def present(run):
    out = f"tmp/{run}.mk"
    for args in (["simple", MDB, f"db/{run}", "intersect", out], ["transform", out, "dump", out + ".txt"]):
        subprocess.run(["kmc_tools", "-t2", *args], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    txt = open(out + ".txt").read().split()
    for ext in (".kmc_pre", ".kmc_suf", ".txt"):
        os.remove(out + ext)
    return np.sort(encode(txt[0::2])) if txt else np.zeros(0, np.uint64)


M = {t: np.load(f"../kmers/groups/{t}.ed2.u64.npy") for t in MAJOR + SUB}
samples = pd.read_csv("samples.tsv", sep="\t")
stats = pd.read_csv("stats.tsv", sep="\t").set_index("run")
rows = []
for _, s in samples.iterrows():
    if not os.path.exists(f"db/{s.run}.kmc_suf"):
        continue
    p = present(s.run)
    sc = {t: np.isin(M[t], p, assume_unique=True).mean() for t in M}
    major = sorted(MAJOR, key=sc.get, reverse=True)
    call = major[0] if sc[major[0]] >= CALL else "none"
    if call == "VNI":
        sub = max(SUB, key=sc.get)
        call = sub if sc[sub] >= CALL else "VNI"
    second = major[1] if sc[major[1]] >= CALL else ""
    rows.append(dict(isolate=s.isolate, label=s.lineage, run=s.run, country=s.country,
                     peak_cov=stats.peak_cov.get(s.run), **{t: round(sc[t], 4) for t in M},
                     call=call, second=second))

res = pd.DataFrame(rows)
res["low_cov"] = res.run.map(stats.low_cov)
res.to_csv("scores.tsv", sep="\t", index=False)
# metrics on isolates with a usable genome peak; low-coverage ones are reported separately
lc = res[(res.low_cov == 1) & ~res.label.str.contains("/")]
lc_top = sum(max(MAJOR, key=r.get) == {"VNIa": "VNI", "VNIb": "VNI", "VNIc": "VNI"}.get(r.label, r.label)
             for _, r in lc.iterrows())
print(f"Low-coverage pure isolates (no metrics): {len(lc)}; top major-lineage score matches the label in {lc_top}")
pure = res[~res.label.str.contains("/") & (res.low_cov == 0)]
ok = pure.call == pure.label
print(f"Pure-lineage isolates: {ok.sum()}/{len(pure)} called as labelled")
off = [max(r[t] for t in MAJOR if t != {"VNIa": "VNI", "VNIb": "VNI", "VNIc": "VNI"}.get(r.label, r.label))
       for _, r in pure.iterrows()]
print(f"Max off-target major-lineage score: {max(off, default=0):.4f}")
print("\nOwn-lineage scores (external sensitivity):")
own = pure.apply(lambda r: r[r.label], axis=1)
print(own.groupby(pure.label).describe()[["count", "min", "50%", "max"]].round(3).to_string())
print("\nMismatches and hybrids:")
cols = ["isolate", "label", "peak_cov", *MAJOR, *SUB, "call", "second"]
print(pd.concat([pure[~ok], res[res.label.str.contains("/")]])[cols + ["low_cov"]].to_string(index=False))
