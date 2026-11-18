"""Leave-one-out sensitivity and lineage scores per isolate.

Markers per target T: M = groups/T.ed2 (current markers), Mx = groups/T.x.ed2 (k-mers in exactly
need-1 target isolates that pass the same filters; loo_extra.sh). Holding out a target isolate X
leaves the out-group unchanged and lowers the threshold by one, so the markers built without X are
M + (Mx - X), and X's leave-one-out sensitivity is |M & X| / (|M| + |Mx - X|). Exact, no rebuild.

Training isolates outside T score 0 against M by construction (markers are absent from every
training out-group isolate), so specificity needs isolates never used in training. The non-training
runs in samples.tsv (low coverage, mixed, conflicting labels) are scored against M as unseen cases.

Presence = the isolate's KMC presence database (error-filtered solid k-mers), so low coverage reads
as absence. Writes loo_sensitivity.tsv (training isolate x own target) and scores.tsv (every run x
every target; LOO markers for its own targets).

Usage: micromamba run -n cladeprint python score_isolates.py
"""
import os
import subprocess
import numpy as np
import pandas as pd

K = 31
TARGETS = ["VNII", "VNBI", "VNBII", "VNI", "VNIa", "VNIb", "VNIc"]
LUT = np.full(256, 255, np.uint8)
for i, b in enumerate(b"ACGT"):
    LUT[b] = i

os.chdir(os.path.dirname(os.path.abspath(__file__)))


def decode(v):
    codes = np.stack([((v >> np.uint64(2 * (K - 1 - i))) & np.uint64(3)).astype(np.uint8) for i in range(K)], 1)
    return np.frombuffer(b"ACGT", np.uint8)[codes].view(f"S{K}").ravel()


def encode(seqs):
    codes = LUT[np.asarray(seqs, dtype=f"S{K}").view(np.uint8).reshape(-1, K)]
    v = np.zeros(len(codes), np.uint64)
    for i in range(K):
        v = (v << np.uint64(2)) | codes[:, i].astype(np.uint64)
    return v


def kmc(*args):
    subprocess.run(args, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


M = {t: np.load(f"groups/{t}.ed2.u64.npy") for t in TARGETS}
MX = {t: np.load(f"groups/{t}.x.ed2.u64.npy") for t in TARGETS}
allm = np.unique(np.concatenate(list(M.values()) + list(MX.values())))

# One KMC database of every marker, so each isolate needs a single small intersection.
mdb = "groups/markers"
if not os.path.exists(mdb + ".kmc_suf"):
    fa = "tmp/markers.fa"
    with open(fa, "wb") as fh:
        for s in decode(allm):
            fh.write(b">m\n" + s + b"\n")
    kmc("kmc", "-k31", "-fm", "-ci1", "-cs255", "-t2", "-m2", fa, mdb, "tmp/")
    os.remove(fa)

samples = pd.read_csv("samples.tsv", sep="\t").set_index("run")
cache = "tmp/present"
os.makedirs(cache, exist_ok=True)


def present(run):
    """Sorted marker k-mers present in the run's presence database."""
    npy = f"{cache}/{run}.npy"
    if os.path.exists(npy):
        return np.load(npy)
    out = f"tmp/{run}.mk"
    kmc("kmc_tools", "-t2", "simple", mdb, f"db/{run}", "intersect", out)
    kmc("kmc_tools", "-t2", "transform", out, "dump", out + ".txt")
    txt = open(out + ".txt").read().split()
    v = np.sort(encode(txt[0::2])) if txt else np.zeros(0, np.uint64)
    for ext in (".kmc_pre", ".kmc_suf", ".txt"):
        os.remove(out + ext)
    np.save(npy, v)
    return v


def member(s, t):
    return s.major == t or s.called == t


loo, scores = [], []
for run, s in samples.iterrows():
    p = present(run)
    row = dict(run=run, called=s.called, use=s.use, peak_cov=s.peak_cov)
    for t in TARGETS:
        m_in = np.isin(M[t], p, assume_unique=True).sum()
        if s.use == "train" and member(s, t):
            x_in = np.isin(MX[t], p, assume_unique=True).sum()
            denom = len(M[t]) + len(MX[t]) - x_in
            row[t] = m_in / denom
            loo.append(dict(run=run, target=t, called=s.called, peak_cov=s.peak_cov, markers=denom,
                            present=int(m_in), sensitivity=round(m_in / denom, 4)))
        else:
            row[t] = m_in / len(M[t])
    scores.append(row)

loo = pd.DataFrame(loo)
loo.to_csv("loo_sensitivity.tsv", sep="\t", index=False)
scores = pd.DataFrame(scores).round(4)
scores.to_csv("scores.tsv", sep="\t", index=False)

print("Leave-one-out sensitivity (fraction of the markers built without the isolate that it carries):")
print(loo.groupby("target").sensitivity.describe(percentiles=[0.1, 0.5])[["count", "min", "10%", "50%", "max"]]
      .loc[TARGETS].round(3).to_string())
print("\nLowest 10:")
print(loo.nsmallest(10, "sensitivity").to_string(index=False))
print("\nUnseen runs (not used in training), score per target:")
print(scores[scores.use != "train"].to_string(index=False))
