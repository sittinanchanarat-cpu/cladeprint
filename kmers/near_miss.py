"""Near-miss exclusion (design step 4): drop candidate markers that any out-group isolate holds
with <= 2 mismatches.

Input per target T (from group_kmers.py): groups/T.strict (candidates absent from the out-group)
and groups/T.out (every k-mer held by at least one out-group isolate). Both are canonical
31-mers. A genomic near-neighbour y of candidate c is stored as y or rc(y), so we query both
c and rc(c) against the canonical out-group set.

Search: pigeonhole. Split 31-mers into blocks of 10/10/11 bases; within 2 mismatches, at least
one block matches exactly. For each block, sort the out-group by that block and pair each query
with out-group k-mers sharing it; Hamming distance on 2-bit codes via xor + popcount. Memory is
a few arrays of ~55M uint64 (~1.5 GB), independent of the number of candidates.

--verify N: brute-force check of N random candidates by enumerating every <= 2-mismatch variant
(~4,300 per orientation) and looking each one up in the out-group. Exits non-zero on mismatch.

Usage: micromamba run -n cladeprint python near_miss.py [--verify N] TARGET [TARGET ...]
"""
import argparse
import os
import subprocess
import time
import numpy as np
import pandas as pd

K = 31
MASK = (1 << 2 * K) - 1
LUT = np.full(256, 255, np.uint8)
for i, b in enumerate(b"ACGT"):
    LUT[b] = i
BLOCKS = [(42, 10), (22, 10), (0, 11)]  # (bit shift, bases): bases 0-9, 10-19, 20-30
MAX_PAIRS = 10_000_000  # ~1 GB of temporaries per chunk; 40M thrashed swap on VNII

os.chdir(os.path.dirname(os.path.abspath(__file__)))


def load(db):
    """Canonical k-mers of a KMC database as sorted uint64 (2 bits/base), cached as .npy."""
    cache = db + ".u64.npy"
    if os.path.exists(cache):
        return np.load(cache)
    txt = db + ".dump.txt"
    subprocess.run(["kmc_tools", "-t2", "transform", db, "dump", txt], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    parts = []
    for ch in pd.read_csv(txt, sep="\t", header=None, usecols=[0], dtype=str, chunksize=4_000_000):
        codes = LUT[ch[0].to_numpy().astype(f"S{K}").view(np.uint8).reshape(-1, K)]
        assert (codes < 4).all(), "non-ACGT k-mer in dump"
        v = np.zeros(len(codes), np.uint64)
        for i in range(K):
            v = (v << np.uint64(2)) | codes[:, i].astype(np.uint64)
        parts.append(v)
    os.remove(txt)
    v = np.sort(np.concatenate(parts)) if parts else np.zeros(0, np.uint64)
    np.save(cache, v)
    return v


def revcomp(v):
    x = v ^ np.uint64(MASK)
    r = np.zeros_like(x)
    for _ in range(K):
        r = (r << np.uint64(2)) | (x & np.uint64(3))
        x = x >> np.uint64(2)
    return r


def hamming(a, b):
    x = a ^ b
    return np.bitwise_count((x | (x >> np.uint64(1))) & np.uint64(0x5555555555555555 & MASK))


def min_dist(queries, out):
    """Min Hamming distance (0-2) from each query to the out-group; 3 means > 2."""
    best = np.full(len(queries), 3, np.uint8)
    for shift, n in BLOCKS:
        m = np.uint64((1 << 2 * n) - 1)
        okey = ((out >> np.uint64(shift)) & m).astype(np.uint32)
        order = np.argsort(okey, kind="stable")
        okey = okey[order]
        qkey = ((queries >> np.uint64(shift)) & m).astype(np.uint32)
        lo = np.searchsorted(okey, qkey, "left")
        cnt = np.searchsorted(okey, qkey, "right") - lo
        start = 0
        while start < len(queries):  # chunk so the expanded pair list stays bounded
            cum = np.cumsum(cnt[start:])
            stop = start + max(1, int(np.searchsorted(cum, MAX_PAIRS, "right")))
            c = cnt[start:stop]
            qi = np.repeat(np.arange(start, stop), c)
            off = np.arange(len(qi)) - np.repeat(np.cumsum(c) - c, c)
            oi = order[np.repeat(lo[start:stop], c) + off]
            d = hamming(queries[qi], out[oi])
            hit = d < 3
            np.minimum.at(best, qi[hit], d[hit].astype(np.uint8))
            start = stop
        del okey, order
    return best


def brute(c, out):
    """Min distance by enumerating all <= 2-mismatch variants of c and rc(c)."""
    best = 3
    for q in (int(c), int(revcomp(np.array([c], np.uint64))[0])):
        var = [q]
        for i in range(K):
            s = 2 * (K - 1 - i)
            for a in range(1, 4):
                v1 = q ^ (a << s)
                var.append(v1)
                for j in range(i + 1, K):
                    t = 2 * (K - 1 - j)
                    for b in range(1, 4):
                        var.append(v1 ^ (b << t))
        var = np.array(var, np.uint64)
        canon = np.minimum(var, revcomp(var))
        pos = np.minimum(np.searchsorted(out, canon), len(out) - 1)
        found = out[pos] == canon
        if found.any():
            d = hamming(np.full(found.sum(), q, np.uint64), var[found])
            best = min(best, int(d.min()))
    return best


def save_row(row, path="near_miss.tsv"):
    """Upsert one target's row, so finished targets survive a later crash."""
    res = pd.DataFrame([row])
    if os.path.exists(path):
        res = pd.concat([pd.read_csv(path, sep="\t").query("target != @row['target']"), res])
    res.to_csv(path, sep="\t", index=False)


ap = argparse.ArgumentParser()
ap.add_argument("targets", nargs="+")
ap.add_argument("--verify", type=int, default=0)
args = ap.parse_args()

rows = []
for t in args.targets:
    t0 = time.time()
    cand = load(f"groups/{t}.strict")
    out = load(f"groups/{t}.out")
    q = np.concatenate([cand, revcomp(cand)])
    d = min_dist(q, out)
    d = np.minimum(d[:len(cand)], d[len(cand):])
    keep = cand[d == 3]
    np.save(f"groups/{t}.nm2.u64.npy", keep)
    rows.append(dict(target=t, strict=len(cand), out_kmers=len(out), d1=int((d == 1).sum()),
                     d2=int((d == 2).sum()), d0=int((d == 0).sum()), survive=len(keep),
                     survive_frac=round(len(keep) / max(1, len(cand)), 4), secs=round(time.time() - t0)))
    print(rows[-1], flush=True)
    save_row(rows[-1])
    if args.verify:
        rng = np.random.default_rng(0)
        # half random, half from those flagged as near misses, so both outcomes are tested
        pick = np.concatenate([rng.choice(len(cand), min(args.verify // 2, len(cand)), replace=False),
                               rng.permutation(np.flatnonzero(d < 3))[:args.verify // 2]])
        checked = [(int(i), int(d[i]), brute(cand[i], out)) for i in pick]
        bad = [x for x in checked if x[1] != x[2]]
        print(f"  verify {t}: {len(pick)} checked, {len(bad)} disagree {bad[:5]}", flush=True)
        if bad:
            raise SystemExit(1)

print(pd.DataFrame(rows).to_string(index=False))
