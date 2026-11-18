"""Edit-distance near-miss exclusion: drop Hamming survivors (groups/T.nm2) that an out-group
k-mer matches within 2 edits once indels are allowed.

The Hamming filter (near_miss.py) misses indels: a 1-bp indel shifts the rest of a 31-mer, so it
looks far in Hamming terms. Here a candidate c is rejected if some out-group 31-mer y aligns to it
with <= 2 substitutions/indels, allowing up to 2 unaligned bases at either end of either k-mer
(overlap alignment within diagonals -2..2). The free ends are needed because y is a fixed 31-base
window: after an indel, the matching out-group sequence runs 1-2 bases past y or stops short. So a
candidate is rejected if an out-group k-mer overlaps it by >= 29 bases with <= 2 edits. This is
slightly conservative: a k-mer whose only differences sit in its last 2 bases is also rejected.

Search: pigeonhole. Two edits break at most two of the three 9-base blocks c[2:11], c[11:20],
c[20:29] (the 27-base core that always aligns), so one block matches y exactly, shifted by -2..2.
For each block and shift, sort the out-group by that 9-mer, pair each query with out-group k-mers
sharing it, and score the pairs with a banded DP vectorised over pairs. Queries are c and rc(c)
against the canonical out-group set, as in near_miss.py.

Prefilter: in any alignment with <= 2 edits, at least 25 of the 27 core bases of c equal y at some
shift in -2..2 (each matched base is equal at its own diagonal). That count takes a few 64-bit ops
per pair and passes only a few % of random pairs, so the DP runs on far fewer pairs.

--test: check the vectorised DP against a plain Python DP on random and edited pairs.
--verify N: brute-force N candidates (half kept, half rejected) against every out-group k-mer.

Usage: micromamba run -n cladeprint python edit_miss.py [--verify N] [--test] TARGET [TARGET ...]
"""
import argparse
import os
import time
import numpy as np
import pandas as pd

K = 31
MASK = (1 << 2 * K) - 1
BLOCKS = [2, 11, 20]  # start of each 9-base block in the candidate
BLEN = 9
SHIFTS = range(-2, 3)
FREE = 2  # unaligned bases allowed at each end
MAX_PAIRS = 10_000_000  # candidate pairs per chunk (~0.8 GB of temporaries)
FAR = 3  # distances are capped here; 3 means "> 2"
EVEN = np.uint64(0x5555555555555555 & MASK)
CORE = np.uint64(sum(1 << 2 * (K - 1 - i) for i in range(BLOCKS[0], BLOCKS[-1] + BLEN)))
CORE_MIN = BLOCKS[-1] + BLEN - BLOCKS[0] - 2  # 27 core bases, at most 2 lost to edits

os.chdir(os.path.dirname(os.path.abspath(__file__)))


def revcomp(v):
    x = v ^ np.uint64(MASK)
    r = np.zeros_like(x)
    for _ in range(K):
        r = (r << np.uint64(2)) | (x & np.uint64(3))
        x = x >> np.uint64(2)
    return r


def bases(v):
    """(K, n) uint8 matrix of bases, row i = base i (most significant first)."""
    return np.stack([((v >> np.uint64(2 * (K - 1 - i))) & np.uint64(3)).astype(np.uint8) for i in range(K)])


def overlap_ed(a, b):
    """Banded overlap edit distance (capped at FAR) between 31-mers a[i] and b[i].

    D[i][j] = edits aligning a[:i] with b[:j]; diagonal d = j - i kept in -2..2. Either sequence may
    skip up to FREE leading bases for free (D[i][0] = 0 for i <= FREE, D[0][j] = 0 for j <= FREE),
    and the alignment may end at D[K][j] or D[i][K] for i, j >= K - FREE.
    """
    A, B = bases(a), bases(b)
    n = len(a)
    far = np.full(n, FAR, np.uint8)
    zero = np.zeros(n, np.uint8)
    prev = [far if d < 0 else zero for d in range(-2, 3)]  # row 0: D[0][d]
    best = far.copy()
    for i in range(1, K + 1):
        cur = []
        for k, d in enumerate(range(-2, 3)):
            j = i + d
            if j < 0 or j > K:
                cur.append(far)
                continue
            if j == 0:
                cur.append(zero if i <= FREE else far)
                continue
            v = prev[k] + (A[i - 1] != B[j - 1])  # diagonal from D[i-1][j-1]
            if k + 1 <= 4:
                v = np.minimum(v, prev[k + 1] + 1)  # D[i-1][j]: gap in b
            if k - 1 >= 0:
                v = np.minimum(v, cur[k - 1] + 1)  # D[i][j-1]: gap in a
            cur.append(np.minimum(v, FAR).astype(np.uint8))
        for k, d in enumerate(range(-2, 3)):
            j = i + d
            if (i == K and j >= K - FREE) or (j == K and i >= K - FREE):
                best = np.minimum(best, cur[k])
        prev = cur
    return best


def maybe_close(a, b):
    """True where a and b could be within 2 edits (core bases equal at some shift -2..2)."""
    hit = np.zeros(len(a), np.uint64)
    for s in range(-2, 3):
        bs = (b << np.uint64(2 * s)) & np.uint64(MASK) if s > 0 else b >> np.uint64(-2 * s)
        x = a ^ bs
        hit |= ~(x | (x >> np.uint64(1))) & EVEN
    return np.bitwise_count(hit & CORE) >= CORE_MIN


def block(v, pos):
    return ((v >> np.uint64(2 * (K - pos - BLEN))) & np.uint64((1 << 2 * BLEN) - 1)).astype(np.uint32)


def min_ed(queries, out, stats):
    """Index the queries by block (a 4^9 bucket table) and stream the out-group through it, so the
    large out-group set is never sorted."""
    best = np.full(len(queries), FAR, np.uint8)
    for s in SHIFTS:
        for p in BLOCKS:
            active = np.flatnonzero(best > 0)  # queries already at distance 0 can't improve
            qk = block(queries[active], p)
            qidx = active[np.argsort(qk, kind="stable")]
            qcount = np.bincount(qk, minlength=1 << 2 * BLEN)
            qstart = np.cumsum(qcount) - qcount
            step = max(1, MAX_PAIRS * len(qcount) // max(1, 4 * len(active)))  # ~MAX_PAIRS/4 pairs expected
            for i in range(0, len(out), step):
                ok_ = block(out[i:i + step], p + s)
                c = qcount[ok_]
                sel = np.flatnonzero(c)
                c = c[sel]
                oi = np.repeat(sel + i, c)
                off = np.arange(len(oi)) - np.repeat(np.cumsum(c) - c, c)
                qi = qidx[np.repeat(qstart[ok_[sel]], c) + off]
                a, b = queries[qi], out[oi]
                close = maybe_close(a, b)
                stats["pairs"] += len(qi)
                stats["dp"] += int(close.sum())
                np.minimum.at(best, qi[close], overlap_ed(a[close], b[close]))
    return best


def brute(c, out, chunk=4_000_000):
    """Min overlap edit distance of c (both orientations) to every out-group k-mer."""
    best = FAR
    for q in (c, revcomp(np.array([c], np.uint64))[0]):
        for i in range(0, len(out), chunk):
            o = out[i:i + chunk]
            o = o[maybe_close(np.full(len(o), q, np.uint64), o)]
            if len(o):
                best = min(best, int(overlap_ed(np.full(len(o), q, np.uint64), o).min()))
    return best


def ref_ed(a, b):
    """Plain DP with the same rules, on base lists, for --test."""
    INF = 99
    D = [[INF] * (K + 1) for _ in range(K + 1)]
    for i in range(K + 1):
        for j in range(K + 1):
            if abs(j - i) > 2:
                continue
            if i == 0 or j == 0:
                D[i][j] = 0 if max(i, j) <= FREE else INF
                continue
            D[i][j] = min(D[i - 1][j - 1] + (a[i - 1] != b[j - 1]),
                          D[i - 1][j] + 1 if abs(j - i + 1) <= 2 else INF,
                          D[i][j - 1] + 1 if abs(j - 1 - i) <= 2 else INF)
    ends = [D[K][j] for j in range(K - FREE, K + 1)] + [D[i][K] for i in range(K - FREE, K + 1)]
    return min(min(ends), FAR)


def self_test(n=3000, seed=1):
    rng = np.random.default_rng(seed)
    enc = lambda s: sum(int(x) << 2 * (K - 1 - i) for i, x in enumerate(s))
    a_list, b_list = [], []
    for t in range(n):
        g = list(rng.integers(0, 4, K + 6))
        a = g[3:3 + K]
        b = g[3:]  # edit a copy of the genome window, then take a 31-base window of it
        for _ in range(rng.integers(0, 4)):
            op, p = rng.integers(0, 3), int(rng.integers(3, len(b) - 3))
            if op == 0:
                b[p] = (b[p] + int(rng.integers(1, 4))) % 4
            elif op == 1:
                del b[p]
            else:
                b.insert(p, int(rng.integers(0, 4)))
        o = int(rng.integers(0, 3))
        b = (b + [0] * K)[o:o + K]
        if t % 5 == 0:
            b = list(rng.integers(0, 4, K))  # unrelated pair
        a_list.append(a)
        b_list.append(b)
    av = np.array([enc(a) for a in a_list], np.uint64)
    bv = np.array([enc(b) for b in b_list], np.uint64)
    got = overlap_ed(av, bv)
    want = np.array([ref_ed(a, b) for a, b in zip(a_list, b_list)])
    bad = np.flatnonzero(got != want)
    print(f"self-test: {n} pairs, distance counts {np.bincount(want, minlength=4).tolist()}, "
          f"{len(bad)} disagree", flush=True)
    close = maybe_close(av, bv)
    lost = np.flatnonzero((want < FAR) & ~close)
    rnd = maybe_close(*np.random.default_rng(seed).integers(0, MASK, (2, 1_000_000), dtype=np.uint64))
    print(f"prefilter: drops {len(lost)} of {(want < FAR).sum()} true hits; passes {rnd.mean():.1%} of random pairs",
          flush=True)
    if len(bad) or len(lost):
        raise SystemExit(1)


def save_row(row, path="edit_miss.tsv"):
    res = pd.DataFrame([row])
    if os.path.exists(path):
        res = pd.concat([pd.read_csv(path, sep="\t").query("target != @row['target']"), res])
    res.to_csv(path, sep="\t", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("targets", nargs="*")
    ap.add_argument("--verify", type=int, default=0)
    ap.add_argument("--test", action="store_true")
    args = ap.parse_args()
    if args.test:
        self_test()
    for t in args.targets:
        t0 = time.time()
        cand = np.load(f"groups/{t}.nm2.u64.npy")
        out = np.load(f"groups/{t}.out.u64.npy")
        stats = dict(pairs=0, dp=0)
        d = min_ed(np.concatenate([cand, revcomp(cand)]), out, stats)
        d = np.minimum(d[:len(cand)], d[len(cand):])
        keep = cand[d == FAR]
        np.save(f"groups/{t}.ed2.u64.npy", keep)
        row = dict(target=t, hamming_survive=len(cand), ed0=int((d == 0).sum()), ed1=int((d == 1).sum()),
                   ed2=int((d == 2).sum()), survive=len(keep),
                   survive_frac=round(len(keep) / max(1, len(cand)), 4), pairs=stats["pairs"], dp_pairs=stats["dp"], secs=round(time.time() - t0))
        print(row, flush=True)
        save_row(row)
        if args.verify:
            rng = np.random.default_rng(0)
            kept, rej = np.flatnonzero(d == FAR), np.flatnonzero(d < FAR)
            pick = np.concatenate([rng.permutation(kept)[:args.verify // 2], rng.permutation(rej)[:args.verify // 2]])
            checked = [(int(i), int(d[i]), brute(cand[i], out)) for i in pick]
            bad = [x for x in checked if x[1] != x[2]]
            print(f"  verify {t}: {len(pick)} checked, {len(bad)} disagree {bad[:5]}", flush=True)
            if bad:
                raise SystemExit(1)
