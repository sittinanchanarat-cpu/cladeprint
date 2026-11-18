"""Place surviving markers (groups/T.nm2, or T.ed2 with --set ed2) on H99 (design step 5, first pass).

1. Unitigs: chain the target's canonical 31-mers through unambiguous (k-1)-overlaps, so each
   unitig is a contiguous run of marker sequence. Counts regions without any reference.
2. Map unitigs to H99 (14 chromosomes + mitochondrion) with minimap2, tuned for short queries.
3. Loci: merge placed unitigs closer than --gap bp on the same contig.

Writes groups/T[.ed2].unitigs.fa, .paf and .loci.tsv, plus one summary row per target and set in
place_markers.tsv.

Usage: micromamba run -n cladeprint python place_markers.py [--gap 500] [--set nm2|ed2] TARGET [TARGET ...]
"""
import argparse
import os
import subprocess
import time
import numpy as np
import pandas as pd

K = 31
MASK = (1 << 2 * K) - 1
TOP = 2 * (K - 1)
REF = "../pilot/ref/H99.fna"
MITO = "NC_018792.1"
END_ZONE = 50_000  # "near a chromosome end" if within this many bp

os.chdir(os.path.dirname(os.path.abspath(__file__)))


def rc_int(x):
    x ^= MASK
    r = 0
    for _ in range(K):
        r = (r << 2) | (x & 3)
        x >>= 2
    return r


def decode(v, n):
    return "".join("ACGT"[(v >> 2 * (n - 1 - i)) & 3] for i in range(n))


def unitigs(kmers):
    """Maximal non-branching paths through the node-centric de Bruijn graph of `kmers`.

    Yields (sequence, n_kmers). Tracks each oriented k-mer as (fw, rc) so extending by one base
    is O(1); canonical form = min(fw, rc).
    """
    kset = set(int(x) for x in kmers)
    seen = set()

    def succ(fw, rc):
        out = []
        for b in range(4):
            f = ((fw << 2) & MASK) | b
            r = (rc >> 2) | ((3 - b) << TOP)
            if min(f, r) in kset:
                out.append((f, r))
        return out

    def pred(fw, rc):  # predecessors of fw = successors of rc, reverse-complemented
        return [(r, f) for f, r in succ(rc, fw)]

    def extend(fw, rc):
        """Walk forward from (fw, rc); returns the appended bases."""
        bases = []
        while True:
            s = succ(fw, rc)
            if len(s) != 1:
                break
            f, r = s[0]
            c = min(f, r)
            if c in seen or len(pred(f, r)) != 1:
                break
            seen.add(c)
            bases.append(f & 3)
            fw, rc = f, r
        return bases

    for x in kset:
        if x in seen:
            continue
        seen.add(x)
        fw, rc = x, rc_int(x)
        right = extend(fw, rc)
        left = extend(rc, fw)  # extending the reverse strand = growing leftwards
        seq = "".join("ACGT"[3 - b] for b in reversed(left)) + decode(fw, K) + "".join("ACGT"[b] for b in right)
        yield seq, 1 + len(left) + len(right)


def place(target, gap, mset):
    t0 = time.time()
    km = np.load(f"groups/{target}.{mset}.u64.npy")
    stem = f"groups/{target}" + ("" if mset == "nm2" else f".{mset}")
    fa = f"{stem}.unitigs.fa"
    lens = []
    with open(fa, "w") as fh:
        for i, (seq, n) in enumerate(unitigs(km)):
            fh.write(f">u{i}_n{n}\n{seq}\n")
            lens.append(len(seq))
    lens = np.array(lens)
    assert (lens - K + 1).sum() == len(km), "every k-mer must land in exactly one unitig"

    # Short queries with several mismatches against H99: small seeds, dense minimizers.
    paf = f"{stem}.paf"
    subprocess.run(["minimap2", "-c", "-x", "sr", "-k", "11", "-w", "3", "--secondary=no", "-t", "2",
                    "-o", paf, REF, fa], check=True, stderr=subprocess.DEVNULL)
    cols = ["q", "qlen", "qs", "qe", "strand", "t", "tlen", "ts", "te", "match", "alen", "mapq"]
    hits = pd.read_csv(paf, sep="\t", header=None, usecols=range(12), names=cols)
    nm = pd.read_csv(paf, sep="\t", header=None, usecols=[12], names=["nm"])["nm"].str[5:].astype(int)
    hits["nm"] = nm
    hits = hits[hits.qe - hits.qs >= 0.8 * hits.qlen].sort_values("mapq", ascending=False)
    hits = hits.drop_duplicates("q")  # --secondary=no still allows supplementary pieces
    hits["nk"] = hits.q.str.split("_n").str[1].astype(int)

    # Loci: unitigs within `gap` bp on the same contig.
    loci = []
    for t, h in hits.sort_values(["t", "ts"]).groupby("t", sort=False):
        start, end, nk, nu, mq0 = None, None, 0, 0, 0
        for r in h.itertuples():
            if start is not None and r.ts - end > gap:
                loci.append((t, start, end, nu, nk, mq0))
                start = None
            if start is None:
                start, end, nk, nu, mq0 = r.ts, r.te, 0, 0, 0
            end = max(end, r.te)
            nk += r.nk
            nu += 1
            mq0 += r.mapq == 0
        loci.append((t, start, end, nu, nk, mq0))
    loci = pd.DataFrame(loci, columns=["contig", "start", "end", "unitigs", "kmers", "mapq0_unitigs"])
    tlen = hits.drop_duplicates("t").set_index("t").tlen
    loci["end_dist"] = np.minimum(loci.start, loci.contig.map(tlen) - loci.end)
    loci.to_csv(f"{stem}.loci.tsv", sep="\t", index=False)

    nk_all = len(km)
    placed_k = hits.nk.sum()
    nuc = loci[loci.contig != MITO]
    ranked = np.sort(loci.kmers.to_numpy())[::-1]
    top10 = ranked[:10].sum() / max(1, placed_k)
    n50 = int(ranked[np.searchsorted(np.cumsum(ranked), ranked.sum() / 2)]) if len(ranked) else 0
    return dict(
        target=target, set=mset, kmers=nk_all, unitigs=len(lens), unitig_med=int(np.median(lens)),
        unitig_max=int(lens.max()), placed_kmer_frac=round(placed_k / nk_all, 3),
        mismatch_med=float(hits.nm.median()), mapq0_kmer_frac=round(hits[hits.mapq == 0].nk.sum() / max(1, placed_k), 3),
        loci=len(loci), loci_ge10k=int((loci.kmers >= 10).sum()), top10_loci_frac=round(top10, 3),
        locus_kmers_n50=n50,
        mito_kmer_frac=round(loci[loci.contig == MITO].kmers.sum() / max(1, placed_k), 3),
        end_zone_kmer_frac=round(nuc[nuc.end_dist < END_ZONE].kmers.sum() / max(1, nuc.kmers.sum()), 3),
        secs=round(time.time() - t0))


ap = argparse.ArgumentParser()
ap.add_argument("targets", nargs="+")
ap.add_argument("--gap", type=int, default=500)
ap.add_argument("--set", default="nm2", choices=["nm2", "ed2"])
args = ap.parse_args()

path = "place_markers.tsv"
for t in args.targets:
    row = {k: (float(v) if isinstance(v, np.floating) else v) for k, v in place(t, args.gap, args.set).items()}
    print(row, flush=True)
    res = pd.DataFrame([row])
    if os.path.exists(path):
        old = pd.read_csv(path, sep="\t")
        if "set" not in old:
            old.insert(1, "set", "nm2")
        res = pd.concat([old[(old.target != t) | (old.set != args.set)], res])
    res.to_csv(path, sep="\t", index=False)
