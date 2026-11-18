"""Local assembly of the HOC1-RGD1 locus per isolate, and where each contig switches between
C. neoformans-type (closer to H99 chr7) and C. gattii-type (closer to WM276 chrF) sequence.

Per run: megahit on the read pairs touching either locus window (map_locus.sh), then each contig is
aligned to both references (minimap2 asm20, best hit only) and cut into 200 bp bins. A bin is N or G
when its differences to one reference are under half those to the other, '.' when only one reference
aligns there or the call is unclear. Prints contig segments with reference coordinates.
Usage: python junctions.py RUN...
"""
import re, subprocess, sys
from pathlib import Path

BIN = 200
REFS = {"N": "h99chr7.fa", "G": "wm276chrF.fa"}


def sh(*a):
    return subprocess.run(a, capture_output=True, text=True, check=True).stdout


def align(ref, fa):
    """Best alignment per contig: query span, ref span, strand, query positions of differences, query->ref map."""
    out = {}
    for line in sh("minimap2", "-c", "--cs", "-x", "asm20", ref, fa).splitlines():
        f = line.split("\t")
        if f[0] in out:
            continue  # minimap2 prints the primary first
        qs, qe, strand, ts = int(f[2]), int(f[3]), f[4], int(f[7])
        cs = next(t for t in f if t.startswith("cs:Z:"))[5:]
        diffs, anchors, q, t = [], [], 0, 0  # q, t: offsets along the aligned query (in ref strand) and ref
        for op, val in re.findall(r"([:*+\-])([0-9a-z]+)", cs):
            if op == ":":
                anchors.append((q, t)); q += int(val); t += int(val)
            elif op == "*":
                diffs.append(q); q += 1; t += 1
            elif op == "+":
                diffs.append(q); q += len(val)
            else:
                diffs.append(q); t += len(val)
        qlen = qe - qs
        # convert offsets to contig coordinates (minus strand: cs runs along the reverse complement)
        tocontig = (lambda o: qs + o) if strand == "+" else (lambda o: qe - 1 - o)
        out[f[0]] = dict(qs=qs, qe=qe, diffs=sorted(tocontig(d) for d in diffs),
                         ref=[(tocontig(a), ts + b) for a, b in anchors], strand=strand, id=int(f[9]) / int(f[10]))
    return out


def refpos(a, qpos):
    """Approximate ref coordinate of a contig position, from the nearest match anchor."""
    return min(a["ref"], key=lambda x: abs(x[0] - qpos))[1] + (qpos - min(a["ref"], key=lambda x: abs(x[0] - qpos))[0]) * (1 if a["strand"] == "+" else -1)


def calls(n_al, g_al, length):
    out = []
    for b in range(0, length - BIN + 1, BIN):
        have = {}
        for k, a in (("N", n_al), ("G", g_al)):
            if a and a["qs"] <= b and b + BIN <= a["qe"]:
                have[k] = sum(b <= d < b + BIN for d in a["diffs"])
        if len(have) == 2:
            c = "N" if have["N"] * 2 < have["G"] else "G" if have["G"] * 2 < have["N"] else "?"
        else:
            c = "".join(have).lower() or "-"  # aligned to one reference only
        out.append((b, c, have))
    return out


for run in sys.argv[1:]:
    asm = Path(f"asm_{run}")
    if not (asm / "final.contigs.fa").exists():
        subprocess.run(["megahit", "-1", f"{run}.locus.1.fq", "-2", f"{run}.locus.2.fq", "-o", str(asm), "-t", "2",
                        "--min-contig-len", "1000"], capture_output=True, check=True)
    fa = str(asm / "final.contigs.fa")
    lens = {l[1:].split()[0]: int(l.split("len=")[1]) for l in open(fa) if l.startswith(">")}
    al = {k: align(r, fa) for k, r in REFS.items()}
    print(f"\n### {run}: {len(lens)} contigs")
    for c in sorted(lens, key=lambda c: -lens[c]):
        n, g = al["N"].get(c), al["G"].get(c)
        if not (n or g):
            continue
        cl = calls(n, g, lens[c])
        track = "".join(x[1] for x in cl)
        print(f"{c} len {lens[c]}  track (1 char = {BIN} bp; N/G = closer to H99/WM276, n/g = only that one aligns):")
        print("  " + track)
        # segment boundaries where a definite call changes N <-> G
        last = None
        for b, ch, _ in cl:
            if ch in "NG":
                if last and ch != last[1]:
                    mid = (last[0] + BIN + b) // 2
                    print(f"  switch {last[1]}->{ch} at contig {last[0] + BIN}-{b}: "
                          f"H99 ~{refpos(n, mid):,}  WM276 ~{refpos(g, mid):,}")
                last = (b, ch)
        for k, a in (("H99", n), ("WM276", g)):
            if a:
                print(f"  {k}: contig {a['qs']}-{a['qe']} -> {min(x[1] for x in a['ref']):,}-{max(x[1] for x in a['ref']):,} ({a['strand']}) id {a['id']:.3f}")
