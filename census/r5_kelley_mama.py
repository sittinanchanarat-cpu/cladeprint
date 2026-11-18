"""Kelley et al. 2014 MAMA assays scored on the census genomes, for Fig. 8C.

These values were previously copied into figures/fig8_assays.py as a literal, with a note that
they were recomputable only from census unitigs that are no longer stored. That was wrong: the
Kelley oligos were part of the R5 audit, so their per-genome binding is in census/out/r5_sites.tsv
and the rates recompute from it exactly.

Scored as published: a genome is called positive for a MAMA assay when the match primer binds with
fewer mismatches than the mismatch primer (Kelley et al. 2014; the universal primer is shared by
both reactions and does not discriminate, so it does not enter the score). Groups and the genome
universe follow census/design3p_screen_summary.py, so the diamonds and the designed-pair circles in
Fig. 8C sit on the same genomes: the training runs, grouped by second-round label.

Writes census/out/r5_kelley_mama.tsv.
Usage: micromamba run -n cladeprint python census/r5_kelley_mama.py
"""

import csv
from collections import defaultdict
from pathlib import Path

OUT = Path("census/out")
# MAMA assay -> the C. gattii-complex species it targets (VGI-VGIV)
TARGET = {"VGI-MPD471": ("VGI", "C.gattii"),
          "VGII-MPD495": ("VGII", "C.deuterogattii"),
          "VGIII-MPD198": ("VGIII", "C.bacillisporus"),
          "VGIV-MPD423": ("VGIV", "C.tetragattii")}


def mm(x):
    """Mismatch count, or None where primer_sites found no site within its search radius."""
    try:
        return int(x)
    except (TypeError, ValueError):
        return None


def main():
    qc = {r["sample"]: r for r in csv.DictReader(open(OUT / "qc.tsv"), delimiter="\t")}
    lab = {r["sample"]: r for r in csv.DictReader(open(OUT / "prop/samples.tsv"), delimiter="\t")}
    universe = {s for s, r in qc.items() if r["train"] == "true"}

    role = {}
    for base, (vg, _) in TARGET.items():
        role[f"{base} match MAMA"] = (vg, "match")
        role[f"{base} mismatch MAMA"] = (vg, "mismatch")

    counts = defaultdict(dict)
    for r in csv.DictReader(open(OUT / "r5_sites.tsv"), delimiter="\t"):
        key = role.get(r["oligo"].split("|", 1)[-1])
        if key:
            counts[(r["sample"], key[0])][key[1]] = mm(r["mismatches"])

    rows = []
    for base, (vg, species) in TARGET.items():
        positive = set()
        for s in universe:
            d = counts.get((s, vg), {})
            a, b = d.get("match"), d.get("mismatch")
            if a is not None and b is not None and a < b:
                positive.add(s)
        member = {s for s in universe if lab.get(s, {}).get("species") == species}
        nonmember = universe - member
        rows.append(dict(assay=base, molecular_type=vg, target_species=species,
                         n_member=len(member), n_nonmember=len(nonmember),
                         member=round(len(positive & member) / len(member), 4),
                         nonmember=round(len(positive & nonmember) / len(nonmember), 6)))

    rows.sort(key=lambda r: r["molecular_type"])
    path = OUT / "r5_kelley_mama.tsv"
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(universe)} training runs; wrote {path}")
    for r in rows:
        print(f"  {r['molecular_type']:6s} {r['target_species']:18s} n={r['n_member']:4d}"
              f"  member={r['member']:.4f}  non-member={r['nonmember']:.6f}")


if __name__ == "__main__":
    main()
