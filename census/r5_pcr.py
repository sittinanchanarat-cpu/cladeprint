"""R5 product-level audit: in-silico PCR of the published assays on public NCBI assemblies.

The census audit (r5_summary.py) asks whether each oligo binds somewhere in a genome's unitigs;
that cannot tell a product from two unrelated sites (Lengeler MFalpha2: in C. neoformans the two
primers bind 3-mismatch sites on different unitigs). Assemblies are long enough for products.

Input
  census/asm/r5_asm_sites.tsv     primer_sites --all on the assemblies (every site <= 3 mismatches)
  census/asm/classify_subclade.tsv  cladeprint classify of the assemblies (subclade bundle): labels
  census/asm/assemblies.tsv       NCBI metadata (organism names for species outside the bundle)
  labels/assays/assay_list_for_review_audited.tsv, r5_oligos.tsv, r5_panel_pairs.tsv

Reactions: every assay with its oligos in one tube (a multiplex gives its band pattern), except the
singleplex panels (Lengeler, Yan, Carvalho), which are taken pair by pair from r5_panel_pairs.tsv.
Skipped: duplicate listings, LAMP, MAMA (scored relatively in r5_summary.py), fingerprinting.

Binding (three rules, all reported; <= 3 mismatches in each): strict = none in the last 5 bases;
soft = the 3'-terminal base may mismatch, bases 2-5 from the 3' end may not; lenient = only bases
2-3 from the 3' end must match. Strict fails validated assays on their own reference strains
(Lengeler MFalpha2 in JEC21: terminal mismatch; Carvalho MATalpha-A in H99: mismatch 4 bases from
the 3' end, 1,213 bp product vs 1,200 published), so it is a lower bound, lenient an upper one.
Product: two binding sites facing each other on one contig, 30-3000 bp (any two oligos of the
reaction, either orientation, since F/R labels are not always stated). Expected: within 10% (at
least 20 bp) of a size the source gives; when it gives none, any product counts. Probe reactions
also need each probe inside the product within 1 mismatch (any position).

Oligos are searched as their binding sequence (binding_5to3 when given: Lengeler JOHE1895/1896
without their NdeI/KpnI cloning tails). The A048+A050v pair is a sensitivity variant (JOHE2189 with
the 3' end that serotype A genomes carry; the printed one mismatches all of them).

Scored against each reaction's intended target (labels/assays/r5_intended_targets.tsv: serotype A
= C. neoformans, D = C. deneoformans, B/C = the gattii complex), with the complex's other species
as non-targets; multiplex typing assays ("pattern") give their band pattern per species instead.

Usage: python census/r5_pcr.py
Writes census/asm/r5_pcr_products.tsv, r5_pcr_by_group.tsv, r5_pcr_intended.tsv,
labels/assays/r5_pcr_reactions.tsv.
"""

import csv
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

ASM = Path("census/asm")
LAB = Path("labels/assays")
MIN_BP, MAX_BP = 30, 3000
SKIP = (
    "QSP1 species-specific qPCR (Mbangiwa et al. 2024)",  # duplicate of the four Mbangiwa rows
    "URA5 PCR for URA5-RFLP molecular typing",  # same primers as URA5-RFLP (Meyer et al. 2003)
    "Meyer 2003 URA5-RFLP original description",  # not retrieved
    "Lucas 2010 CAP59 allele-specific LAMP",  # LAMP: no simple product
    "Lucas 2010 CAP59 allele-specific LAMP (duplicate of A103)",
    "PCR-fingerprinting (M13)",
    "PCR-fingerprinting (microsatellite)",
)
PANELS = ("Lengeler 2001", "Yan 2002", "Carvalho 2007")
GATTII = {"C.gattii", "C.deuterogattii", "C.bacillisporus", "C.tetragattii", "C.decagattii"}
COMPLEX = GATTII | {"C.neoformans", "C.deneoformans"}


def target_of(name, assay, targets):
    """Intended species of a reaction: a set, "pattern", or None (not judged)."""
    key = next((k for k in targets if assay.startswith(k)), None)
    if key is None:
        return None
    t = targets[key]
    if t == "by_pair":
        txt = PANEL_TARGET.get(name, "")
        a = re.search(r"serotypes? A\b|\bA and D|D and A", txt)
        d = re.search(r"serotypes? D\b|\bA and D|D and A", txt)
        return ({"C.neoformans"} if a else set()) | ({"C.deneoformans"} if d else set()) or None
    return {"complex": COMPLEX, "gattii_complex": GATTII}.get(t, {t} if t != "pattern" else "pattern")


PANEL_TARGET = {}


def sizes(text):
    """Product sizes (bp) stated in a free-text field: integers >= 30 not part of a name."""
    return sorted({int(x) for x in re.findall(r"(?<![\w.-])(\d{2,4})(?![\w.])", text or "") if int(x) >= 30})


def reactions(alist, seq_rep):
    out = []  # (reaction, assay, oligo ids, probe ids, expected sizes)
    by = defaultdict(list)
    for r in alist:
        by[r["assay"]].append(r)
    for assay, rs in by.items():
        if assay in SKIP or "MAMA" in assay or assay.startswith(PANELS):
            continue
        rs = [r for r in rs if r["sequence_5to3"] and r["id"] in seq_rep]
        prim = [r["id"] for r in rs if not r["role"].lower().startswith("probe")]
        probes = [r["id"] for r in rs if r["role"].lower().startswith("probe")]
        if len(prim) < 2:
            continue
        exp = sorted({s for r in rs for s in sizes(r["expected_product_bp"])})
        out.append((assay, assay, prim, probes, exp))
    meta = {r["id"]: r for r in alist}
    for p in csv.DictReader(open(LAB / "r5_panel_pairs.tsv"), delimiter="\t"):
        ids = [p["f_id"], p["r_id"]]
        exp = sorted({s for i in ids for s in sizes(meta[i]["expected_product_bp"])})
        out.append((f"{p['assay']}: {p['pair']}", p["assay"], ids, [], exp))
        # the reverse primer's target names the reaction (Lengeler reuses JOHE2926/2596 forwards
        # for the D and A reactions, labelled "serotypes D and A")
        rt = meta[p["r_id"]]["intended_target"]
        single = len(set(re.findall(r"serotypes? ([AD])\b|\b([AD]) and|and ([AD])\b", rt))) == 1 and " and " not in rt
        PANEL_TARGET[f"{p['assay']}: {p['pair']}"] = rt if single else p["target"]
    return out


def main():
    alist = list(csv.DictReader(open(LAB / "assay_list_for_review_audited.tsv"), delimiter="\t"))
    # binding_5to3 = the part that anneals (cloning tails trimmed), when it differs from the print
    for r in alist:
        r["sequence_5to3"] = r["binding_5to3"] or r["sequence_5to3"]
    # the deduplicated oligo name primer_sites used for each assay oligo (shared sequences)
    rep_of_seq = {r["sequence"].upper(): r["oligo"] for r in csv.DictReader(open(LAB / "r5_oligos.tsv"), delimiter="\t")}
    seq_rep = {r["id"]: rep_of_seq[r["sequence_5to3"].upper()] for r in alist if r["sequence_5to3"].upper() in rep_of_seq}
    length = {r["id"]: len(r["sequence_5to3"]) for r in alist}
    rx = reactions(alist, seq_rep)
    with open(LAB / "r5_pcr_reactions.tsv", "w") as f:
        f.write("reaction\tassay\tprimers\tprobes\texpected_bp\n")
        for name, assay, prim, probes, exp in rx:
            f.write(f"{name}\t{assay}\t{','.join(prim)}\t{','.join(probes)}\t{','.join(map(str, exp))}\n")

    # sites per genome and representative oligo
    sites = defaultdict(list)  # (genome, rep) -> [(contig, pos, strand, mm, tail)]
    with open(ASM / "r5_asm_sites.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sites[(r["sample"], r["oligo"])].append((r["contig"], int(r["pos"]), r["strand"], int(r["mismatches"]), r["tail"]))

    def binds(tail, rule):
        if rule == "strict":
            return "x" not in tail
        if rule == "soft":
            return "x" not in tail[:-1]  # only the terminal base may mismatch
        return "x" not in tail[-3:-1]  # lenient: bases 2-3 from the 3' end must match

    genomes = [l.split("\t")[0] for l in open(ASM / "download.list").read().split()]
    rows = []
    for g in genomes:
        for name, assay, prim, probes, exp in rx:
            for rule in ("strict", "soft", "lenient"):
                plus, minus = defaultdict(list), defaultdict(list)
                for o in prim:
                    for c, pos, st, mm, tail in sites.get((g, seq_rep[o]), []):
                        if binds(tail, rule):
                            (plus if st == "+" else minus)[c].append((pos, o))
                prods = []
                for c, ps in plus.items():
                    for x, a in ps:
                        for y, b in minus.get(c, []):
                            size = y + length[b] - x
                            if MIN_BP <= size <= MAX_BP:
                                prods.append((c, x, y + length[b], a, b, size))
                for c, s0, e0, a, b, size in prods:
                    hit = [e for e in exp if abs(size - e) <= max(20, 0.1 * e)]
                    probe_ok = all(
                        any(pc == c and s0 <= pp and pp + length[p] <= e0 and mm <= 1
                            for pc, pp, _, mm, _ in sites.get((g, seq_rep[p]), []))
                        for p in probes
                    )
                    rows.append(dict(genome=g, reaction=name, rule=rule, a=a, b=b, size=size,
                                     expected=bool(hit) if exp else True, probe_ok=probe_ok,
                                     contig=c, start=s0))
    prod = pd.DataFrame(rows)
    prod.to_csv(ASM / "r5_pcr_products.tsv", sep="\t", index=False)

    # labels: cladeprint calls inside the complex, NCBI names outside it
    cl = pd.read_csv(ASM / "classify_subclade.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("sample")
    meta = pd.read_csv(ASM / "assemblies.tsv", sep="\t", dtype=str).set_index("accession")
    groups = defaultdict(set)
    for g in genomes:
        sp = cl.at[g, "species"] if g in cl.index else "NA"
        if sp == "NA":
            org = meta.at[g, "organism"] if g in meta.index else "?"
            groups["outside:" + " ".join(org.split()[:2])].add(g)
            continue
        for lv in ("species", "lineage", "sublineage", "subclade"):
            v = cl.at[g, lv]
            if v != "NA":
                groups[v].add(g)
    allg = set(genomes)
    out = []
    ok = prod[prod.expected & prod.probe_ok]
    for name, *_ in rx:
        for rule in ("strict", "soft", "lenient"):
            have = set(ok[(ok.reaction == name) & (ok.rule == rule)].genome)
            anyp = set(prod[(prod.reaction == name) & (prod.rule == rule)].genome)
            for grp, mem in sorted(groups.items()):
                non = allg - mem
                out.append(dict(reaction=name, rule=rule, group=grp, n=len(mem),
                                member_expected=round(len(have & mem) / len(mem), 3),
                                nonmember_expected=round(len(have & non) / max(1, len(non)), 3),
                                member_any_product=round(len(anyp & mem) / len(mem), 3)))
    pd.DataFrame(out).to_csv(ASM / "r5_pcr_by_group.tsv", sep="\t", index=False)

    # against the intended target
    targets = {r["key"]: r["target"] for r in csv.DictReader(open(LAB / "r5_intended_targets.tsv"), delimiter="\t")}
    species = {g: (cl.at[g, "species"] if g in cl.index else "NA") for g in genomes}
    # AD hybrids carry both parents' alleles: reported apart, not as target or non-target
    hybrid = {g for g in genomes if g in meta.index and " x " in meta.at[g, "organism"]}
    inside = {g for g, s in species.items() if s in COMPLEX} - hybrid
    outside = set(genomes) - inside - hybrid
    res = []
    for name, assay, *_ in rx:
        tg = target_of(name, assay, targets)
        for rule in ("strict", "soft", "lenient"):
            sel = ok[(ok.reaction == name) & (ok.rule == rule)]
            have = set(sel.genome)
            row = dict(reaction=name, rule=rule)
            if tg == "pattern":
                # band pattern per species: product sizes to the nearest 10 bp, most common pattern
                pat = sel.groupby("genome")["size"].apply(lambda s: "+".join(map(str, sorted({int(round(x, -1)) for x in s}))))
                for sp in sorted(COMPLEX):
                    mem = [g for g in inside if species[g] == sp]
                    if not mem:
                        continue
                    vc = pat.reindex(mem).fillna("none").value_counts()
                    top = sorted(vc.index[vc == vc.iloc[0]])  # ties: all of them, in a fixed order
                    row[sp] = f"{' | '.join(top)} ({vc.iloc[0]}/{len(mem)}{' each' if len(top) > 1 else ''})"
                row["target"] = "pattern"
            elif tg:
                mem = [g for g in inside if species[g] in tg]
                non = [g for g in inside if species[g] not in tg]
                row.update(target="+".join(sorted(tg)) if tg != COMPLEX else "complex",
                           n_target=len(mem), target_hit=round(len(have & set(mem)) / max(1, len(mem)), 3),
                           n_nontarget=len(non), nontarget_hit=round(len(have & set(non)) / max(1, len(non)), 3) if non else None)
            else:
                row["target"] = "not judged"
            row["ad_hybrid_hit"] = f"{len(have & hybrid)}/{len(hybrid)}"
            row["outside_complex_hit"] = f"{len(have & outside)}/{len(outside)}"
            res.append(row)
    pd.DataFrame(res).to_csv(ASM / "r5_pcr_intended.tsv", sep="\t", index=False)
    print(f"{len(rx)} reactions x {len(genomes)} assemblies; {len(prod)} products; groups:",
          {k: len(v) for k, v in sorted(groups.items())})


if __name__ == "__main__":
    main()
