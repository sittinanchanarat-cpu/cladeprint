#!/usr/bin/env bash
# Gene models on the VNII long-read genome (GCA_022832995.1) for align.sh / align_fast.sh:
# miniprot (spliced protein-to-genome) of the 11 WM276 block + flank proteins and of the H99 proteome;
# H99 models overlapping a block model, or a better H99 model, are dropped. Writes genes.saf (CDS
# exons, GeneID = protein; block genes prefixed BLOCK:) and the hisat2 index vnii_idx.
set -euo pipefail
cd "$(dirname "$0")"
G=../junction/vnii_ont.fna
micromamba run -n rnaseq hisat2-build -q -p 4 $G vnii_idx
micromamba run -n cladeprint miniprot -t 4 -d vnii.mpi $G
micromamba run -n cladeprint miniprot -t 4 --gff --outn=1 vnii.mpi ../genes/block_proteins.faa > block.gff
micromamba run -n cladeprint miniprot -t 4 --gff --outn=1 vnii.mpi ../genes/h99_protein.faa > h99.gff
python3 - <<'PY'
def models(f, tag):
    m, cds = {}, {}
    for l in open(f):
        x = l.rstrip("\n").split("\t")
        if l.startswith("#") or len(x) < 9:
            continue
        a = dict(kv.split("=", 1) for kv in x[8].split(";") if "=" in kv)
        if x[2] == "mRNA":
            m[a["ID"]] = dict(id=tag + a["Target"].split()[0], chr=x[0], s=int(x[3]), e=int(x[4]), st=x[6],
                              ident=float(a.get("Identity", 0)))
        elif x[2] == "CDS":
            cds.setdefault(a["Parent"], []).append((int(x[3]), int(x[4])))
    for k in m:
        m[k]["cds"] = cds.get(k, [])
    return list(m.values())

B, H = models("block.gff", "BLOCK:"), models("h99.gff", "")
keep = sorted((h for h in H if not any(h["chr"] == b["chr"] and h["s"] < b["e"] and b["s"] < h["e"] for b in B)),
              key=lambda h: -h["ident"])
final, occ = [], {}
for h in keep:
    if any(h["s"] < o[1] and o[0] < h["e"] for o in occ.get(h["chr"], [])):
        continue
    occ.setdefault(h["chr"], []).append((h["s"], h["e"]))
    final.append(h)
with open("genes.saf", "w") as f:
    f.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    for g in B + final:
        for s, e in g["cds"]:
            f.write(f"{g['id']}\t{g['chr']}\t{s}\t{e}\t{g['st']}\n")
print(len(B), "block models,", len(H), "H99 models ->", len(final), "non-overlapping kept")
PY
