#!/usr/bin/env bash
# Species screen: map a read subsample competitively against the 7 species of the
# C. neoformans/gattii complex; report the share of reads each species wins (MAPQ>=20,
# nuclear contigs >=100 kb only) and the share that maps nowhere.
# Usage: micromamba run -n cladeprint bash species_screen.sh RUN [RUN...]   (OUT=file, default screen/species.tsv: rewritten)
set -euo pipefail
cd "$(dirname "$0")"
PANEL=ref/panel/complex.mmi
NPAIRS=${NPAIRS:-200000}
OUT=${OUT:-screen/species.tsv}
mkdir -p screen
awk '$2>=100000{print $1}' ref/panel/complex.fna.fai > screen/nuclear.txt 2>/dev/null \
  || { samtools faidx ref/panel/complex.fna; awk '$2>=100000{print $1}' ref/panel/complex.fna.fai > screen/nuclear.txt; }

echo -e "run\treads\tunmapped\tlowmapq\tneoformans\tdeneoformans\tgattii\tdeuterogattii\tbacillisporus\ttetragattii\tdecagattii" > $OUT
for run in "$@"; do
  urls=$(curl -sS --retry 3 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$run&result=read_run&fields=fastq_ftp&format=tsv" | tail -1 | cut -f2 | tr ';' ' ')
  r1=$(tr ' ' '\n' <<<"$urls" | grep -E '_1\.fastq' || true)
  r2=$(tr ' ' '\n' <<<"$urls" | grep -E '_2\.fastq' || true)
  [[ -n $r1 && -n $r2 ]] || { echo "FAIL $run (no ENA paired fastq)" >&2; continue; }
  n=$((NPAIRS * 4))
  (curl -sS --retry 3 "https://$r1" | zcat 2>/dev/null | head -n $n || true) > screen/${run}_1.fq
  (curl -sS --retry 3 "https://$r2" | zcat 2>/dev/null | head -n $n || true) > screen/${run}_2.fq
  minimap2 -ax sr -t 4 $PANEL screen/${run}_1.fq screen/${run}_2.fq 2>/dev/null \
    | samtools view -F 0x900 - \
    | awk -v run=$run 'BEGIN{OFS="\t"; while((getline l < "screen/nuclear.txt")>0) nuc[l]=1;
        split("neoformans deneoformans gattii deuterogattii bacillisporus tetragattii decagattii",sp," ")}
      {tot++; if(and($2,4)){un++; next} if($5<20 || !($3 in nuc)){lq++; next} split($3,a,"|"); c[a[1]]++}
      END{printf "%s\t%d\t%.3f\t%.3f", run, tot, un/tot, lq/tot; for(i=1;i<=7;i++) printf "\t%.3f", c[sp[i]]/tot; print ""}' \
    | tee -a $OUT
  rm -f screen/${run}_[12].fq
done
