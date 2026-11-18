#!/usr/bin/env bash
# Per-isolate k-mer presence databases (k=31, canonical) from the genotyping BAMs.
#
# The BAMs keep every read (mapped or not), so no new download is needed. Per isolate:
#   1. count all k-mers (KMC)
#   2. error cutoff = the valley between the error peak (count 1) and the genome peak,
#      read from that isolate's histogram (design step 2)
#   3. keep k-mers at or above the cutoff, with counts set to 1 -> presence database
# stats.tsv records the cutoff, the genome-peak coverage and the solid k-mer count, so
# low-coverage isolates (false absences) can be flagged downstream.
#
# Usage: micromamba run -n cladeprint bash build_dbs.sh RUN [RUN...]   (default: all BAMs)
set -euo pipefail
cd "$(dirname "$0")"
K=31
mkdir -p db tmp hist
[[ -f stats.tsv ]] || echo -e "run\tcutoff\tpeak_cov\tsolid_kmers\terror_kmers\tlow_cov" > stats.tsv
runs=("$@")
[[ ${#runs[@]} -gt 0 ]] || runs=($(ls ../pilot/bam/*.bam | xargs -n1 basename | sed 's/\.bam$//'))

for r in "${runs[@]}"; do
  [[ -s db/$r.kmc_suf ]] && continue
  samtools fastq -F 0x900 -@2 ../pilot/bam/$r.bam 2>/dev/null > tmp/$r.fq
  kmc -k$K -m2 -t6 -ci1 -cs1000 tmp/$r.fq tmp/$r.all tmp/ > /dev/null 2>&1
  rm -f tmp/$r.fq
  kmc_tools transform tmp/$r.all histogram hist/$r.txt -cx1000 > /dev/null 2>&1
  # cutoff: first local minimum after count 1 (at least 2); peak: highest bin above the cutoff.
  # No valley within 20 (histogram only falls: coverage too low for a genome peak) -> cutoff 2,
  # flagged low_cov. Also low_cov if the genome peak is below 5x.
  read cut peak solid err low < <(awk '{h[$1]=$2} END{
      c=2; while (c<20 && h[c+1]<h[c]) c++;
      low=0; if (c>=20) {c=2; low=1}
      p=c; for (i=c; i<=1000; i++) if (h[i]>h[p]) p=i;
      if (low) p=0; else if (p<5) low=1;
      s=0; for (i=c; i<=1000; i++) s+=h[i]; e=0; for (i=1; i<c; i++) e+=h[i];
      print c, p, s, e, low}' hist/$r.txt)
  kmc_tools transform tmp/$r.all -ci$cut set_counts 1 db/$r > /dev/null 2>&1
  rm -f tmp/$r.all.kmc_pre tmp/$r.all.kmc_suf
  echo -e "$r\t$cut\t$peak\t$solid\t$err\t$low" | tee -a stats.tsv
done
