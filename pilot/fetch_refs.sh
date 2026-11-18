#!/usr/bin/env bash
# Download references: H99 (mapping) and the two species-screen panels:
#   complex   - the 7 species of the Cryptococcus complex; nuclear contigs = >= 100 kb
#   sporothrix - S. brasiliensis, schenckii, globosa, luriei (clinical clade), pallida, mexicana;
#                several are fragmented, so nuclear contigs = >= 5 kb and GC >= 35% (the
#                mitochondrial contigs are 25-36 kb at ~25% GC)
# Each panel: ref/panel/NAME.fna (contigs named species|accession), NAME.mmi, NAME.species.txt
# (species in column order of the screen output), NAME.nuclear.txt. Panels already built are kept.
# Usage: micromamba run -n cladeprint bash fetch_refs.sh
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p ref/panel
get() {  # accession -> FASTA on stdout
  curl -sSL --retry 5 "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$1/download?include_annotation_type=GENOME_FASTA" -o ref/$1.zip
  unzip -p ref/$1.zip '*.fna'; rm -f ref/$1.zip
}
panel() {  # name min_len min_gc species:accession...
  local name=$1 min_len=$2 min_gc=$3; shift 3
  local P=ref/panel/$name
  [[ -s $P.mmi && -s $P.nuclear.txt ]] && return 0
  rm -f $P.fna
  for x in "$@"; do
    [[ -s ref/panel/${x%%:*}.fna ]] || get ${x#*:} | sed "s/^>/>${x%%:*}|/" > ref/panel/${x%%:*}.fna
    cat ref/panel/${x%%:*}.fna >> $P.fna
    echo ${x%%:*}
  done > $P.species.txt
  samtools faidx $P.fna
  awk -v L=$min_len -v G=$min_gc '/^>/{if(id!="" && n>=L && 100*g/n>=G) print id; id=substr($1,2); n=0; g=0; next}
    {n+=length($0); g+=gsub(/[GCgc]/,"")} END{if(n>=L && 100*g/n>=G) print id}' $P.fna > $P.nuclear.txt
  minimap2 -x sr -d $P.mmi $P.fna 2>/dev/null
  echo "panel $name: $(wc -l < $P.species.txt) species, $(wc -l < $P.nuclear.txt) of $(wc -l < $P.fna.fai) contigs nuclear"
}
if [[ ! -s ref/H99.fna.mmi ]]; then
  get GCF_000149245.1 > ref/H99.fna   # C. neoformans H99
  samtools faidx ref/H99.fna
  minimap2 -d ref/H99.fna.mmi ref/H99.fna 2>/dev/null
fi
[[ -s ref/panel/neoformans.fna ]] || sed 's/^>/>neoformans|/' ref/H99.fna > ref/panel/neoformans.fna
panel complex 100000 0 neoformans:GCF_000149245.1 deneoformans:GCF_000091045.1 gattii:GCF_000185945.1 \
  deuterogattii:GCF_002954075.1 bacillisporus:GCF_000836335.1 tetragattii:GCF_000835755.1 decagattii:GCF_036417295.1
panel sporothrix 5000 35 brasiliensis:GCF_000820605.1 schenckii:GCF_000961545.1 globosa:GCA_001630435.1 \
  luriei:GCA_021398005.1 pallida:GCA_021396235.1 mexicana:GCA_021396375.1
echo "refs ready"
