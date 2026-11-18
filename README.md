# cladeprint

Species- and lineage-specific k-mer markers from public population genomes, for every group of a labelled
hierarchy at once. cladeprint down-weights clonal cohorts, requires each marker in every child group of its
group, removes k-mers that any out-group genome holds within two substitutions or two edits, and estimates
marker sensitivity by leave-one-out and leave-one-cluster-out without rebuilding. A marker bundle classifies
new samples (reads or assemblies) with a call per level, a confidence and flags for mixtures, excess k-mers,
unexplained sequence and low coverage.

Paper: *cladeprint: species and lineage markers from public population genomes of the Cryptococcus species
complex* (S. Chanarat; under review).
Archived release, marker bundles and per-run results: https://doi.org/10.5281/zenodo.23985538 (all versions).

## Build

Rust (edition 2024) and [KMC 3](https://github.com/refresh-bio/KMC) on the `PATH`:

```
cargo build --release
./target/release/cladeprint --help
```

The analysis scripts use Python 3.12 and the tools listed in `envs/genotype.yml`
(`conda env create -f envs/genotype.yml`).

## Use

Classify a sample against a released bundle (the bundles are in the Zenodo record):

```
cladeprint classify bundles/census --reads SAMPLE_1.fastq.gz SAMPLE_2.fastq.gz --name SAMPLE -o calls.tsv
cladeprint classify bundles/census --assembly genome.fasta --name GENOME
cladeprint classify bundles/census --samples sheet.tsv -o calls.tsv
```

Build a bundle from a configuration (counting, QC, markers for every group at every level, sensitivity
estimates, scores and the bundle):

```
cladeprint run config.toml
```

The configuration names a sample sheet, an output directory and the label levels. Every other key is optional;
the values below are the defaults, except `max_out`, which was 1 for the paper (default 0):

```toml
samples = "samples.tsv"
output  = "out"
levels  = ["species", "lineage", "sublineage"]

[kmers]
k = 31
min_peak_coverage = 4          # below: LOW_COVERAGE, not used for training

[markers]
f_in = 0.9                     # share of a group's (cluster-weighted) training runs holding a marker
min_group_size = 3             # training runs a group needs to get markers
max_out = 1                    # out-group training runs allowed to hold a marker (default 0)
balance_children = true        # hierarchy-aware prevalence: f_in also within every child group
near_miss = true               # Hamming and edit-distance (<= 2) near-miss filters
max_candidates = 1000000       # hash-sample cap per group (0 = no cap)

[training]
min_kmer_excess = 0.80         # below: partial or foreign genome; does not train
max_kmer_excess = 1.10         # above: two genomes (mixture or hybrid); does not train

[resources]
threads = 0                    # 0 = all cores
memory_gb = 4                  # KMC memory
score_batch = 256              # samples scored at once
# tmp = "out/tmp"              # default <output>/tmp, on disk
```

The sample sheet is a TSV with a `sample` column, either `reads` (one or two comma-separated FASTQ files) or
`assembly` (FASTA), one column per label level named as in `levels` (empty = unknown at that level), and
optional `train` (`auto`, `yes` or `no`) and `cluster` (genotype cluster, from `cladeprint clusters`) columns.
Relative paths are resolved from the configuration file's directory.

The individual steps (`count`, `subsample`, `clusters`, `groups`, `near-miss`, `edit-miss`, `unitigs`,
`score`) are also available as subcommands; `cladeprint <subcommand> --help` lists their options.

## Repository layout

| Path | Contents |
|---|---|
| `crates/` | Rust library (`cladeprint-core`) and command-line tool (`cladeprint-cli`) |
| `census/` | Scripts for the public-data census, species screen, holdouts, assay audit and primer design |
| `labels/` | Literature labels, corrections and the published assay oligos; sources in `labels/*/SOURCES.md` |
| `vnii/` | Analysis of the VNII gene cluster |
| `pilot/`, `kmers/`, `external/` | Pilot analyses and tool comparisons |
| `figures/` | Figure scripts; `figures/msystems/README.md` maps every figure of the paper to its script and inputs |

Large intermediate data (per-run k-mer databases, about 0.5 TB) are not included; they can be regenerated
from the public SRA runs listed in Supplementary Table S3.

## Licence

Dual-licensed under the MIT licence or the Apache License 2.0, at your option (`LICENSE-MIT`,
`LICENSE-APACHE`).
