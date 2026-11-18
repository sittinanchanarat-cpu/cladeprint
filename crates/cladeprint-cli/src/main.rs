//! cladeprint command line. Reads and writes the sorted uint64 `.npy` k-mer files of the Python
//! prototype (`kmers/groups/`), so each step can replace its Python counterpart.

use anyhow::{Context, Result, bail};
use cladeprint_core::{
    classify, clusters, config, count, edit, groups, hamming, kmc, npy, pipeline, score, subsample,
    unitig,
};
use clap::{Parser, Subcommand};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Parser)]
#[command(
    version,
    about = "Lineage- and species-specific k-mer markers for fungal pathogens"
)]
struct Cli {
    /// Worker threads (default: all cores)
    #[arg(long, global = true)]
    threads: Option<usize>,
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// The whole pipeline from a config: counting, QC, markers for every group at every level,
    /// leave-one-out sensitivity, scores, and the marker bundle (see README.md)
    Run {
        /// Config (TOML)
        config: PathBuf,
    },
    /// Classify samples against a marker bundle: a call per level with confidence and flags
    Classify {
        /// Bundle directory (output of run: <output>/bundle)
        bundle: PathBuf,
        /// One sample: FASTQ file(s)
        #[arg(long, num_args = 1..=2, conflicts_with_all = ["assembly", "samples"])]
        reads: Vec<PathBuf>,
        /// One sample: assembly FASTA
        #[arg(long, conflicts_with = "samples")]
        assembly: Option<PathBuf>,
        /// Many samples: a sample sheet (sample, reads or assembly; other columns ignored)
        #[arg(long)]
        samples: Option<PathBuf>,
        /// Sample name for --reads/--assembly
        #[arg(long, default_value = "sample")]
        name: String,
        /// With --samples: use existing presence databases <dir>/<sample> (and <sample>.stats.toml
        /// when present) instead of counting the reads
        #[arg(long, requires = "samples")]
        kmc_db: Option<PathBuf>,
        /// Output TSV (default: stdout)
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// One column per group score as well
        #[arg(long)]
        long: bool,
        /// Temporary directory (on disk)
        #[arg(long, default_value = "tmp")]
        tmp: PathBuf,
        #[arg(long, default_value_t = 4)]
        memory_gb: u32,
    },
    /// Drop candidates that an out-group k-mer matches with <= 2 mismatches (near_miss.py)
    NearMiss(Filter),
    /// Drop candidates that an out-group k-mer matches within 2 edits, indels allowed (edit_miss.py)
    EditMiss(Filter),
    /// Chain k-mers into unitigs (FASTA, headers u<i>_n<kmers>)
    Unitigs {
        /// Canonical k-mers (.npy)
        kmers: PathBuf,
        /// Output FASTA
        #[arg(short, long)]
        output: PathBuf,
    },
    /// Genotype clusters for clade weighting: split-k-mer SNP distances between presence
    /// databases, and single-linkage clusters at each cutoff (the `cluster` column of `run`)
    Clusters {
        /// Sample sheet (TSV with a `sample` column; other columns ignored)
        #[arg(long)]
        samples: PathBuf,
        /// Directory with each sample's KMC presence database <sample>.kmc_pre/.kmc_suf (KMC 1 layout)
        #[arg(long)]
        kmc_db: PathBuf,
        /// Cluster cutoff in SNPs per Mb of shared genome, repeatable (one cluster column each)
        #[arg(long = "cutoff", required = true)]
        cutoffs: Vec<f64>,
        /// Exact distances for pairs that may be within this many SNPs per Mb (default: 4 x the
        /// largest cutoff, at least 50)
        #[arg(long)]
        exact_max: Option<f64>,
        /// Pairs sharing less covered genome than this (bp, estimated) get no distance
        #[arg(long, default_value_t = 1e6)]
        min_positions: f64,
        /// Hash-sampling scale for the screened sites and each sample's covered genome
        #[arg(long, default_value_t = 64)]
        scale: u64,
        /// Coarser hash-sampling scale for each pair's shared genome (a multiple of --scale)
        #[arg(long, default_value_t = 1024)]
        overlap_scale: u64,
        /// Ignore alleles held by fewer samples (2 drops run-private alleles: mostly sequencing
        /// errors of old error-rich runs, and private SNPs)
        #[arg(long, default_value_t = 1)]
        min_allele: u32,
        /// Output directory: clusters.tsv, pairs.tsv, screen_hist.tsv
        #[arg(short, long)]
        output: PathBuf,
    },
    /// Group-specific k-mers from per-isolate presence databases, one streaming pass (group_kmers.py)
    Groups {
        /// Sample sheet (TSV with a `run` column)
        #[arg(long)]
        samples: PathBuf,
        /// Directory with each run's KMC presence database <run>.kmc_pre/.kmc_suf (KMC 1 layout)
        #[arg(long)]
        kmc_db: PathBuf,
        /// Use only rows where COLUMN=VALUE (training isolates)
        #[arg(long, default_value = "use=train")]
        filter: String,
        /// Target groups as NAME=COLUMN:VALUE, repeatable
        #[arg(long = "target", default_values = ["VNII=major:VNII", "VNBI=major:VNBI", "VNBII=major:VNBII",
            "VNI=major:VNI", "VNIa=called:VNIa", "VNIb=called:VNIb", "VNIc=called:VNIc"])]
        targets: Vec<String>,
        /// Share of target isolates a marker must be in
        #[arg(long, default_value_t = 0.9)]
        f_in: f64,
        /// Output directory: <T>.strict.u64.npy, <T>.x.strict.u64.npy, [<T>.out.u64.npy], groups.tsv
        #[arg(short, long)]
        output: PathBuf,
        /// Also write each target's out-group k-mer set (input of near-miss/edit-miss; ~0.45 GB each)
        #[arg(long)]
        write_out: bool,
        /// Pieces k-mer space is processed in (more = less memory; 0 = ~32M records per piece)
        #[arg(long, default_value_t = 0)]
        batches: u64,
    },
    /// Presence database of one sample: reads via KMC (error cutoff from the histogram) or an assembly
    Count {
        /// FASTQ file(s), gzipped or not
        #[arg(long, num_args = 1..=2, conflicts_with = "assembly", required_unless_present = "assembly")]
        reads: Vec<PathBuf>,
        /// Assembly FASTA
        #[arg(long)]
        assembly: Option<PathBuf>,
        /// Database path without .kmc_pre/.kmc_suf
        #[arg(short, long)]
        output: PathBuf,
        /// Temporary directory (on disk)
        #[arg(long, default_value = "tmp")]
        tmp: PathBuf,
        #[arg(long, default_value_t = 4)]
        memory_gb: u32,
        #[arg(long, default_value_t = 4)]
        min_peak: u32,
    },
    /// Random spot subsample of a FASTQ stream (stdin -> stdout), deterministic per spot name, so
    /// chunks of a run processed separately give the same sample as the whole run
    Subsample {
        /// Share of spots to keep
        #[arg(long)]
        fraction: f64,
        #[arg(long, default_value_t = 20261012)]
        seed: u64,
        /// Write "spots kept_spots records" here
        #[arg(long)]
        counts: Option<PathBuf>,
    },
    /// Write a sorted .npy k-mer set as a KMC presence database (KMC 1 layout, as kmc_tools writes)
    KmcWrite {
        /// Sorted, unique k-mers (.npy)
        kmers: PathBuf,
        /// Database path without .kmc_pre/.kmc_suf
        #[arg(short, long)]
        output: PathBuf,
        #[arg(long, default_value_t = 31)]
        k: u32,
        /// LUT prefix length (k minus it must be a multiple of 4)
        #[arg(long, default_value_t = 7)]
        prefix_len: u32,
    },
    /// Write the k-mers of a KMC database (KMC 1 or 2/3 format) as a sorted .npy
    KmcDump {
        /// Database path without .kmc_pre/.kmc_suf
        db: PathBuf,
        /// Output .npy
        #[arg(short, long)]
        output: PathBuf,
    },
    /// Lineage scores and leave-one-out sensitivity per isolate (score_isolates.py)
    Score {
        /// Directory with <T>.ed2.u64.npy (markers) and <T>.x.ed2.u64.npy (need-1 extras)
        #[arg(long)]
        groups: PathBuf,
        /// Sample sheet with columns run, called, major, use, peak_cov
        #[arg(long)]
        samples: PathBuf,
        /// Directory with <run>.npy: the marker k-mers present in each run
        #[arg(long, required_unless_present = "kmc_db")]
        present: Option<PathBuf>,
        /// Or: directory with each run's KMC presence database <run>.kmc_pre/.kmc_suf
        #[arg(long, conflicts_with = "present")]
        kmc_db: Option<PathBuf>,
        /// Targets, in output column order
        #[arg(
            long,
            value_delimiter = ',',
            default_value = "VNII,VNBI,VNBII,VNI,VNIa,VNIb,VNIc"
        )]
        targets: Vec<String>,
        /// Output prefix: <prefix>scores.tsv and <prefix>loo_sensitivity.tsv
        #[arg(long, default_value = "")]
        prefix: String,
    },
}

#[derive(clap::Args)]
struct Filter {
    /// Candidate k-mers, canonical (.npy)
    #[arg(long)]
    candidates: PathBuf,
    /// Out-group k-mers, canonical (.npy)
    #[arg(long)]
    outgroup: PathBuf,
    /// Surviving candidates (.npy), in input order
    #[arg(short, long)]
    output: PathBuf,
}

/// Raises the soft open-file limit to the hard limit: group counting and clusters keep two files
/// open per database, and the common default soft limit (1,024) is reached at ~500 samples.
#[cfg(unix)]
fn raise_open_file_limit() {
    // SAFETY: getrlimit/setrlimit with a valid pointer to a local struct
    unsafe {
        let mut r = libc::rlimit {
            rlim_cur: 0,
            rlim_max: 0,
        };
        if libc::getrlimit(libc::RLIMIT_NOFILE, &mut r) == 0 && r.rlim_cur < r.rlim_max {
            r.rlim_cur = r.rlim_max;
            libc::setrlimit(libc::RLIMIT_NOFILE, &r);
        }
    }
}

#[cfg(not(unix))]
fn raise_open_file_limit() {}

fn main() -> Result<()> {
    raise_open_file_limit();
    let cli = Cli::parse();
    if let Some(n) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()?;
    }
    match cli.cmd {
        Cmd::Run { config } => {
            let cfg = config::Config::load(&config)?;
            if cfg.resources.threads > 0 && cli.threads.is_none() {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(cfg.resources.threads)
                    .build_global()?;
            }
            let b = pipeline::run(&cfg)?;
            eprintln!(
                "bundle: {} groups ({} skipped) in {}",
                b.groups.len(),
                b.skipped.len(),
                cfg.output.join("bundle").display()
            );
            Ok(())
        }
        Cmd::Classify {
            bundle,
            reads,
            assembly,
            samples,
            name,
            kmc_db,
            output,
            long,
            tmp,
            memory_gb,
        } => {
            let samples: Vec<(String, config::Input)> = match (samples, assembly) {
                (Some(sheet), _) => config::load_samples(&sheet, &[])?
                    .into_iter()
                    .map(|s| (s.name, s.input))
                    .collect(),
                (None, Some(a)) => vec![(name, config::Input::Assembly(a))],
                (None, None) if !reads.is_empty() => vec![(name, config::Input::Reads(reads))],
                _ => bail!("give --reads, --assembly or --samples"),
            };
            run_classify(
                &bundle,
                &samples,
                &ClassifyOpts {
                    output: output.as_deref(),
                    long,
                    tmp: &tmp,
                    memory_gb,
                    threads: cli.threads.unwrap_or(0),
                    kmc_db: kmc_db.as_deref(),
                },
            )
        }
        Cmd::NearMiss(f) => filter(&f, "hamming", hamming::candidate_dist, hamming::FAR),
        Cmd::EditMiss(f) => filter(&f, "edit", edit::candidate_dist, edit::FAR),
        Cmd::Unitigs { kmers, output } => {
            let km = npy::read_u64(&kmers)?;
            let u = unitig::unitigs(&km);
            let mut w = BufWriter::new(File::create(&output)?);
            for (i, x) in u.iter().enumerate() {
                writeln!(w, ">u{i}_n{}\n{}", x.n_kmers, x.seq)?;
            }
            w.flush()?;
            eprintln!("{} k-mers -> {} unitigs", km.len(), u.len());
            Ok(())
        }
        Cmd::Clusters {
            samples,
            kmc_db,
            cutoffs,
            exact_max,
            min_positions,
            scale,
            overlap_scale,
            min_allele,
            output,
        } => run_clusters(
            &samples,
            &kmc_db,
            &cutoffs,
            exact_max,
            min_positions,
            scale,
            overlap_scale,
            min_allele,
            &output,
        ),
        Cmd::Groups {
            samples,
            kmc_db,
            filter,
            targets,
            f_in,
            output,
            write_out,
            batches,
        } => run_groups(
            &samples, &kmc_db, &filter, &targets, f_in, &output, write_out, batches,
        ),
        Cmd::Count {
            reads,
            assembly,
            output,
            tmp,
            memory_gb,
            min_peak,
        } => {
            let s = match assembly {
                Some(a) => count::count_assembly(&a, &output)?,
                None => count::count_reads(
                    &reads,
                    &output,
                    &tmp,
                    cli.threads.unwrap_or(0),
                    memory_gb,
                    min_peak,
                    0, // `run` reuses these databases: same rule as training
                )?,
            };
            // <db>.stats.toml, written last: marks the database complete, and `run` reuses it
            let mut sf = output.clone().into_os_string();
            sf.push(".stats.toml");
            fs::write(&sf, toml::to_string(&s)?)?;
            println!(
                "cutoff\t{}\npeak\t{}\nsolid_kmers\t{}\nerror_kmers\t{}\nlow_coverage\t{}",
                s.cutoff,
                s.peak.map_or("-".into(), |p| p.to_string()),
                s.solid_kmers,
                s.error_kmers,
                s.low_coverage
            );
            Ok(())
        }
        Cmd::Subsample {
            fraction,
            seed,
            counts,
        } => {
            let stdin = std::io::stdin().lock();
            let input = std::io::BufReader::with_capacity(1 << 20, stdin);
            let out = BufWriter::with_capacity(1 << 20, std::io::stdout().lock());
            let c = subsample::subsample(input, out, fraction, seed)?;
            if let Some(p) = counts {
                fs::write(p, format!("{} {} {}\n", c.spots, c.kept_spots, c.records))?;
            }
            Ok(())
        }
        Cmd::KmcWrite {
            kmers,
            output,
            k,
            prefix_len,
        } => kmc::write_kmc1(&output, &npy::read_u64(&kmers)?, k, prefix_len),
        Cmd::KmcDump { db, output } => {
            let v = kmc::read_sorted(&db)?;
            npy::write_u64(&output, &v)?;
            eprintln!("{} k-mers", v.len());
            Ok(())
        }
        Cmd::Score {
            groups,
            samples,
            present,
            kmc_db,
            targets,
            prefix,
        } => {
            let src = match (present, kmc_db) {
                (Some(p), _) => Present::Npy(p),
                (None, Some(d)) => Present::Kmc(d),
                (None, None) => bail!("give --present or --kmc-db"),
            };
            run_score(&groups, &samples, &src, &targets, &prefix)
        }
    }
}

fn filter(f: &Filter, name: &str, dist: fn(&[u64], &[u64]) -> Vec<u8>, far: u8) -> Result<()> {
    let t0 = Instant::now();
    let cand = npy::read_u64(&f.candidates)?;
    let out = npy::read_u64(&f.outgroup)?;
    let d = dist(&cand, &out);
    let keep: Vec<u64> = cand
        .iter()
        .zip(&d)
        .filter(|&(_, &x)| x == far)
        .map(|(&c, _)| c)
        .collect();
    npy::write_u64(&f.output, &keep)?;
    let n = |v: u8| d.iter().filter(|&&x| x == v).count();
    println!(
        "{name}\tcandidates={}\tout_kmers={}\td0={}\td1={}\td2={}\tsurvive={}\tsecs={:.1}",
        cand.len(),
        out.len(),
        n(0),
        n(1),
        n(2),
        keep.len(),
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

struct ClassifyOpts<'a> {
    output: Option<&'a Path>,
    long: bool,
    tmp: &'a Path,
    memory_gb: u32,
    threads: usize,
    kmc_db: Option<&'a Path>,
}

fn run_classify(
    bundle: &Path,
    samples: &[(String, config::Input)],
    o: &ClassifyOpts,
) -> Result<()> {
    let (output, long, tmp, memory_gb, threads) = (o.output, o.long, o.tmp, o.memory_gb, o.threads);
    let b = classify::load(bundle)?;
    let levels = &b.bundle.levels;
    let mut out: Box<dyn Write> = match output {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(std::io::stdout().lock()),
    };
    let mut head = vec![
        "sample".to_string(),
        "peak".into(),
        "low_coverage".into(),
        "unexplained".into(),
        "kmer_excess".into(),
    ];
    for l in levels {
        head.extend([l.clone(), format!("{l}_score"), format!("{l}_confidence")]);
    }
    head.push("flags".into());
    if long {
        head.extend(b.bundle.groups.iter().map(|g| format!("score:{}", g.id)));
    }
    writeln!(out, "{}", head.join("\t"))?;
    fs::create_dir_all(tmp)?;
    // a valley with fewer solid k-mers above it than this share of the smallest reference genome
    // did not separate a genome from errors (see count::cutoff_from_histogram)
    let min_solid = b
        .bundle
        .solid_reference
        .values()
        .fold(f64::INFINITY, |m, &r| m.min(r));
    let min_solid = if min_solid.is_finite() {
        (classify::MIN_GENOME_SHARE * min_solid) as u64
    } else {
        0
    };
    for (name, input) in samples {
        let ((present, unexplained), st) = if let Some(dir) = o.kmc_db {
            // existing presence database; stats when a .stats.toml sits next to it
            let stats = match fs::read_to_string(dir.join(format!("{name}.stats.toml")))
                .ok()
                .and_then(|s| toml::from_str::<count::CountStats>(&s).ok())
            {
                Some(s) => s,
                None => count::CountStats {
                    cutoff: 0,
                    peak: None,
                    solid_kmers: kmc::header(&dir.join(name))?.total_kmers,
                    error_kmers: 0,
                    low_coverage: false,
                },
            };
            (b.read(&dir.join(name))?, stats)
        } else {
            let db = tmp.join(format!("classify_{}", pipeline::sanitize(name)));
            let st = match input {
                config::Input::Reads(files) => count::count_reads(
                    files,
                    &db,
                    &tmp.join("count"),
                    threads,
                    memory_gb,
                    b.bundle.min_peak_coverage,
                    min_solid,
                )?,
                config::Input::Assembly(fa) => count::count_assembly(fa, &db)?,
            };
            let r = b.read(&db)?;
            for ext in [".kmc_pre", ".kmc_suf"] {
                let mut p = db.clone().into_os_string();
                p.push(ext);
                fs::remove_file(p).ok();
            }
            (r, st)
        };
        let sc = classify::scores(&b, &present);
        let calls = classify::call(&b.bundle.groups, levels, &sc, st.low_coverage);
        let top = calls.first().and_then(|c| c.call.as_deref());
        let excess = pipeline::kmer_excess(&b.bundle.solid_reference, top, st.solid_kmers);
        let f = |x: Option<f64>| x.map_or("NA".to_string(), |v| format!("{v:.4}"));
        let mut row = vec![
            name.clone(),
            st.peak.map_or("NA".into(), |p| p.to_string()),
            st.low_coverage.to_string(),
            format!("{unexplained:.4}"),
            format!("{excess:.3}"),
        ];
        let mut flags: Vec<String> = Vec::new();
        // at low coverage, sequencing errors that pass the error cutoff also count as unexplained
        // (up to 84% in pure low-coverage test isolates), so the flag is qualified: errors or a
        // foreign genome cannot be told apart without more reads
        // two genomes (mixture, heterozygous hybrid) or contamination. At low coverage missing
        // k-mers deflate it (pure test isolates 0.22-0.67) and error-rich runs can inflate it, so
        // the flag is qualified there, as UNEXPLAINED
        if excess > b.bundle.max_kmer_excess {
            let q = if st.low_coverage { "?" } else { "" };
            flags.push(format!("EXCESS_KMERS{q}={excess:.2}"));
        }
        if unexplained > b.bundle.unexplained_flag {
            let q = if st.low_coverage { "?" } else { "" };
            flags.push(format!("UNEXPLAINED{q}={unexplained:.2}"));
        }
        for c in &calls {
            row.extend([
                c.call.clone().unwrap_or_else(|| "NA".into()),
                f(c.score),
                c.confidence.to_string(),
            ]);
            for fl in &c.flags {
                let tagged = if fl == "LOW_COVERAGE" {
                    fl.clone()
                } else {
                    format!("{}:{fl}", c.level)
                };
                if !flags.contains(&tagged) {
                    flags.push(tagged);
                }
            }
        }
        row.push(if flags.is_empty() {
            "-".into()
        } else {
            flags.join(",")
        });
        if long {
            row.extend(sc.iter().map(|s| f(s.score)));
        }
        writeln!(out, "{}", row.join("\t"))?;
    }
    out.flush()?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_clusters(
    samples: &Path,
    kmc_db: &Path,
    cutoffs: &[f64],
    exact_max: Option<f64>,
    min_positions: f64,
    scale: u64,
    overlap_scale: u64,
    min_allele: u32,
    output: &Path,
) -> Result<()> {
    let t0 = Instant::now();
    let text =
        fs::read_to_string(samples).with_context(|| format!("reading {}", samples.display()))?;
    let mut lines = text
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'));
    let head: Vec<&str> = lines
        .next()
        .context("empty sample sheet")?
        .split('\t')
        .collect();
    let c = head
        .iter()
        .position(|&h| h.trim() == "sample")
        .context("sample sheet needs a 'sample' column")?;
    let names: Vec<String> = lines
        .filter_map(|l| l.split('\t').nth(c).map(|x| x.trim().to_string()))
        .filter(|x| !x.is_empty())
        .collect();
    let dbs: Vec<PathBuf> = names.iter().map(|n| kmc_db.join(n)).collect();
    if cutoffs.iter().any(|x| x.is_nan() || *x < 0.0) {
        bail!("cutoffs must be >= 0");
    }
    let largest = cutoffs.iter().copied().fold(0.0, f64::max);
    let opt = clusters::Options {
        scale,
        overlap_scale,
        min_positions,
        exact_max: exact_max.unwrap_or((4.0 * largest).max(50.0)),
        min_allele: min_allele.max(1),
        batches: 0,
    };
    if opt.exact_max < largest {
        bail!("--exact-max must be at least the largest cutoff");
    }
    eprintln!("{} samples", names.len());
    let r = clusters::run(&dbs, &opt)?;
    eprintln!(
        "{} sites ({} in the screen), {} exact pairs, {} pairs with too little shared genome",
        r.sites,
        r.screen_sites,
        r.pairs.len(),
        r.too_little_overlap
    );
    fs::create_dir_all(output)?;
    let ids: Vec<Vec<usize>> = cutoffs
        .iter()
        .map(|&cut| clusters::single_linkage(names.len(), &r.pairs, cut))
        .collect();
    let mut w = BufWriter::new(File::create(output.join("clusters.tsv"))?);
    write!(w, "sample\tpositions_mb\tsites\tmixed\tmixed_share")?;
    for cut in cutoffs {
        write!(w, "\tcluster_{cut}")?;
    }
    writeln!(w)?;
    for (i, n) in names.iter().enumerate() {
        let s = &r.samples[i];
        let share = s.mixed as f64 / (s.sites + s.mixed).max(1) as f64;
        write!(
            w,
            "{n}\t{:.3}\t{}\t{}\t{share:.5}",
            s.positions / 1e6,
            s.sites,
            s.mixed
        )?;
        for id in &ids {
            write!(w, "\tc{}", id[i] + 1)?;
        }
        writeln!(w)?;
    }
    w.flush()?;
    let mut w = BufWriter::new(File::create(output.join("pairs.tsv"))?);
    writeln!(
        w,
        "sample_a\tsample_b\tpositions_mb\tsites\tdiffs\tsnps_per_mb"
    )?;
    for p in &r.pairs {
        writeln!(
            w,
            "{}\t{}\t{:.3}\t{}\t{}\t{:.3}",
            names[p.i as usize],
            names[p.j as usize],
            p.positions / 1e6,
            p.sites,
            p.diffs,
            p.per_mb
        )?;
    }
    w.flush()?;
    let mut w = BufWriter::new(File::create(output.join("screen_hist.tsv"))?);
    writeln!(w, "snps_per_mb_from\tpairs")?;
    for (edge, n) in &r.hist {
        writeln!(w, "{edge:.4}\t{n}")?;
    }
    w.flush()?;
    for (cut, id) in cutoffs.iter().zip(&ids) {
        let k = id.iter().max().map_or(0, |m| m + 1);
        let multi = (0..k)
            .filter(|&c| id.iter().filter(|&&x| x == c).count() > 1)
            .count();
        eprintln!("cutoff {cut}: {k} clusters, {multi} with more than one sample");
    }
    eprintln!("{:.0} s", t0.elapsed().as_secs_f64());
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_groups(
    samples: &Path,
    kmc_db: &Path,
    filter: &str,
    target_specs: &[String],
    f_in: f64,
    output: &Path,
    write_out: bool,
    batches: u64,
) -> Result<()> {
    let t0 = Instant::now();
    let sheet =
        fs::read_to_string(samples).with_context(|| format!("reading {}", samples.display()))?;
    let mut lines = sheet.lines().filter(|l| !l.is_empty());
    let head: Vec<&str> = lines
        .next()
        .context("empty sample sheet")?
        .split('\t')
        .collect();
    let col = |name: &str| {
        head.iter()
            .position(|&h| h == name)
            .with_context(|| format!("no column {name}"))
    };
    let (fcol, fval) = filter
        .split_once('=')
        .context("--filter must be COLUMN=VALUE")?;
    let (c_run, c_f) = (col("run")?, col(fcol)?);
    let rows: Vec<Vec<&str>> = lines
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .filter(|f| f[c_f] == fval)
        .collect();
    let dbs: Vec<PathBuf> = rows.iter().map(|f| kmc_db.join(f[c_run])).collect();
    let mut targets = Vec::new();
    for spec in target_specs {
        let (name, rest) = spec
            .split_once('=')
            .context("--target must be NAME=COLUMN:VALUE")?;
        let (c, v) = rest
            .split_once(':')
            .context("--target must be NAME=COLUMN:VALUE")?;
        let c = col(c)?;
        let members: Vec<usize> = (0..rows.len()).filter(|&i| rows[i][c] == v).collect();
        if members.len() < 3 {
            eprintln!("{name}: only {} training isolates, skipped", members.len());
            continue;
        }
        targets.push(groups::Target::new(name, members, f_in));
    }
    eprintln!("{} isolates, {} targets", dbs.len(), targets.len());
    fs::create_dir_all(output)?;
    let mk = |t: &groups::Target, what: &str| {
        npy::Writer::create(&output.join(format!("{}.{what}.u64.npy", t.name)))
    };
    let mut strict = targets
        .iter()
        .map(|t| mk(t, "strict"))
        .collect::<Result<Vec<_>>>()?;
    let mut extra = targets
        .iter()
        .map(|t| mk(t, "x.strict"))
        .collect::<Result<Vec<_>>>()?;
    let mut out = if write_out {
        targets
            .iter()
            .map(|t| mk(t, "out"))
            .collect::<Result<Vec<_>>>()?
    } else {
        Vec::new()
    };
    let stats = groups::run(&dbs, &targets, write_out, batches, |b| {
        for (t, v) in b.strict.iter().enumerate() {
            strict[t].write(v)?;
            extra[t].write(&b.extra[t])?;
            if write_out {
                out[t].write(&b.out[t])?;
            }
        }
        Ok(())
    })?;
    for w in strict.into_iter().chain(extra).chain(out) {
        w.finish()?;
    }
    let mut tsv = BufWriter::new(File::create(output.join("groups.tsv"))?);
    writeln!(
        tsv,
        "target\tn_in\tn_out\tmin_in\tcandidates\tout0\tout_le1\textra\tout_kmers"
    )?;
    for (t, s) in targets.iter().zip(&stats) {
        let line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            t.name,
            t.members.len(),
            dbs.len() - t.members.len(),
            t.need,
            s.candidates,
            s.out0,
            s.out_le1,
            s.extra,
            s.out_kmers
        );
        println!("{line}");
        writeln!(tsv, "{line}")?;
    }
    tsv.flush()?;
    eprintln!("{:.0} s", t0.elapsed().as_secs_f64());
    Ok(())
}

/// Where each run's marker k-mers come from.
enum Present {
    Npy(PathBuf),
    Kmc(PathBuf),
}

fn run_score(
    groups: &Path,
    samples: &Path,
    present: &Present,
    targets: &[String],
    prefix: &str,
) -> Result<()> {
    let load = |p: PathBuf| npy::read_u64(&p);
    let mut m = HashMap::new();
    let mut mx = HashMap::new();
    for t in targets {
        m.insert(t.as_str(), load(groups.join(format!("{t}.ed2.u64.npy")))?);
        mx.insert(t.as_str(), load(groups.join(format!("{t}.x.ed2.u64.npy")))?);
    }
    // every marker of every target, for intersecting with a run's database
    let mut all: Vec<u64> = m.values().chain(mx.values()).flatten().copied().collect();
    all.sort_unstable();
    all.dedup();
    let sheet =
        fs::read_to_string(samples).with_context(|| format!("reading {}", samples.display()))?;
    let mut lines = sheet.lines();
    let head: Vec<&str> = lines
        .next()
        .context("empty sample sheet")?
        .split('\t')
        .collect();
    let col = |name: &str| {
        head.iter()
            .position(|&h| h == name)
            .with_context(|| format!("no column {name}"))
    };
    let (c_run, c_called, c_major, c_use, c_cov) = (
        col("run")?,
        col("called")?,
        col("major")?,
        col("use")?,
        col("peak_cov")?,
    );

    let mut scores = BufWriter::new(File::create(format!("{prefix}scores.tsv"))?);
    let mut loo = BufWriter::new(File::create(format!("{prefix}loo_sensitivity.tsv"))?);
    writeln!(scores, "run\tcalled\tuse\tpeak_cov\t{}", targets.join("\t"))?;
    writeln!(
        loo,
        "run\ttarget\tcalled\tpeak_cov\tmarkers\tpresent\tsensitivity"
    )?;
    for line in lines.filter(|l| !l.is_empty()) {
        let f: Vec<&str> = line.split('\t').collect();
        let (run, called, major, use_, cov) =
            (f[c_run], f[c_called], f[c_major], f[c_use], f[c_cov]);
        let p = match present {
            Present::Npy(d) => load(d.join(format!("{run}.npy")))?,
            Present::Kmc(d) => kmc::read_present(&d.join(run), &all)?,
        };
        if p.windows(2).any(|w| w[0] >= w[1]) {
            bail!("{run}: present k-mers not sorted and unique");
        }
        let mut row = vec![
            run.to_string(),
            called.to_string(),
            use_.to_string(),
            cov.to_string(),
        ];
        for t in targets {
            let member = use_ == "train" && (major == t || called == t);
            let s = score::score(&m[t.as_str()], &mx[t.as_str()], &p, member);
            if member {
                writeln!(
                    loo,
                    "{run}\t{t}\t{called}\t{cov}\t{}\t{}\t{}",
                    s.markers,
                    s.present,
                    round4(s.score)
                )?;
            }
            row.push(round4(s.score));
        }
        writeln!(scores, "{}", row.join("\t"))?;
    }
    scores.flush()?;
    loo.flush()?;
    Ok(())
}

/// Like pandas `round(4)` then default float printing (0.9574, 0.8, 0.0).
fn round4(x: f64) -> String {
    let r = (x * 1e4).round() / 1e4;
    let s = format!("{r}");
    if s.contains('.') { s } else { format!("{s}.0") }
}
