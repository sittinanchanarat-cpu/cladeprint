//! Per-sample presence databases (KMC 1 layout): reads through KMC with an error cutoff from the
//! sample's own k-mer histogram (as `kmers/build_dbs.sh`), assemblies counted natively (every
//! canonical 31-mer, no error filter).

use crate::kmc;
use crate::kmer::{K, canonical, encode};
use anyhow::{Context, Result, bail, ensure};
use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct CountStats {
    /// presence threshold (count >= cutoff); 1 for assemblies
    pub cutoff: u32,
    /// genome-peak coverage (0 = no peak found; None for assemblies)
    pub peak: Option<u32>,
    pub solid_kmers: u64,
    pub error_kmers: u64,
    pub low_coverage: bool,
}

/// Error cutoff and genome peak from a k-mer count histogram (index = count, up to 1000), as the
/// prototype: the cutoff is the first local minimum after count 1 (at least 2); if the histogram
/// only falls up to 20 there is no genome peak (cutoff 2, low coverage). The peak is the highest
/// bin at or above the cutoff; below `min_peak` is low coverage too. Below about 1x the genome
/// k-mers sit at counts 1-3 and a dip in the repeat tail (10-20) passes for the valley, which
/// would drop the genome as errors and take a repeat bump for the peak: fewer than `min_solid`
/// k-mers at or above the cutoff (0 = no check) is treated as no genome peak.
pub fn cutoff_from_histogram(h: &[u64], min_peak: u32, min_solid: u64) -> CountStats {
    let at = |i: usize| h.get(i).copied().unwrap_or(0);
    let mut c = 2;
    while c < 20 && at(c + 1) < at(c) {
        c += 1;
    }
    let no_valley = c >= 20 || (c..h.len()).map(at).sum::<u64>() < min_solid;
    if no_valley {
        c = 2;
    }
    let mut p = c;
    for i in c..h.len() {
        if at(i) > at(p) {
            p = i;
        }
    }
    let peak = if no_valley { 0 } else { p as u32 };
    CountStats {
        cutoff: c as u32,
        peak: Some(peak),
        solid_kmers: (c..h.len()).map(at).sum(),
        error_kmers: (1..c).map(at).sum(),
        low_coverage: no_valley || peak < min_peak,
    }
}

fn run(cmd: &mut Command) -> Result<()> {
    let out = cmd
        .output()
        .with_context(|| format!("running {:?} (is it installed?)", cmd.get_program()))?;
    if !out.status.success() {
        bail!(
            "{:?} failed: {}{}",
            cmd.get_program(),
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
    }
    Ok(())
}

/// Reads (FASTQ, gzipped or not) -> presence database `db` via KMC. `tmp` must be on disk.
pub fn count_reads(
    files: &[PathBuf],
    db: &Path,
    tmp: &Path,
    threads: usize,
    memory_gb: u32,
    min_peak: u32,
    min_solid: u64,
) -> Result<CountStats> {
    for f in files {
        ensure!(f.exists(), "reads file {} not found", f.display());
    }
    std::fs::create_dir_all(tmp)?;
    let name = db
        .file_name()
        .context("database path")?
        .to_string_lossy()
        .to_string();
    let list = tmp.join(format!("{name}.files"));
    std::fs::write(
        &list,
        files
            .iter()
            .map(|f| format!("{}\n", f.display()))
            .collect::<String>(),
    )?;
    let all = tmp.join(format!("{name}.all"));
    let work = tmp.join(format!("{name}.kmc"));
    std::fs::create_dir_all(&work)?;
    let threads = if threads == 0 {
        rayon::current_num_threads()
    } else {
        threads
    };
    run(Command::new("kmc")
        .arg(format!("-k{K}"))
        .arg(format!("-m{}", memory_gb.max(2)))
        .arg(format!("-t{threads}"))
        .args(["-ci1", "-cs1000"])
        .arg(format!("@{}", list.display()))
        .arg(&all)
        .arg(&work))?;
    let hist = tmp.join(format!("{name}.hist"));
    run(Command::new("kmc_tools")
        .args(["transform"])
        .arg(&all)
        .arg("histogram")
        .arg(&hist)
        .arg("-cx1000"))?;
    let mut h = vec![0u64; 1001];
    for line in std::fs::read_to_string(&hist)?.lines() {
        let mut f = line.split_whitespace();
        if let (Some(c), Some(n)) = (f.next(), f.next()) {
            let c: usize = c.parse()?;
            if c < h.len() {
                h[c] = n.parse()?;
            }
        }
    }
    let stats = cutoff_from_histogram(&h, min_peak, min_solid);
    run(Command::new("kmc_tools")
        .arg("transform")
        .arg(&all)
        .arg(format!("-ci{}", stats.cutoff))
        .arg("set_counts")
        .arg("1")
        .arg(db))?;
    for ext in [".kmc_pre", ".kmc_suf"] {
        let mut p = all.clone().into_os_string();
        p.push(ext);
        std::fs::remove_file(p).ok();
    }
    std::fs::remove_dir_all(&work).ok();
    std::fs::remove_file(&list).ok();
    // keep the k-mer count histogram next to the database (for later mixture or ploidy checks)
    let mut keep = db.as_os_str().to_owned();
    keep.push(".hist");
    std::fs::rename(&hist, &keep).or_else(|_| std::fs::copy(&hist, &keep).map(|_| ()))?;
    Ok(stats)
}

/// Every canonical 31-mer of an assembly (FASTA, gzipped or not), sorted and unique; k-mers with
/// a non-ACGT base are skipped.
pub fn assembly_kmers(fasta: &Path) -> Result<Vec<u64>> {
    let mut reader = needletail::parse_fastx_file(fasta)
        .with_context(|| format!("reading {}", fasta.display()))?;
    let mut v = Vec::new();
    while let Some(rec) = reader.next() {
        let rec = rec.with_context(|| format!("parsing {}", fasta.display()))?;
        let seq = rec.seq();
        for w in seq.windows(K) {
            if let Some(x) = encode(w) {
                v.push(canonical(x));
            }
        }
    }
    v.sort_unstable();
    v.dedup();
    Ok(v)
}

/// Assembly -> presence database `db` (no error filter).
pub fn count_assembly(fasta: &Path, db: &Path) -> Result<CountStats> {
    let v = assembly_kmers(fasta)?;
    ensure!(!v.is_empty(), "{}: no 31-mers", fasta.display());
    kmc::write_kmc1(db, &v, K as u32, 7)?;
    Ok(CountStats {
        cutoff: 1,
        peak: None,
        solid_kmers: v.len() as u64,
        error_kmers: 0,
        low_coverage: false,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hist(pairs: &[(usize, u64)]) -> Vec<u64> {
        let mut h = vec![0u64; 1001];
        for &(c, n) in pairs {
            h[c] = n;
        }
        h
    }

    #[test]
    fn cutoff_rules_match_the_prototype() {
        // errors at 1-3, valley at 4, genome peak at 15
        let mut h = hist(&[(1, 9000), (2, 2000), (3, 500), (4, 100), (5, 150)]);
        for (c, v) in h.iter_mut().enumerate().take(31).skip(6) {
            *v = 1000 - (c as i64 - 15).unsigned_abs() * 40;
        }
        let s = cutoff_from_histogram(&h, 5, 0);
        assert_eq!((s.cutoff, s.peak, s.low_coverage), (4, Some(15), false));
        assert_eq!(s.error_kmers, 9000 + 2000 + 500);
        // only falling: no genome peak
        let falling: Vec<u64> = (0..1001)
            .map(|c| if c == 0 { 0 } else { 100_000 / c as u64 })
            .collect();
        let s = cutoff_from_histogram(&falling, 5, 0);
        assert_eq!((s.cutoff, s.peak, s.low_coverage), (2, Some(0), true));
        // a peak, but below min_peak
        let low = hist(&[(1, 5000), (2, 900), (3, 1000), (4, 800)]);
        let s = cutoff_from_histogram(&low, 5, 0);
        assert_eq!((s.cutoff, s.peak, s.low_coverage), (2, Some(3), true));
        // below 1x: genome k-mers at 1-3, a dip in the repeat tail at 11, a repeat bump at 12
        // (thinned from a 6x run); only the solid k-mer floor rejects that valley
        let sub1x = hist(&[
            (1, 6_000_000),
            (2, 524_305),
            (3, 61_122),
            (4, 11_700),
            (5, 5_848),
            (6, 4_318),
            (7, 3_146),
            (8, 2_270),
            (9, 1_676),
            (10, 1_303),
            (11, 1_136),
            (12, 1_171),
            (13, 1_171),
            (14, 1_071),
            (15, 1_005),
            (16, 849),
            (17, 721),
            (18, 546),
        ]);
        let s = cutoff_from_histogram(&sub1x, 5, 0);
        assert_eq!((s.cutoff, s.peak, s.low_coverage), (11, Some(12), false));
        let s = cutoff_from_histogram(&sub1x, 5, 1_000_000);
        assert_eq!((s.cutoff, s.peak, s.low_coverage), (2, Some(0), true));
        assert_eq!(s.solid_kmers, sub1x[2..].iter().sum::<u64>());
        // a real valley clears the floor: same result as without it
        let s = cutoff_from_histogram(&h, 5, 10_000);
        assert_eq!((s.cutoff, s.peak, s.low_coverage), (4, Some(15), false));
    }

    #[test]
    fn assembly_counting() {
        let d = std::env::temp_dir().join(format!("cladeprint-count-{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        let fa = d.join("a.fa");
        let s1 = "ACGTTGCAAGGCTTAACCGGTTAAGCTTGCAAGT"; // 34 bp: 4 k-mers
        std::fs::write(&fa, format!(">x\n{s1}\n>y\nACGTNACGT\n>z\n{}\n", &s1[..31])).unwrap();
        let v = assembly_kmers(&fa).unwrap();
        assert_eq!(v.len(), 4, "z repeats the first k-mer of x; y is too short");
        assert!(v.windows(2).all(|w| w[0] < w[1]));
        let db = d.join("a");
        let st = count_assembly(&fa, &db).unwrap();
        assert_eq!(st.solid_kmers, 4);
        assert_eq!(kmc::read_sorted(&db).unwrap(), v);
        std::fs::remove_dir_all(&d).unwrap();
    }
}
