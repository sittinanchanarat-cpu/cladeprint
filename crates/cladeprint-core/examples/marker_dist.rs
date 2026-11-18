//! Benchmark helper (R4): distance from every query k-mer to each genome, the common yardstick for
//! comparing marker sets from different tools.
//!
//! For each genome (a presence database) and each query 31-mer: 0 = present, 1-2 = nearest k-mer
//! at that Hamming or overlap edit distance (the near-miss rules of `run`), 3 = farther.
//!
//! Usage: marker_dist <kmers.txt> <genomes.txt> <db_dir> <out.u8>
//!   kmers.txt   one 31-mer per line (any orientation)
//!   genomes.txt one sample name per line (database <db_dir>/<sample>)
//!   out.u8      genomes x k-mers bytes, row-major (row order = genomes.txt)

use anyhow::{Context, Result, bail};
use cladeprint_core::{edit, hamming, kmc, kmer};
use std::io::Write;
use std::path::Path;

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    if a.len() != 5 {
        bail!("usage: marker_dist <kmers.txt> <genomes.txt> <db_dir> <out.u8>");
    }
    let queries: Vec<u64> = std::fs::read_to_string(&a[1])?
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            kmer::encode(l.trim().as_bytes())
                .filter(|_| l.trim().len() == kmer::K)
                .map(kmer::canonical)
                .with_context(|| format!("not a 31-mer: {l}"))
        })
        .collect::<Result<_>>()?;
    let genomes: Vec<String> = std::fs::read_to_string(&a[2])?
        .lines()
        .filter(|l| !l.is_empty())
        .map(str::to_string)
        .collect();
    let mut out = std::io::BufWriter::new(std::fs::File::create(&a[4])?);
    for (i, g) in genomes.iter().enumerate() {
        let db = kmc::read_sorted(&Path::new(&a[3]).join(g))?;
        let mut d = hamming::candidate_dist(&queries, &db);
        // exact presence is Hamming 0; indel neighbours only for those Hamming calls far
        let far: Vec<usize> = (0..queries.len())
            .filter(|&j| d[j] == hamming::FAR)
            .collect();
        let fq: Vec<u64> = far.iter().map(|&j| queries[j]).collect();
        for (&j, e) in far.iter().zip(edit::candidate_dist(&fq, &db)) {
            d[j] = e.min(edit::FAR);
        }
        out.write_all(&d)?;
        eprintln!("[{}/{}] {g}", i + 1, genomes.len());
    }
    out.flush()?;
    Ok(())
}
