//! Genotype clusters for clade weighting, from split k-mers in the presence databases (the idea of
//! SKA, Harris 2018): no reads, no reference, and low coverage costs sites rather than adding
//! differences.
//!
//! A canonical 31-mer = left flank (bases 0-14) + middle base (15) + right flank (16-30). K-mers
//! with the same flanks are the alleles of one site. Reverse complementing keeps the middle base in
//! the middle, and canonical orientation is decided by the flanks alone unless the left flank is
//! the reverse complement of the right (those k-mers are dropped), so every allele of a site is
//! stored in the same orientation. At each site a sample has one allele (called), none (missing)
//! or several (mixed: a mixture, a heterozygous hybrid or a sequencing error; counted per sample,
//! treated as missing). A site is kept when at least two alleles are each called in some sample.
//! Optionally alleles held by fewer than `min_allele` samples are ignored first.
//!
//! Distance between two samples = differing sites per Mb of genome both cover: a SNP density,
//! which does not shrink as more diverse samples join the set (a share of segregating sites would).
//! Covered genome is estimated from a hash sample of every split k-mer (1/scale per sample; a
//! coarser 1/overlap_scale for each pair's shared genome, which is all a 1% estimate needs and keeps
//! the all-pairs stage cheap). The same hash samples the screened sites, so ratios are unbiased. A SNP is seen only if no other difference lies
//! within 15 bases, so large distances (between species) are underestimated; close pairs, which
//! clustering is about, are not.
//!
//! All pairs are screened on the hash-sampled sites (also giving a histogram for choosing cutoffs);
//! pairs that could be within `exact_max` get exact distances over all sites. Clusters are single
//! linkage at a cutoff.
//!
//! The pass reads every database in step, one range of k-mer space at a time, like
//! `groups::run_sketch`; the alleles of a site share their first 15 bases, so never cross a batch.

use crate::kmc::SortedReader;
use crate::kmer::{K, mix, sketch_threshold};
use anyhow::{Result, ensure};
use rayon::prelude::*;
use std::path::PathBuf;

const MID: usize = K / 2;
const FLANK_BITS: u32 = 2 * MID as u32;
const FLANK_MASK: u64 = (1 << FLANK_BITS) - 1;
/// batches are ranges of at most this many leading bases (fewer if a database's LUT is shorter)
const MAX_PART_BASES: u32 = 7;

/// Reverse complement of a packed 15-base flank.
#[inline]
fn revcomp_flank(x: u64) -> u64 {
    let mut x = !x;
    x = ((x >> 2) & 0x3333_3333_3333_3333) | ((x & 0x3333_3333_3333_3333) << 2);
    x = ((x >> 4) & 0x0F0F_0F0F_0F0F_0F0F) | ((x & 0x0F0F_0F0F_0F0F_0F0F) << 4);
    x.swap_bytes() >> (64 - FLANK_BITS)
}

/// Split form of a canonical 31-mer: left flank << 32 | right flank << 2 | middle base, so that
/// sorting brings the alleles of a site (value >> 2) together. None when the left flank is the
/// reverse complement of the right (the stored orientation would depend on the middle base).
#[inline]
pub fn split(x: u64) -> Option<u64> {
    let left = x >> (FLANK_BITS + 2);
    let mid = (x >> FLANK_BITS) & 3;
    let right = x & FLANK_MASK;
    (left != revcomp_flank(right)).then_some((left << 32) | (right << 2) | mid)
}

#[derive(Debug, Clone)]
pub struct Options {
    /// hash-sampling scale of the screened sites and of each sample's covered genome
    pub scale: u64,
    /// coarser hash-sampling scale (a multiple of `scale`) for each pair's shared covered genome:
    /// ~18k sampled k-mers per genome at 1024, a 1% estimate, and the all-pairs stage's main cost
    pub overlap_scale: u64,
    /// pairs whose shared covered genome is estimated below this get no distance
    pub min_positions: f64,
    /// exact distances for every pair the screen cannot rule out below this (SNPs per Mb)
    pub exact_max: f64,
    /// alleles held by fewer samples are ignored (1 = keep all; 2 drops run-private alleles,
    /// mostly sequencing errors in old error-rich runs, at the cost of private SNPs)
    pub min_allele: u32,
    /// pieces k-mer space is cut into (0 = automatic, ~32M records per batch)
    pub batches: u64,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            scale: 64,
            overlap_scale: 1024,
            min_positions: 1e6,
            exact_max: 200.0,
            min_allele: 1,
            batches: 0,
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct SampleStats {
    /// estimated genome positions covered (split k-mers called, from the hash sample)
    pub positions: f64,
    /// sites called
    pub sites: u64,
    /// sites with more than one allele in this sample
    pub mixed: u64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Pair {
    pub i: u32,
    pub j: u32,
    /// estimated genome positions both cover
    pub positions: f64,
    /// sites both call
    pub sites: u64,
    /// sites where they differ
    pub diffs: u64,
    pub per_mb: f64,
}

#[derive(Debug, Clone)]
pub struct Report {
    pub sites: usize,
    pub screen_sites: usize,
    pub samples: Vec<SampleStats>,
    /// exact distances, i < j, for the pairs the screen kept (covering >= min_positions)
    pub pairs: Vec<Pair>,
    /// screened distances of all pairs covering >= min_positions: (lower bin edge in SNPs per Mb,
    /// pairs); the first bin is exactly 0, then 10 bins per decade from 0.1
    pub hist: Vec<(f64, u64)>,
    /// pairs sharing too little covered genome for a distance
    pub too_little_overlap: u64,
}

/// Genotypes of one sample: per 64 sites, [called, allele bit 0, allele bit 1].
type Plane = Vec<[u64; 3]>;

fn push_site(planes: &mut [Plane], n: usize, calls: &[(u32, u8)]) {
    let bit = 1u64 << (n % 64);
    if n.is_multiple_of(64) {
        planes.iter_mut().for_each(|p| p.push([0; 3]));
    }
    for &(s, a) in calls {
        let w = &mut planes[s as usize][n / 64];
        w[0] |= bit;
        if a & 1 != 0 {
            w[1] |= bit;
        }
        if a & 2 != 0 {
            w[2] |= bit;
        }
    }
}

/// (sites both call, sites where they differ)
#[inline]
fn compare(a: &[[u64; 3]], b: &[[u64; 3]]) -> (u64, u64) {
    let (mut shared, mut diff) = (0u64, 0u64);
    for (x, y) in a.iter().zip(b) {
        let both = x[0] & y[0];
        shared += both.count_ones() as u64;
        diff += (both & ((x[1] ^ y[1]) | (x[2] ^ y[2]))).count_ones() as u64;
    }
    (shared, diff)
}

#[inline]
fn shared_bits(a: &[u64], b: &[u64]) -> u64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x & y).count_ones() as u64)
        .sum()
}

/// Results of one chunk of k-mer space, in order.
#[derive(Default)]
struct Chunk {
    /// sites: calls[site_off[s]..site_off[s + 1]] as (sample, allele)
    site_off: Vec<u32>,
    calls: Vec<(u32, u8)>,
    site_sampled: Vec<bool>,
    /// hash-sampled split k-mers called by >= 2 samples (all, not only sites): samples calling each
    key_off: Vec<u32>,
    key_calls: Vec<u32>,
    stats: Vec<SampleStats>,
}

fn count_chunk(slices: &[&[u64]], thr: u64, thr_overlap: u64, min_allele: u32) -> Chunk {
    use std::cmp::Reverse;
    use std::collections::BinaryHeap;
    let n = slices.len();
    let mut c = Chunk {
        stats: vec![SampleStats::default(); n],
        ..Default::default()
    };
    let mut pos = vec![0usize; n];
    let mut mask = vec![0u8; n];
    let mut touched: Vec<u32> = Vec::new();
    let mut heap: BinaryHeap<Reverse<(u64, u32)>> = slices
        .iter()
        .enumerate()
        .filter(|(_, s)| !s.is_empty())
        .map(|(i, s)| Reverse((s[0], i as u32)))
        .collect();
    while let Some(&Reverse((v, _))) = heap.peek() {
        let key = v >> 2;
        touched.clear();
        while let Some(&Reverse((w, i))) = heap.peek() {
            if w >> 2 != key {
                break;
            }
            heap.pop();
            let iu = i as usize;
            if mask[iu] == 0 {
                touched.push(i);
            }
            mask[iu] |= 1 << (w & 3);
            pos[iu] += 1;
            if let Some(&next) = slices[iu].get(pos[iu]) {
                heap.push(Reverse((next, i)));
            }
        }
        // alleles held by fewer than min_allele samples are ignored (sequencing errors are
        // almost always private to one run)
        let mut support = [0u32; 4];
        for &i in &touched {
            let m = mask[i as usize];
            for (a, s) in support.iter_mut().enumerate() {
                *s += (m >> a & 1) as u32;
            }
        }
        let common = (0..4).fold(0u8, |acc, a| {
            acc | (((support[a] >= min_allele) as u8) << a)
        });
        let mut clean = [0u32; 4];
        for &i in &touched {
            let m = mask[i as usize] & common;
            mask[i as usize] = m;
            if m.count_ones() == 1 {
                clean[m.trailing_zeros() as usize] += 1;
            }
        }
        let h = mix(key);
        let sampled = h < thr;
        if sampled {
            for &i in &touched {
                if mask[i as usize].count_ones() == 1 {
                    c.stats[i as usize].positions += 1.0;
                }
            }
        }
        if h < thr_overlap {
            // only k-mers called by two or more samples can count for a pair's shared genome, so
            // only those are stored (most k-mers held by one run are sequencing errors)
            let start = c.key_calls.len();
            for &i in &touched {
                if mask[i as usize].count_ones() == 1 {
                    c.key_calls.push(i);
                }
            }
            if c.key_calls.len() - start >= 2 {
                c.key_off.push(start as u32);
            } else {
                c.key_calls.truncate(start);
            }
        }
        if clean.iter().filter(|&&x| x > 0).count() >= 2 {
            c.site_off.push(c.calls.len() as u32);
            c.site_sampled.push(sampled);
            for &i in &touched {
                let m = mask[i as usize];
                if m.count_ones() == 1 {
                    c.calls.push((i, m.trailing_zeros() as u8));
                    c.stats[i as usize].sites += 1;
                } else if m != 0 {
                    c.stats[i as usize].mixed += 1;
                }
            }
        }
        for &i in &touched {
            mask[i as usize] = 0;
        }
    }
    c
}

/// Chunk boundaries for sorted split-form lists: quantiles of the longest list, rounded down to a
/// site boundary (multiple of 4), from 0 to u64::MAX (above every split value).
fn split_points(lists: &[Vec<u64>], n: usize) -> Vec<u64> {
    let longest = lists
        .iter()
        .max_by_key(|l| l.len())
        .map(Vec::as_slice)
        .unwrap_or(&[]);
    let mut b = vec![0u64];
    for c in 1..n {
        if let Some(&x) = longest.get(longest.len() * c / n) {
            let x = x & !3;
            if x > *b.last().unwrap() {
                b.push(x);
            }
        }
    }
    b.push(u64::MAX);
    b
}

/// Genotype planes for every sample: all sites, the hash-sampled sites, and the hash-sampled
/// split k-mers each sample calls.
struct Genotypes {
    full: Vec<Plane>,
    screen: Vec<Plane>,
    keys: Vec<Vec<u64>>,
    sites: usize,
    screen_sites: usize,
    n_keys: usize,
    stats: Vec<SampleStats>,
}

fn genotype(dbs: &[PathBuf], opt: &Options) -> Result<Genotypes> {
    ensure!(!dbs.is_empty(), "no databases");
    ensure!(dbs.len() < u32::MAX as usize, "too many databases");
    let mut readers = dbs
        .iter()
        .map(|d| SortedReader::open(d))
        .collect::<Result<Vec<_>>>()?;
    ensure!(
        readers.iter().all(|r| r.header().k == K as u32),
        "clusters needs k = 31 databases"
    );
    let part_bases = readers
        .iter()
        .map(|r| r.header().lut_prefix_len)
        .min()
        .unwrap()
        .min(MAX_PART_BASES);
    let parts = 1u64 << (2 * part_bases);
    let total: u64 = readers.iter().map(|r| r.header().total_kmers).sum();
    let batches = if opt.batches == 0 {
        total.div_ceil(32 << 20)
    } else {
        opt.batches
    }
    .clamp(1, parts);
    ensure!(
        opt.scale >= 1
            && opt.overlap_scale >= opt.scale
            && opt.overlap_scale.is_multiple_of(opt.scale),
        "overlap_scale must be a multiple of scale"
    );
    let thr = sketch_threshold(opt.scale);
    let thr_overlap = sketch_threshold(opt.overlap_scale);
    let n = dbs.len();
    let mut g = Genotypes {
        full: vec![Vec::new(); n],
        screen: vec![Vec::new(); n],
        keys: vec![Vec::new(); n],
        sites: 0,
        screen_sites: 0,
        n_keys: 0,
        stats: vec![SampleStats::default(); n],
    };
    for b in 0..batches {
        let end = parts * (b + 1) / batches;
        let lists: Vec<Vec<u64>> = readers
            .par_iter_mut()
            .map(|r| {
                let mut v = Vec::new();
                r.read_until(part_bases, end, &mut v)?;
                let mut s: Vec<u64> = v.into_iter().filter_map(split).collect();
                s.sort_unstable();
                Ok(s)
            })
            .collect::<Result<_>>()?;
        let bounds = split_points(&lists, 8 * rayon::current_num_threads());
        let chunks: Vec<Chunk> = bounds
            .par_windows(2)
            .map(|w| {
                let slices: Vec<&[u64]> = lists
                    .iter()
                    .map(|l| &l[l.partition_point(|&x| x < w[0])..l.partition_point(|&x| x < w[1])])
                    .collect();
                count_chunk(&slices, thr, thr_overlap, opt.min_allele)
            })
            .collect();
        drop(lists);
        for c in chunks {
            for (s, t) in g.stats.iter_mut().zip(&c.stats) {
                s.positions += t.positions;
                s.sites += t.sites;
                s.mixed += t.mixed;
            }
            for s in 0..c.site_off.len() {
                let end = c.site_off.get(s + 1).map_or(c.calls.len(), |&e| e as usize);
                let calls = &c.calls[c.site_off[s] as usize..end];
                push_site(&mut g.full, g.sites, calls);
                g.sites += 1;
                if c.site_sampled[s] {
                    push_site(&mut g.screen, g.screen_sites, calls);
                    g.screen_sites += 1;
                }
            }
            for k in 0..c.key_off.len() {
                let end = c
                    .key_off
                    .get(k + 1)
                    .map_or(c.key_calls.len(), |&e| e as usize);
                if g.n_keys.is_multiple_of(64) {
                    g.keys.iter_mut().for_each(|p| p.push(0));
                }
                for &s in &c.key_calls[c.key_off[k] as usize..end] {
                    g.keys[s as usize][g.n_keys / 64] |= 1 << (g.n_keys % 64);
                }
                g.n_keys += 1;
            }
        }
    }
    // planes of samples never touched by a later push are short: pad (all-zero words)
    let (wf, ws, wk) = (
        g.sites.div_ceil(64),
        g.screen_sites.div_ceil(64),
        g.n_keys.div_ceil(64),
    );
    for i in 0..n {
        g.full[i].resize(wf, [0; 3]);
        g.screen[i].resize(ws, [0; 3]);
        g.keys[i].resize(wk, 0);
        g.stats[i].positions *= opt.scale as f64;
    }
    Ok(g)
}

/// Histogram bin of a distance in SNPs per Mb: 0 for exactly 0, then 10 per decade from 0.1.
fn bin(per_mb: f64) -> usize {
    if per_mb <= 0.0 {
        0
    } else {
        1 + ((per_mb.log10() + 1.0) * 10.0).floor().max(0.0) as usize
    }
}

fn bin_edge(b: usize) -> f64 {
    if b == 0 {
        0.0
    } else {
        10f64.powf((b - 1) as f64 / 10.0 - 1.0)
    }
}

/// Genotypes every database, screens all pairs, and computes exact distances for close ones.
pub fn run(dbs: &[PathBuf], opt: &Options) -> Result<Report> {
    let g = genotype(dbs, opt)?;
    let n = dbs.len();
    let bytes =
        n * (24 * (g.sites.div_ceil(64) + g.screen_sites.div_ceil(64)) + 8 * g.n_keys.div_ceil(64));
    eprintln!(
        "[clusters] {} sites, {} screened sites, {} shared sampled k-mers: {:.2} GB of genotypes",
        g.sites,
        g.screen_sites,
        g.n_keys,
        bytes as f64 / 1e9
    );
    let (scale, oscale) = (opt.scale as f64, opt.overlap_scale as f64);
    let total_pairs = n * n.saturating_sub(1) / 2;
    let done = std::sync::atomic::AtomicUsize::new(0);
    let t0 = std::time::Instant::now();
    // per row: (exact pairs, histogram, pairs with too little overlap)
    let rows: Vec<(Vec<Pair>, Vec<u64>, u64)> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut pairs = Vec::new();
            let mut hist = Vec::new();
            let mut little = 0u64;
            for j in i + 1..n {
                let positions = shared_bits(&g.keys[i], &g.keys[j]) as f64 * oscale;
                if positions < opt.min_positions {
                    little += 1;
                    continue;
                }
                let (_, ds) = compare(&g.screen[i], &g.screen[j]);
                let b = bin(ds as f64 * scale / positions * 1e6);
                if hist.len() <= b {
                    hist.resize(b + 1, 0);
                }
                hist[b] += 1;
                // screen: sampled differences expected at exact_max, with generous slack
                let e = opt.exact_max * positions / 1e6 / scale;
                if (ds as f64) <= 2.0 * e + 4.0 {
                    let (sites, diffs) = compare(&g.full[i], &g.full[j]);
                    pairs.push(Pair {
                        i: i as u32,
                        j: j as u32,
                        positions,
                        sites,
                        diffs,
                        per_mb: diffs as f64 / positions * 1e6,
                    });
                }
            }
            // progress, every ~10% of the pairs
            let before = done.fetch_add(n - 1 - i, std::sync::atomic::Ordering::Relaxed);
            let tenth = total_pairs.div_ceil(10).max(1);
            if (before + n - 1 - i) / tenth > before / tenth {
                eprintln!(
                    "[clusters] pairs {:.0}% ({:.0} s)",
                    100.0 * (before + n - 1 - i) as f64 / total_pairs as f64,
                    t0.elapsed().as_secs_f64()
                );
            }
            (pairs, hist, little)
        })
        .collect();
    let mut pairs = Vec::new();
    let mut hist: Vec<u64> = Vec::new();
    let mut too_little_overlap = 0;
    for (p, h, l) in rows {
        pairs.extend(p);
        if hist.len() < h.len() {
            hist.resize(h.len(), 0);
        }
        for (a, b) in hist.iter_mut().zip(h) {
            *a += b;
        }
        too_little_overlap += l;
    }
    Ok(Report {
        sites: g.sites,
        screen_sites: g.screen_sites,
        samples: g.stats,
        pairs,
        hist: hist
            .into_iter()
            .enumerate()
            .map(|(b, c)| (bin_edge(b), c))
            .collect(),
        too_little_overlap,
    })
}

/// Single-linkage clusters of `n` samples joining pairs at or below `cutoff` SNPs per Mb. Cluster
/// numbers start at 0 and go by size (largest first), then by first member.
pub fn single_linkage(n: usize, pairs: &[Pair], cutoff: f64) -> Vec<usize> {
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    let mut parent: Vec<usize> = (0..n).collect();
    for p in pairs.iter().filter(|p| p.per_mb <= cutoff) {
        let (a, b) = (
            find(&mut parent, p.i as usize),
            find(&mut parent, p.j as usize),
        );
        if a != b {
            parent[a.max(b)] = a.min(b);
        }
    }
    let roots: Vec<usize> = (0..n).map(|i| find(&mut parent, i)).collect();
    let mut size = vec![0usize; n];
    for &r in &roots {
        size[r] += 1;
    }
    let mut order: Vec<usize> = (0..n).filter(|&r| size[r] > 0).collect();
    order.sort_by_key(|&r| (std::cmp::Reverse(size[r]), r));
    let mut id = vec![0usize; n];
    for (k, &r) in order.iter().enumerate() {
        id[r] = k;
    }
    roots.iter().map(|&r| id[r]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::canonical;
    use crate::testutil::Rng;

    #[test]
    fn alleles_of_a_site_share_flanks_in_both_orientations() {
        let mut rng = Rng::new(3);
        let mut seen = 0;
        for _ in 0..20_000 {
            let x = rng.next_u64() & crate::kmer::MASK;
            let sib: Vec<u64> = (0..4u64)
                .map(|m| (x & !(3 << FLANK_BITS)) | (m << FLANK_BITS))
                .collect();
            // whichever strand each allele is stored on, they land on one site key
            let fwd: Vec<Option<u64>> = sib.iter().map(|&s| split(canonical(s))).collect();
            if let Some(k0) = fwd[0] {
                seen += 1;
                for f in &fwd {
                    assert_eq!(f.unwrap() >> 2, k0 >> 2, "alleles share the site key");
                }
                let mids: std::collections::HashSet<u64> =
                    fwd.iter().map(|f| f.unwrap() & 3).collect();
                assert_eq!(mids.len(), 4, "four distinct alleles");
            }
        }
        assert!(seen > 19_900);
        // left flank = revcomp(right flank): dropped
        let right = rng.next_u64() & FLANK_MASK;
        let x = (revcomp_flank(right) << 32) | (2 << FLANK_BITS) | right;
        assert_eq!(split(x), None);
    }

    #[test]
    fn linkage_and_bins() {
        let p = |i, j, d| Pair {
            i,
            j,
            positions: 1e6,
            sites: 100,
            diffs: d,
            per_mb: d as f64,
        };
        // 0-1-2 chained at 5, 3-4 at 20, 5 alone
        let pairs = [p(0, 1, 5), p(1, 2, 5), p(3, 4, 20), p(0, 3, 500)];
        assert_eq!(single_linkage(6, &pairs, 10.0), vec![0, 0, 0, 1, 2, 3]);
        assert_eq!(single_linkage(6, &pairs, 50.0), vec![0, 0, 0, 1, 1, 2]);
        assert_eq!(single_linkage(6, &pairs, 1000.0), vec![0, 0, 0, 0, 0, 1]);
        assert_eq!((bin(0.0), bin(0.1), bin(1.0), bin(10.0)), (0, 1, 11, 21));
        assert!((bin_edge(21) - 10.0).abs() < 1e-9 && bin(9.99) == 20);
    }
}
