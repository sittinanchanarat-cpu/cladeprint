//! Group-specific k-mers in one streaming pass over the per-isolate presence databases (port of
//! `kmers/group_kmers.py` and the candidate step of `kmers/loo_extra.sh`).
//!
//! For each target group T, every member isolate i has a weight w_i (1 unless clade weighting is
//! on: then 1 / size of i's genotype cluster within T, so a clonal cohort counts as one isolate).
//! With W_T = sum of member weights and need_T = ceil(f_in * W_T) (unweighted) or f_in * W_T
//! (weighted, to a small tolerance):
//!   in(k)  = summed weight of members holding k, out(k) = number of non-members holding k
//!   strict = in >= need and out <= max_out   (markers before the near-miss filters)
//!   extra  = extra_lo <= in < need and out <= max_out (for leave-one-out sensitivity; extra_lo is
//!            need - 1 unweighted, the threshold without the heaviest member weighted)
//!   out    = out > max_out                   (out-group set, input of the near-miss filters)
//! plus the counts of the prototype's summary table: candidates (in >= need), out0 (in >= need and
//! out = 0) and out_le1 (in >= need and out <= 1).
//!
//! Child constraints (`children`, hierarchy-aware prevalence): a marker of a parent group must
//! also reach f_in inside every child group that was built (a C. neoformans marker must be in >= 90%
//! of VNI, of VNII and of VNB). Without this, a large child outvotes a small one: in the census VNII
//! is 19 of 2,700 C. neoformans training runs, so a third of the species markers were absent from
//! every VNII isolate. Children are given as indices into the same target slice, so their weights
//! are already counted; a child that was skipped (too few training samples) constrains nothing.
//!
//! Ignored isolates (`ignore`): training isolates whose label at the target's level is unknown but
//! whose known labels match the target's ancestors (a C. neoformans run without a lineage, for a
//! VNI target) count neither in nor out: they may well belong to the group, so holding its k-mers
//! is no evidence against a marker.
//!
//! Out-group tolerance (`max_out`, default 0) is a count of out-group isolates, not a share: a
//! share of a large out-group would tolerate a whole small species. It is applied per k-mer, also
//! in the near-miss filters (an out-group k-mer held by <= max_out isolates blocks nothing), so it
//! survives a few mislabelled or contaminated out-group runs. Out-group counts are never weighted:
//! a whole clonal cluster carrying a k-mer is real cross-reaction, not one mislabel.
//!
//! The databases are globally sorted (KMC 1 layout), so all of them are read in step, one batch of
//! k-mer space at a time (k-mers whose first few bases fall in a range). Memory is bounded
//! by the batch, not by the number of isolates; each database is read once, sequentially.

use crate::kmc::SortedReader;
use crate::kmer::mix;
use anyhow::{Result, ensure};
use rayon::prelude::*;
use std::path::PathBuf;

/// k-mer space is cut by the first bases (at most this many; fewer if a database's LUT prefix is
/// shorter): batches are ranges of these prefixes
const MAX_PART_BASES: u32 = 7;

/// Weighted thresholds are compared with this slack, so float rounding in the sums never decides.
const EPS: f64 = 1e-9;

/// Threshold on the in-group weight for a group of total weight `total`: ceil(f_in * total) when
/// unweighted (the integer rule, same float arithmetic as `score::need_for`), else f_in * total.
pub fn threshold(f_in: f64, total: f64, weighted: bool) -> f64 {
    if weighted {
        f_in * total - EPS
    } else {
        (f_in * total).ceil()
    }
}

#[derive(Clone, Debug)]
pub struct Target {
    pub name: String,
    /// isolates (indices into the database list) in the target group
    pub members: Vec<usize>,
    /// weight of each member (parallel to `members`)
    pub weights: Vec<f64>,
    pub f_in: f64,
    /// some weight differs from 1
    pub weighted: bool,
    /// summed member weight (= number of members when unweighted)
    pub total: f64,
    /// in-group weight a strict k-mer needs
    pub need: f64,
    /// lowest in-group weight kept as a leave-one-out extra
    pub extra_lo: f64,
    /// out-group isolates allowed to hold a marker
    pub max_out: u32,
    /// isolates counted neither in nor out (label unknown at this level)
    pub ignore: Vec<usize>,
    /// indices into the same target slice: child groups a marker must also reach f_in within
    pub children: Vec<usize>,
}

impl Target {
    /// Unweighted target with no out-group tolerance (the prototype's rule).
    pub fn new(name: &str, members: Vec<usize>, f_in: f64) -> Self {
        let w = vec![1.0; members.len()];
        Self::weighted(name, members, w, f_in, 0)
    }

    /// Target with member weights (parallel to `members`, each in (0, 1]) and out-group tolerance.
    pub fn weighted(
        name: &str,
        members: Vec<usize>,
        weights: Vec<f64>,
        f_in: f64,
        max_out: u32,
    ) -> Self {
        assert_eq!(members.len(), weights.len(), "one weight per member");
        assert!(
            weights.iter().all(|&w| w > 0.0 && w <= 1.0),
            "weights must be in (0, 1]"
        );
        let weighted = weights.iter().any(|&w| w != 1.0);
        let total: f64 = weights.iter().sum();
        let need = threshold(f_in, total, weighted);
        let extra_lo = if weighted {
            let wmax = weights.iter().copied().fold(0.0, f64::max);
            threshold(f_in, total - wmax, true)
        } else {
            need - 1.0
        };
        Target {
            name: name.to_string(),
            members,
            weights,
            f_in,
            weighted,
            total,
            need,
            extra_lo,
            max_out,
            ignore: Vec::new(),
            children: Vec::new(),
        }
    }

    /// The same target with child groups (indices into the same target slice) whose own `need` a
    /// marker must also reach.
    pub fn with_children(mut self, children: Vec<usize>) -> Self {
        self.children = children;
        self
    }

    /// The same target with isolates that count neither in nor out (must not be members).
    pub fn with_ignore(mut self, ignore: Vec<usize>) -> Self {
        assert!(
            ignore.iter().all(|i| !self.members.contains(i)),
            "an ignored isolate cannot be a member"
        );
        self.ignore = ignore;
        self
    }

    /// Threshold of the group with member `j` (index into `members`) held out. Weights of the
    /// others are kept as they are (holding out one member of a cluster does not reweight its
    /// clone-mates; they have the same k-mers, so this barely matters).
    pub fn need_without(&self, j: usize) -> f64 {
        threshold(self.f_in, self.total - self.weights[j], self.weighted)
    }
}

/// Clade weights: 1 / (members of the same cluster), for members given by their cluster labels
/// (None = a cluster of its own).
pub fn cluster_weights(clusters: &[Option<&str>]) -> Vec<f64> {
    let mut size = std::collections::HashMap::new();
    for c in clusters.iter().flatten() {
        *size.entry(*c).or_insert(0usize) += 1;
    }
    clusters
        .iter()
        .map(|c| c.map_or(1.0, |c| 1.0 / size[c] as f64))
        .collect()
}

#[derive(Default, Clone, Debug, PartialEq)]
pub struct Stats {
    pub candidates: u64,
    pub out0: u64,
    /// in >= need and out <= max_out (= out0 without tolerance)
    pub strict: u64,
    pub out_le1: u64,
    pub extra: u64,
    pub out_kmers: u64,
}

/// One batch of results, in k-mer order; `out[t]` is empty unless out-group sets were requested.
pub struct Batch {
    pub strict: Vec<Vec<u64>>,
    /// per strict k-mer, `1 + children.len()` weights: the in-group weight (the member count when
    /// unweighted) then each child's, for leave-one-out
    pub strict_count: Vec<Vec<f64>>,
    pub extra: Vec<Vec<u64>>,
    /// the same `1 + children.len()` weights per extra
    pub extra_count: Vec<Vec<f64>>,
    pub out: Vec<Vec<u64>>,
    /// hash-sampled k-mers held by any database (only when a sketch threshold is given)
    pub sketch: Vec<u64>,
}

/// Per database: sampled k-mers it holds, and how many of them no other database holds (its
/// leave-one-out "unexplained" k-mers).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SketchCounts {
    pub total: Vec<u64>,
    pub unique: Vec<u64>,
}

/// Runs the pass. `batches`: how many pieces k-mer space is cut into (memory ~ 1/batches of all
/// records; 0 = automatic). `sink` receives each batch in order.
pub fn run(
    dbs: &[PathBuf],
    targets: &[Target],
    want_out: bool,
    batches: u64,
    sink: impl FnMut(Batch) -> Result<()>,
) -> Result<Vec<Stats>> {
    Ok(run_sketch(dbs, targets, want_out, batches, None, sink)?.0)
}

/// As `run`, and with `sketch = Some(threshold)` also the hash-sampled k-mer set of all databases
/// (in each batch) and per-database sample counts.
pub fn run_sketch(
    dbs: &[PathBuf],
    targets: &[Target],
    want_out: bool,
    batches: u64,
    sketch: Option<u64>,
    mut sink: impl FnMut(Batch) -> Result<()>,
) -> Result<(Vec<Stats>, SketchCounts)> {
    ensure!(targets.len() <= 32, "at most 32 targets");
    for tg in targets {
        for &c in &tg.children {
            ensure!(
                c < targets.len(),
                "target {}: child index {c} out of range",
                tg.name
            );
            ensure!(
                !targets[c].children.contains(&c),
                "target {}: child {} is its own child",
                tg.name,
                targets[c].name
            );
        }
    }
    ensure!(dbs.len() < u32::MAX as usize, "too many databases");
    let mut mask = vec![0u32; dbs.len()];
    // ign[i]: targets for which isolate i counts neither in nor out
    let mut ign = vec![0u32; dbs.len()];
    for (t, tg) in targets.iter().enumerate() {
        for &i in &tg.ignore {
            ensure!(
                i < dbs.len(),
                "target {}: ignored {i} out of range",
                tg.name
            );
            ign[i] |= 1 << t;
        }
    }
    // weight[i * nt + t]: weight of isolate i in target t (0 if not a member)
    let mut weight = vec![0f64; dbs.len() * targets.len()];
    for (t, tg) in targets.iter().enumerate() {
        for (&i, &w) in tg.members.iter().zip(&tg.weights) {
            ensure!(i < dbs.len(), "target {}: member {i} out of range", tg.name);
            mask[i] |= 1 << t;
            weight[i * targets.len() + t] = w;
        }
    }
    ensure!(!dbs.is_empty(), "no databases");
    let mut readers = dbs
        .iter()
        .map(|d| SortedReader::open(d))
        .collect::<Result<Vec<_>>>()?;
    let k = readers[0].header().k;
    ensure!(
        readers.iter().all(|r| r.header().k == k),
        "databases differ in k"
    );
    ensure!(k <= 31, "group counting needs k <= 31");
    let part_bases = readers
        .iter()
        .map(|r| r.header().lut_prefix_len)
        .min()
        .unwrap()
        .min(MAX_PART_BASES);
    let parts = 1u64 << (2 * part_bases);
    // batches = 0: enough that one batch holds ~32M records (~256 MB of k-mers)
    let total: u64 = readers.iter().map(|r| r.header().total_kmers).sum();
    let batches = if batches == 0 {
        total.div_ceil(32 << 20)
    } else {
        batches
    }
    .clamp(1, parts);
    let mut stats = vec![Stats::default(); targets.len()];
    let nt = targets.len();
    let mut sk = SketchCounts {
        total: vec![0; dbs.len()],
        unique: vec![0; dbs.len()],
    };

    for b in 0..batches {
        let end = parts * (b + 1) / batches;
        // read this batch from every database, in parallel
        let lists: Vec<Vec<u64>> = readers
            .par_iter_mut()
            .map(|r| {
                let mut v = Vec::new();
                r.read_until(part_bases, end, &mut v).map(|_| v)
            })
            .collect::<Result<_>>()?;
        // count in parallel chunks of k-mer space, each a k-way merge of the (sorted) lists
        let chunks = split_points(&lists, 8 * rayon::current_num_threads());
        let parts: Vec<(Batch, Vec<Stats>, SketchCounts)> = chunks
            .par_windows(2)
            .map(|w| {
                let slices: Vec<&[u64]> = lists
                    .iter()
                    .map(|l| &l[l.partition_point(|&x| x < w[0])..l.partition_point(|&x| x < w[1])])
                    .collect();
                count_chunk(&slices, &mask, &ign, &weight, targets, want_out, sketch)
            })
            .collect();
        drop(lists);
        let mut batch = Batch {
            strict: vec![Vec::new(); nt],
            strict_count: vec![Vec::new(); nt],
            extra: vec![Vec::new(); nt],
            extra_count: vec![Vec::new(); nt],
            out: vec![Vec::new(); nt],
            sketch: Vec::new(),
        };
        for (b, st, sc) in parts {
            batch.sketch.extend(&b.sketch);
            for (i, (t, u)) in sc.total.iter().zip(&sc.unique).enumerate() {
                sk.total[i] += t;
                sk.unique[i] += u;
            }
            for t in 0..nt {
                batch.strict[t].extend(&b.strict[t]);
                batch.strict_count[t].extend(&b.strict_count[t]);
                batch.extra[t].extend(&b.extra[t]);
                batch.extra_count[t].extend(&b.extra_count[t]);
                batch.out[t].extend(&b.out[t]);
                let s = &mut stats[t];
                s.candidates += st[t].candidates;
                s.out0 += st[t].out0;
                s.strict += st[t].strict;
                s.out_le1 += st[t].out_le1;
                s.extra += st[t].extra;
                s.out_kmers += st[t].out_kmers;
            }
        }
        sink(batch)?;
    }
    Ok((stats, sk))
}

/// Chunk boundaries for `lists` (each sorted): quantiles of the longest list, from 0 to
/// u64::MAX, which is above every k-mer for k <= 31 (checked in `run`).
fn split_points(lists: &[Vec<u64>], n: usize) -> Vec<u64> {
    let longest = lists
        .iter()
        .max_by_key(|l| l.len())
        .map(Vec::as_slice)
        .unwrap_or(&[]);
    let mut b = vec![0u64];
    for c in 1..n {
        if let Some(&x) = longest.get(longest.len() * c / n)
            && x > *b.last().unwrap()
        {
            b.push(x);
        }
    }
    b.push(u64::MAX);
    b
}

/// Counts one chunk: a k-way merge over the per-isolate slices.
fn count_chunk(
    slices: &[&[u64]],
    mask: &[u32],
    ign: &[u32],
    weight: &[f64],
    targets: &[Target],
    want_out: bool,
    sketch: Option<u64>,
) -> (Batch, Vec<Stats>, SketchCounts) {
    use std::cmp::Reverse;
    use std::collections::BinaryHeap;
    let nt = targets.len();
    let mut batch = Batch {
        strict: vec![Vec::new(); nt],
        strict_count: vec![Vec::new(); nt],
        extra: vec![Vec::new(); nt],
        extra_count: vec![Vec::new(); nt],
        out: vec![Vec::new(); nt],
        sketch: Vec::new(),
    };
    let mut stats = vec![Stats::default(); nt];
    let mut sk = SketchCounts {
        total: vec![0; slices.len()],
        unique: vec![0; slices.len()],
    };
    let mut holders: Vec<u32> = Vec::new();
    let mut pos = vec![0usize; slices.len()];
    let mut heap: BinaryHeap<Reverse<(u64, u32)>> = slices
        .iter()
        .enumerate()
        .filter(|(_, s)| !s.is_empty())
        .map(|(i, s)| Reverse((s[0], i as u32)))
        .collect();
    let mut inc = vec![0u32; nt];
    let mut ignc = vec![0u32; nt];
    let mut inw = vec![0f64; nt];
    while let Some(&Reverse((kmer, _))) = heap.peek() {
        let mut total = 0u32;
        inc.iter_mut().for_each(|c| *c = 0);
        ignc.iter_mut().for_each(|c| *c = 0);
        inw.iter_mut().for_each(|c| *c = 0.0);
        holders.clear();
        while let Some(&Reverse((k, i))) = heap.peek() {
            if k != kmer {
                break;
            }
            heap.pop();
            holders.push(i);
            let i = i as usize;
            let mut m = mask[i];
            while m != 0 {
                let t = m.trailing_zeros() as usize;
                inc[t] += 1;
                inw[t] += weight[i * nt + t];
                m &= m - 1;
            }
            let mut m = ign[i];
            while m != 0 {
                ignc[m.trailing_zeros() as usize] += 1;
                m &= m - 1;
            }
            total += 1;
            pos[i] += 1;
            if let Some(&next) = slices[i].get(pos[i]) {
                heap.push(Reverse((next, i as u32)));
            }
        }
        if let Some(thr) = sketch
            && mix(kmer) < thr
        {
            batch.sketch.push(kmer);
            for &h in &holders {
                sk.total[h as usize] += 1;
            }
            if holders.len() == 1 {
                sk.unique[holders[0] as usize] += 1;
            }
        }
        for (t, tg) in targets.iter().enumerate() {
            let (i, o) = (inw[t], total - inc[t] - ignc[t]);
            let s = &mut stats[t];
            let tolerated = o <= tg.max_out;
            // hierarchy-aware prevalence: every built child must reach its own threshold, and for
            // leave-one-out extras its lowest possible threshold (without its heaviest member)
            let child_ok = tg.children.iter().all(|&c| inw[c] >= targets[c].need);
            let child_ok_x = tg.children.iter().all(|&c| inw[c] >= targets[c].extra_lo);
            if i >= tg.need {
                s.candidates += 1;
                if o == 0 {
                    s.out0 += 1;
                }
                if o <= 1 {
                    s.out_le1 += 1;
                }
            }
            if tolerated {
                let push = |v: &mut Vec<f64>| {
                    v.push(i);
                    for &c in &tg.children {
                        v.push(inw[c]);
                    }
                };
                if i >= tg.need && child_ok {
                    s.strict += 1;
                    batch.strict[t].push(kmer);
                    push(&mut batch.strict_count[t]);
                } else if i >= tg.extra_lo && child_ok_x {
                    s.extra += 1;
                    batch.extra[t].push(kmer);
                    push(&mut batch.extra_count[t]);
                }
            } else {
                s.out_kmers += 1;
                if want_out {
                    batch.out[t].push(kmer);
                }
            }
        }
    }
    (batch, stats, sk)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn need_rounds_up() {
        assert_eq!(Target::new("a", (0..37).collect(), 0.9).need, 34.0);
        assert_eq!(Target::new("a", (0..8).collect(), 0.9).need, 8.0);
        assert_eq!(Target::new("a", (0..10).collect(), 0.9).need, 9.0);
        let t = Target::new("a", (0..37).collect(), 0.9);
        assert_eq!((t.extra_lo, t.need_without(0)), (33.0, 33.0));
    }

    #[test]
    fn cluster_weights_and_weighted_thresholds() {
        // a clone of 4, a pair, and 2 singletons (one unlabelled): weight 1 per cluster
        let c = [
            Some("x"),
            Some("x"),
            Some("x"),
            Some("x"),
            Some("y"),
            Some("y"),
            Some("z"),
            None,
        ];
        let w = cluster_weights(&c);
        assert_eq!(w, vec![0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 1.0, 1.0]);
        let t = Target::weighted("a", (0..8).collect(), w, 0.9, 0);
        assert!(t.weighted && (t.total - 4.0).abs() < 1e-12);
        assert!((t.need - 3.6).abs() < 1e-6);
        // without a singleton (weight 1): 0.9 * 3 = 2.7; extras start there
        assert!((t.need_without(7) - 2.7).abs() < 1e-6 && (t.extra_lo - 2.7).abs() < 1e-6);
        // exact boundary despite float sums: 10 members of weight 0.1 at f_in 1 need all ten
        let t = Target::weighted("b", (0..10).collect(), vec![0.1; 10], 1.0, 0);
        let sum10: f64 = (0..10).map(|_| 0.1).sum();
        let sum9: f64 = (0..9).map(|_| 0.1).sum();
        assert!(sum10 >= t.need && sum9 < t.need);
        // all weights 1 is the unweighted rule
        let t = Target::weighted("c", (0..37).collect(), vec![1.0; 37], 0.9, 2);
        assert!(!t.weighted && t.need == 34.0 && t.max_out == 2);
    }
}
