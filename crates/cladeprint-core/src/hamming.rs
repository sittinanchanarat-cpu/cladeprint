//! Hamming near-miss exclusion: drop candidate markers that some out-group k-mer matches with
//! <= 2 mismatches (port of `kmers/near_miss.py`).
//!
//! Candidates and the out-group are canonical 31-mers. A genomic near-neighbour y of candidate c is
//! stored as y or rc(y), so both c and rc(c) are queried. Pigeonhole: split the 31-mer into blocks
//! of 10/10/11 bases; within 2 mismatches, at least one block matches exactly.

use crate::index::{BlockIndex, atomics, into_plain, lower, stream};
use crate::kmer::{K, hamming, revcomp};

/// Distances are capped here: `FAR` means > 2.
pub const FAR: u8 = 3;
/// (bit shift, bases) of the three blocks: bases 0-9, 10-19, 20-30
const BLOCKS: [(u32, usize); 3] = [(42, 10), (22, 10), (0, 11)];
const _: () = assert!(K == 31);

/// Minimum Hamming distance (0-2, or `FAR`) from each query to any out-group k-mer.
pub fn min_dist(queries: &[u64], out: &[u64]) -> Vec<u8> {
    let best = atomics(queries.len(), FAR);
    for &(shift, n) in &BLOCKS {
        let m = (1u64 << (2 * n)) - 1;
        let key = |v: u64| ((v >> shift) & m) as u32;
        let idx = BlockIndex::new(queries.iter().map(|&q| key(q)), n);
        stream(out, |o| {
            for &qi in idx.bucket(key(o)) {
                let d = hamming(queries[qi as usize], o);
                if d < FAR as u32 {
                    lower(&best, qi, d as u8);
                }
            }
        });
    }
    into_plain(best)
}

/// Distance of each candidate, over both orientations.
pub fn candidate_dist(cand: &[u64], out: &[u64]) -> Vec<u8> {
    let q: Vec<u64> = cand
        .iter()
        .copied()
        .chain(cand.iter().map(|&c| revcomp(c)))
        .collect();
    let d = min_dist(&q, out);
    (0..cand.len())
        .map(|i| d[i].min(d[i + cand.len()]))
        .collect()
}

/// Brute force: enumerate every <= 2-mismatch variant of c and rc(c) and look it up in the sorted
/// canonical out-group (the `--verify` check of the prototype).
pub fn brute(c: u64, out_sorted: &[u64]) -> u8 {
    let mut best = FAR;
    for q in [c, revcomp(c)] {
        let mut check = |v: u64| {
            let canon = v.min(revcomp(v));
            if out_sorted.binary_search(&canon).is_ok() {
                best = best.min(hamming(q, v) as u8);
            }
        };
        check(q);
        for i in 0..K {
            let s = 2 * (K - 1 - i);
            for a in 1..4u64 {
                let v1 = q ^ (a << s);
                check(v1);
                for j in i + 1..K {
                    let t = 2 * (K - 1 - j);
                    for b in 1..4u64 {
                        check(v1 ^ (b << t));
                    }
                }
            }
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::{MASK, canonical};
    use crate::testutil::{Rng, mutate};

    #[test]
    fn matches_brute_force() {
        let mut rng = Rng::new(11);
        // out-group: random k-mers plus 1- and 2-mismatch variants of some candidates (either strand)
        let cand: Vec<u64> = (0..400).map(|_| canonical(rng.next_u64() & MASK)).collect();
        let mut out: Vec<u64> = (0..20_000)
            .map(|_| canonical(rng.next_u64() & MASK))
            .collect();
        for (i, &c) in cand.iter().enumerate() {
            match i % 4 {
                0 => out.push(canonical(mutate(&mut rng, c, 1))),
                1 => out.push(canonical(revcomp(mutate(&mut rng, c, 2)))),
                2 => out.push(canonical(mutate(&mut rng, c, 3))),
                _ => {}
            }
        }
        out.sort_unstable();
        out.dedup();
        let d = candidate_dist(&cand, &out);
        let mut counts = [0usize; 4];
        for (i, &c) in cand.iter().enumerate() {
            assert_eq!(d[i], brute(c, &out), "candidate {i}");
            counts[d[i] as usize] += 1;
        }
        assert!(
            counts[1] > 50 && counts[2] > 50 && counts[3] > 50,
            "{counts:?}"
        );
    }
}
