//! Unitigs: maximal non-branching paths through the node-centric de Bruijn graph of a set of
//! canonical 31-mers (port of `unitigs()` in `kmers/place_markers.py`).
//!
//! Each oriented k-mer is tracked as (fw, rc), so extending by one base is O(1). Every k-mer lands
//! in exactly one unitig. Starts are taken in the order of the input (sorted, so the output is
//! deterministic); the set of unitigs does not depend on it, except where a start falls inside an
//! isolated cycle, which is cut at that start.

use crate::kmer::{K, MASK, decode, revcomp};
use std::collections::HashSet;

const TOP: usize = 2 * (K - 1);

pub struct Unitig {
    pub seq: String,
    pub n_kmers: usize,
}

pub fn unitigs(kmers: &[u64]) -> Vec<Unitig> {
    let kset: HashSet<u64> = kmers.iter().copied().collect();
    let mut seen: HashSet<u64> = HashSet::with_capacity(kset.len());

    let succ = |fw: u64, rc: u64| -> Vec<(u64, u64)> {
        (0..4u64)
            .filter_map(|b| {
                let f = ((fw << 2) & MASK) | b;
                let r = (rc >> 2) | ((3 - b) << TOP);
                kset.contains(&f.min(r)).then_some((f, r))
            })
            .collect()
    };
    // predecessors of fw = successors of rc, reverse-complemented
    let n_pred = |fw: u64, rc: u64| succ(rc, fw).len();

    let extend = |mut fw: u64, mut rc: u64, seen: &mut HashSet<u64>| -> Vec<u8> {
        let mut bases = Vec::new();
        loop {
            let s = succ(fw, rc);
            if s.len() != 1 {
                break;
            }
            let (f, r) = s[0];
            if seen.contains(&f.min(r)) || n_pred(f, r) != 1 {
                break;
            }
            seen.insert(f.min(r));
            bases.push((f & 3) as u8);
            (fw, rc) = (f, r);
        }
        bases
    };

    let mut out = Vec::new();
    for &x in kmers {
        if !seen.insert(x) {
            continue;
        }
        let (fw, rc) = (x, revcomp(x));
        let right = extend(fw, rc, &mut seen);
        let left = extend(rc, fw, &mut seen); // extending the reverse strand = growing leftwards
        let mut seq = String::with_capacity(K + left.len() + right.len());
        seq.extend(left.iter().rev().map(|&b| b"ACGT"[3 - b as usize] as char));
        seq.push_str(&decode(fw, K));
        seq.extend(right.iter().map(|&b| b"ACGT"[b as usize] as char));
        out.push(Unitig {
            seq,
            n_kmers: 1 + left.len() + right.len(),
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::{canonical, encode, revcomp_str};
    use crate::testutil::Rng;

    fn kmers_of(seq: &str) -> Vec<u64> {
        (0..=seq.len() - K)
            .map(|i| canonical(encode(&seq.as_bytes()[i..i + K]).unwrap()))
            .collect()
    }

    #[test]
    fn reconstructs_separate_sequences() {
        let mut rng = Rng::new(3);
        let seqs: Vec<String> = (0..50)
            .map(|_| {
                (0..40 + rng.below(200))
                    .map(|_| b"ACGT"[rng.below(4) as usize] as char)
                    .collect()
            })
            .collect();
        let mut km: Vec<u64> = seqs.iter().flat_map(|s| kmers_of(s)).collect();
        km.sort_unstable();
        km.dedup();
        let u = unitigs(&km);
        assert_eq!(u.iter().map(|x| x.n_kmers).sum::<usize>(), km.len());
        let norm = |s: &str| s.to_string().min(revcomp_str(s));
        let mut got: Vec<String> = u.iter().map(|x| norm(&x.seq)).collect();
        let mut want: Vec<String> = seqs.iter().map(|s| norm(s)).collect();
        got.sort();
        want.sort();
        assert_eq!(got, want);
    }

    #[test]
    fn splits_at_branches() {
        // two sequences sharing a 40-base middle: 3 + 2 unitigs around the branch points
        let mut rng = Rng::new(9);
        let mut r = |n: u64| -> String {
            (0..n)
                .map(|_| b"ACGT"[rng.below(4) as usize] as char)
                .collect()
        };
        let (a1, a2, mid, b1, b2) = (r(49), r(49), r(40), r(49), r(49));
        // the flanks must differ right next to the shared middle, or the branch moves by a base
        let s1 = format!("{a1}A{mid}C{b1}");
        let s2 = format!("{a2}G{mid}T{b2}");
        let mut km: Vec<u64> = kmers_of(&s1).into_iter().chain(kmers_of(&s2)).collect();
        km.sort_unstable();
        km.dedup();
        let u = unitigs(&km);
        assert_eq!(u.iter().map(|x| x.n_kmers).sum::<usize>(), km.len());
        assert_eq!(u.len(), 5);
        let lens: Vec<usize> = u.iter().map(|x| x.seq.len()).collect();
        assert!(
            u.iter().any(|x| x.seq == mid || x.seq == revcomp_str(&mid)),
            "lengths {lens:?}"
        );
    }
}
