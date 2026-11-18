//! Bucket index of query k-mers by an exact block of bases, shared by the near-miss filters.
//!
//! Pigeonhole search: queries are bucketed by one block (a 4^n table of offsets); the large
//! out-group set is streamed through in parallel, and each out-group k-mer is compared with the
//! queries in the bucket of its own block. Memory: the queries, the table and one array of query
//! indices, independent of the out-group size.

use rayon::prelude::*;
use std::sync::atomic::{AtomicU8, Ordering};

pub struct BlockIndex {
    /// queries in `order[start[key]..start[key + 1]]` share block `key`
    start: Vec<u32>,
    order: Vec<u32>,
}

impl BlockIndex {
    pub fn new(keys: impl Iterator<Item = u32> + Clone, nbases: usize) -> Self {
        let nb = 1usize << (2 * nbases);
        let mut start = vec![0u32; nb + 1];
        for k in keys.clone() {
            start[k as usize + 1] += 1;
        }
        for i in 0..nb {
            start[i + 1] += start[i];
        }
        let mut fill = start.clone();
        let mut order = vec![0u32; start[nb] as usize];
        for (i, k) in keys.enumerate() {
            order[fill[k as usize] as usize] = i as u32;
            fill[k as usize] += 1;
        }
        BlockIndex { start, order }
    }

    #[inline]
    pub fn bucket(&self, key: u32) -> &[u32] {
        &self.order[self.start[key as usize] as usize..self.start[key as usize + 1] as usize]
    }
}

/// Lowers `best[q]` to `d` (atomic min).
#[inline]
pub fn lower(best: &[AtomicU8], q: u32, d: u8) {
    best[q as usize].fetch_min(d, Ordering::Relaxed);
}

/// Streams `out` through `f(out_kmer)` in parallel chunks.
pub fn stream<F: Fn(u64) + Sync>(out: &[u64], f: F) {
    out.par_chunks(1 << 16)
        .for_each(|ch| ch.iter().for_each(|&o| f(o)));
}

pub fn atomics(n: usize, init: u8) -> Vec<AtomicU8> {
    (0..n).map(|_| AtomicU8::new(init)).collect()
}

pub fn into_plain(v: Vec<AtomicU8>) -> Vec<u8> {
    v.into_iter().map(AtomicU8::into_inner).collect()
}
