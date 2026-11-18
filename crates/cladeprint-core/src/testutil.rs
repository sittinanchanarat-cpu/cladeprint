//! Deterministic helpers for tests (no external RNG dependency).

use crate::kmer::K;

/// splitmix64
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng(seed)
    }
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    pub fn below(&mut self, n: u64) -> u64 {
        self.next_u64() % n
    }
}

/// `n` substitutions at distinct positions of a packed 31-mer.
pub fn mutate(rng: &mut Rng, x: u64, n: usize) -> u64 {
    let mut pos: Vec<usize> = Vec::new();
    while pos.len() < n {
        let p = rng.below(K as u64) as usize;
        if !pos.contains(&p) {
            pos.push(p);
        }
    }
    pos.iter()
        .fold(x, |v, &p| v ^ ((1 + rng.below(3)) << (2 * (K - 1 - p))))
}

pub fn pack(b: &[u8]) -> u64 {
    b.iter().fold(0, |v, &x| (v << 2) | x as u64)
}

/// As the prototype's self-test: a 31-base window `a` of a random genome, and a 31-base window of
/// a copy of the genome with 0-3 random edits (substitution, deletion, insertion), offset 0-2.
/// `unrelated` returns a random `b` instead.
pub fn edited_pair(rng: &mut Rng, unrelated: bool) -> (Vec<u8>, Vec<u8>) {
    let g: Vec<u8> = (0..K + 6).map(|_| rng.below(4) as u8).collect();
    let a = g[3..3 + K].to_vec();
    if unrelated {
        return (a, (0..K).map(|_| rng.below(4) as u8).collect());
    }
    let mut b = g[3..].to_vec();
    for _ in 0..rng.below(4) {
        let p = 3 + rng.below(b.len() as u64 - 6) as usize;
        match rng.below(3) {
            0 => b[p] = (b[p] + 1 + rng.below(3) as u8) % 4,
            1 => {
                b.remove(p);
            }
            _ => b.insert(p, rng.below(4) as u8),
        }
    }
    let o = rng.below(3) as usize;
    b.resize(b.len().max(o + K), 0);
    (a, b[o..o + K].to_vec())
}
