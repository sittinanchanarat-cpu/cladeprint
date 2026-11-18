//! 31-mers packed 2 bits per base (A=0, C=1, G=2, T=3), first base in the most significant bits.
//! Same encoding as the Python prototype (`kmers/near_miss.py`), so `.npy` files interchange.

pub const K: usize = 31;
pub const MASK: u64 = (1 << (2 * K)) - 1;
/// The low bit of every 2-bit base.
pub const EVEN: u64 = 0x5555_5555_5555_5555 & MASK;

/// Reverse complement of a packed 31-mer.
#[inline]
pub fn revcomp(x: u64) -> u64 {
    // complement (3 - b = !b on 2 bits), then reverse the 32 two-bit groups of the word
    let mut x = !x;
    x = ((x >> 2) & 0x3333_3333_3333_3333) | ((x & 0x3333_3333_3333_3333) << 2);
    x = ((x >> 4) & 0x0F0F_0F0F_0F0F_0F0F) | ((x & 0x0F0F_0F0F_0F0F_0F0F) << 4);
    x.swap_bytes() >> (64 - 2 * K)
}

#[inline]
pub fn canonical(x: u64) -> u64 {
    x.min(revcomp(x))
}

/// Number of bases that differ between two packed 31-mers.
#[inline]
pub fn hamming(a: u64, b: u64) -> u32 {
    let x = a ^ b;
    ((x | (x >> 1)) & EVEN).count_ones()
}

/// Base `i` (0 = first) of a packed k-mer of length `n`.
#[inline]
pub fn base(x: u64, n: usize, i: usize) -> u8 {
    ((x >> (2 * (n - 1 - i))) & 3) as u8
}

pub fn encode(seq: &[u8]) -> Option<u64> {
    let mut v = 0u64;
    for &c in seq {
        let b = match c {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => return None,
        };
        v = (v << 2) | b;
    }
    Some(v)
}

pub fn decode(x: u64, n: usize) -> String {
    (0..n)
        .map(|i| b"ACGT"[base(x, n, i) as usize] as char)
        .collect()
}

/// Reverse complement of a DNA string (ACGT only).
pub fn revcomp_str(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|c| match c {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}

/// 64-bit mix of a packed k-mer (splitmix64 finalizer), for hash sampling of k-mer sets.
#[inline]
pub fn mix(x: u64) -> u64 {
    let mut z = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Hash-sampling threshold keeping about 1/`scale` of k-mers (FracMinHash style): a k-mer is in
/// the sample when `mix(kmer) < threshold`. The same k-mers are kept in every set, so shares
/// measured on the samples estimate shares of the full sets.
pub fn sketch_threshold(scale: u64) -> u64 {
    u64::MAX / scale.max(1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::Rng;

    fn slow_revcomp(x: u64) -> u64 {
        let mut r = 0;
        let mut x = x ^ MASK;
        for _ in 0..K {
            r = (r << 2) | (x & 3);
            x >>= 2;
        }
        r
    }

    #[test]
    fn revcomp_matches_definition() {
        let mut rng = Rng::new(7);
        for _ in 0..10_000 {
            let x = rng.next_u64() & MASK;
            assert_eq!(revcomp(x), slow_revcomp(x));
            assert_eq!(revcomp(revcomp(x)), x);
        }
        let s = "ACGTTGCAAAAACCCCCGGGGGTTTTTACGTA";
        let s = &s[..K];
        assert_eq!(
            decode(revcomp(encode(s.as_bytes()).unwrap()), K),
            revcomp_str(s)
        );
    }

    #[test]
    fn sketch_keeps_about_one_in_scale() {
        let mut rng = Rng::new(3);
        let thr = sketch_threshold(64);
        let kept = (0..640_000)
            .filter(|_| mix(rng.next_u64() & MASK) < thr)
            .count();
        assert!((9_000..11_000).contains(&kept), "{kept}");
    }

    #[test]
    fn hamming_counts_bases() {
        let a = encode(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").unwrap();
        let b = encode(b"TAAAAAAAAAAAAAAAAAAAAAAAAAAAACA").unwrap();
        assert_eq!(hamming(a, b), 2);
        assert_eq!(hamming(a, a), 0);
        assert_eq!(hamming(0, MASK), 31);
    }
}
