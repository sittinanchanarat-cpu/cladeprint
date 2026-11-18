//! Random spot subsampling of FASTQ streams, deterministic per spot name.
//!
//! A spot (one read, or a read pair) is kept when a 64-bit hash of its name and a seed falls below
//! `fraction * 2^64`. The decision depends only on the spot itself, so a run decoded in chunks (in
//! parallel, e.g. `fastq-dump -N/-X`) gives the same sample as the whole run in one stream; the
//! prototype's awk sampler drew one random stream in order and could not be split.
//! Mates are the consecutive records whose names agree after dropping a trailing "/1" or "/2"
//! (fastq-dump --split-spot writes both mates under the same name).

use anyhow::{Result, bail};
use std::io::{BufRead, Write};

/// FNV-1a over the bytes, then the splitmix64 finalizer, seeded.
pub fn spot_hash(name: &[u8], seed: u64) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64 ^ seed.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    for &b in name {
        h ^= b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01B3);
    }
    let mut z = h.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Spot key of a FASTQ header line: first token without '@' and without a trailing /1 or /2.
fn spot_key(header: &[u8]) -> &[u8] {
    let h = header.strip_prefix(b"@").unwrap_or(header);
    let end = h
        .iter()
        .position(|c| c.is_ascii_whitespace())
        .unwrap_or(h.len());
    let t = &h[..end];
    if t.len() > 2 && t[t.len() - 2] == b'/' && matches!(t[t.len() - 1], b'1' | b'2') {
        &t[..t.len() - 2]
    } else {
        t
    }
}

pub struct Counts {
    pub spots: u64,
    pub kept_spots: u64,
    pub records: u64,
}

/// Copies the kept spots from `input` (FASTQ) to `output`.
pub fn subsample(
    input: impl BufRead,
    mut output: impl Write,
    fraction: f64,
    seed: u64,
) -> Result<Counts> {
    if !(0.0..=1.0).contains(&fraction) {
        bail!("fraction must be in [0, 1]");
    }
    let threshold = if fraction >= 1.0 {
        u64::MAX
    } else {
        (fraction * 2f64.powi(64)) as u64
    };
    let mut input = input;
    let mut c = Counts {
        spots: 0,
        kept_spots: 0,
        records: 0,
    };
    let mut last: Vec<u8> = Vec::new();
    let mut keep = false;
    // four reused line buffers (header, sequence, '+', qualities), newlines included
    let mut rec: [Vec<u8>; 4] = Default::default();
    loop {
        rec[0].clear();
        if input.read_until(b'\n', &mut rec[0])? == 0 {
            break;
        }
        if rec[0].trim_ascii().is_empty() {
            continue;
        }
        if rec[0][0] != b'@' {
            bail!("record {}: header does not start with '@'", c.records + 1);
        }
        for r in rec[1..].iter_mut() {
            r.clear();
            if input.read_until(b'\n', r)? == 0 {
                bail!("truncated FASTQ record {}", c.records + 1);
            }
        }
        c.records += 1;
        let key = spot_key(rec[0].trim_ascii_end());
        if key != last.as_slice() {
            last.clear();
            last.extend_from_slice(key);
            c.spots += 1;
            keep = fraction >= 1.0 || spot_hash(key, seed) < threshold;
            c.kept_spots += keep as u64;
        }
        if keep {
            for r in &rec {
                output.write_all(r)?;
                if r.last() != Some(&b'\n') {
                    output.write_all(b"\n")?;
                }
            }
        }
    }
    output.flush()?;
    Ok(c)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fastq(spots: usize, paired: bool, slash: bool) -> Vec<u8> {
        let mut s = String::new();
        for i in 1..=spots {
            for m in 1..=(1 + paired as usize) {
                let name = if slash {
                    format!("SRR1.{i}/{m}")
                } else {
                    format!("SRR1.{i} HWI:1:{i} length=4")
                };
                s += &format!("@{name}\nACGT\n+\nIIII\n");
            }
        }
        s.into_bytes()
    }

    fn run(input: &[u8], f: f64) -> (Vec<u8>, Counts) {
        let mut out = Vec::new();
        let c = subsample(input, &mut out, f, 20261012).unwrap();
        (out, c)
    }

    #[test]
    fn keeps_mates_together_and_hits_the_fraction() {
        for slash in [false, true] {
            let input = fastq(20_000, true, slash);
            let (out, c) = run(&input, 0.1);
            assert_eq!((c.spots, c.records), (20_000, 40_000));
            assert!((1_800..2_200).contains(&c.kept_spots), "{}", c.kept_spots);
            let text = String::from_utf8(out).unwrap();
            let names: Vec<&str> = text
                .lines()
                .step_by(4)
                .map(|h| spot_key(h.as_bytes()))
                .map(|k| std::str::from_utf8(k).unwrap())
                .collect();
            assert_eq!(names.len() as u64, 2 * c.kept_spots);
            assert!(
                names.chunks(2).all(|p| p[0] == p[1]),
                "a mate was split from its pair"
            );
        }
    }

    #[test]
    fn chunked_equals_whole() {
        let input = fastq(9_999, true, false);
        let (whole, _) = run(&input, 0.05);
        // split at spot boundaries (8 lines per spot), as fastq-dump -N/-X chunks would
        let lines: Vec<&[u8]> = input.split(|&b| b == b'\n').collect();
        let mut joined = Vec::new();
        for chunk in lines.chunks(8 * 1_234) {
            let mut part = chunk.join(&b'\n');
            part.push(b'\n');
            joined.extend(run(&part, 0.05).0);
        }
        assert_eq!(joined, whole);
    }

    #[test]
    fn edge_fractions_and_errors() {
        let input = fastq(100, false, false);
        assert_eq!(run(&input, 1.0).1.kept_spots, 100);
        assert_eq!(run(&input, 0.0).1.kept_spots, 0);
        assert!(subsample(&b"@a\nAC\n+\n"[..], Vec::new(), 0.5, 1).is_err());
        assert!(subsample(&b"a\nAC\n+\nII\n"[..], Vec::new(), 0.5, 1).is_err());
    }
}
