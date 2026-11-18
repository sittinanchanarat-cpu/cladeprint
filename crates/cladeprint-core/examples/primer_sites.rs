//! Benchmark helper (R5, published-assay audit): best binding site of each primer/probe in each
//! genome, with up to MAX_MM mismatches.
//!
//! Each oligo (IUPAC allowed) is expanded to its concrete sequences, then every variant within
//! MAX_MM substitutions goes into a hash table keyed by the 2-bit encoded L-mer. Each genome FASTA
//! (here: BCALM2 unitigs of a presence database; unitigs overlap by 30 bases, so every site of
//! length <= 31 lies inside one) is streamed once per oligo length, both strands.
//!
//! Usage: primer_sites <oligos.tsv> <genomes.tsv> <out.tsv> [--all]
//!   oligos.tsv   name <tab> sequence (5'->3', IUPAC), header line skipped
//!   genomes.tsv  sample <tab> fasta path
//!   out.tsv      sample, oligo, mismatches (0..MAX_MM, or "none"), mm_3prime (last 5 bases),
//!                n_sites at the best mismatch count, contig, pos (0-based, genome strand), strand
//!   --all        every site within MAX_MM instead of the best one (for in-silico PCR on
//!                assemblies): sample, oligo, mismatches, mm_3prime, tail (last 5 bases 5'->3',
//!                '.' match, 'x' mismatch), contig, pos, strand

use anyhow::{Context, Result, bail};
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;

const MAX_MM: u8 = 3;
const TAIL: usize = 5;
const MAX_EXPANSION: usize = 4096;

fn code(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn iupac(b: u8) -> Result<&'static [u8]> {
    Ok(match b.to_ascii_uppercase() {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' | b'U' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"CG",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' => b"ACGT",
        c => bail!("not an IUPAC base: {}", c as char),
    })
}

fn expand(seq: &[u8]) -> Result<Vec<Vec<u8>>> {
    let mut out: Vec<Vec<u8>> = vec![Vec::new()];
    for &b in seq {
        let alts = iupac(b)?;
        out = out
            .iter()
            .flat_map(|p| {
                alts.iter().map(move |&a| {
                    let mut q = p.clone();
                    q.push(a);
                    q
                })
            })
            .collect();
        if out.len() > MAX_EXPANSION {
            bail!("oligo expands to more than {MAX_EXPANSION} sequences");
        }
    }
    Ok(out)
}

fn enc(s: &[u8]) -> u64 {
    s.iter().fold(0, |v, &b| (v << 2) | code(b).unwrap())
}

/// (mismatches, mismatches in the 3' tail, tail mask) of each variant, keyed by encoded L-mer;
/// keeps the best. Mask bit j = mismatch j bases from the 3' end (bit 0: the terminal base).
type Table = HashMap<u64, (u8, u8, u8)>;

/// One length's merged table: encoded L-mer -> (oligo index, mismatches, 3' mismatches, mask).
type Merged = HashMap<u64, Vec<(usize, u8, u8, u8)>>;

fn add_variants(t: &mut Table, s: &[u8]) {
    let l = s.len();
    let mut rec = |v: &[u8], mm: u8| {
        let mask = (0..TAIL.min(l))
            .filter(|&j| v[l - 1 - j] != s[l - 1 - j])
            .fold(0u8, |m, j| m | (1 << j));
        let tail = mask.count_ones() as u8;
        let e = t.entry(enc(v)).or_insert((mm, tail, mask));
        if (mm, tail, mask) < *e {
            *e = (mm, tail, mask);
        }
    };
    let bases = *b"ACGT";
    let mut v = s.to_vec();
    fn go(
        v: &mut Vec<u8>,
        s: &[u8],
        start: usize,
        mm: u8,
        bases: &[u8; 4],
        rec: &mut dyn FnMut(&[u8], u8),
    ) {
        rec(v, mm);
        if mm == MAX_MM {
            return;
        }
        for i in start..s.len() {
            for &b in bases {
                if b != s[i] {
                    v[i] = b;
                    go(v, s, i + 1, mm + 1, bases, rec);
                    v[i] = s[i];
                }
            }
        }
    }
    go(&mut v, s, 0, 0, &bases, &mut rec);
}

struct Oligo {
    name: String,
    len: usize,
    table: Table,
}

#[derive(Clone)]
struct Hit {
    mm: u8,
    tail: u8,
    n: u32,
    contig: String,
    pos: usize,
    strand: char,
}

/// All oligos of one length in one table: encoded L-mer -> (oligo index, mismatches, 3' mismatches),
/// so each genome position costs one lookup per strand whatever the number of oligos.
fn merged(oligos: &[Oligo]) -> Vec<(usize, Merged)> {
    let mut by_len: HashMap<usize, Merged> = HashMap::new();
    for (oi, o) in oligos.iter().enumerate() {
        let m = by_len.entry(o.len).or_default();
        for (&k, &(mm, tail, mask)) in &o.table {
            m.entry(k).or_default().push((oi, mm, tail, mask));
        }
    }
    let mut v: Vec<_> = by_len.into_iter().collect();
    v.sort_by_key(|(l, _)| *l);
    v
}

/// Every site (with `all`): oligo index, mismatches, tail mask, contig, pos, strand.
type Site = (usize, u8, u8, String, usize, char);

fn scan(
    fasta: &str,
    n_oligos: usize,
    tables: &[(usize, Merged)],
    all: bool,
) -> Result<(Vec<Option<Hit>>, Vec<Site>)> {
    let mut best: Vec<Option<Hit>> = vec![None; n_oligos];
    let mut sites = Vec::new();
    let mut r = parse_fastx_file(fasta).with_context(|| format!("reading {fasta}"))?;
    while let Some(rec) = r.next() {
        let rec = rec?;
        let seq = rec.seq();
        let id = String::from_utf8_lossy(rec.id())
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        for (l, table) in tables {
            let l = *l;
            let mask = if l == 32 {
                u64::MAX
            } else {
                (1u64 << (2 * l)) - 1
            };
            let (mut fw, mut rc, mut valid) = (0u64, 0u64, 0usize);
            for (i, &b) in seq.iter().enumerate() {
                let Some(c) = code(b) else {
                    valid = 0;
                    continue;
                };
                fw = ((fw << 2) | c) & mask;
                rc = (rc >> 2) | ((3 - c) << (2 * (l - 1)));
                valid += 1;
                if valid < l {
                    continue;
                }
                let start = i + 1 - l;
                for (key, strand) in [(fw, '+'), (rc, '-')] {
                    let Some(entries) = table.get(&key) else {
                        continue;
                    };
                    for &(oi, mm, tail, mask) in entries {
                        if all {
                            sites.push((oi, mm, mask, id.clone(), start, strand));
                        }
                        {
                            let h = &mut best[oi];
                            match h {
                                Some(b) if (mm, tail) == (b.mm, b.tail) => b.n += 1,
                                Some(b) if (mm, tail) > (b.mm, b.tail) => {}
                                _ => {
                                    *h = Some(Hit {
                                        mm,
                                        tail,
                                        n: 1,
                                        contig: id.clone(),
                                        pos: start,
                                        strand,
                                    })
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    Ok((best, sites))
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let all = a.len() == 5 && a[4] == "--all";
    if a.len() != 4 && !all {
        bail!("usage: primer_sites <oligos.tsv> <genomes.tsv> <out.tsv> [--all]");
    }
    let mut oligos = Vec::new();
    for line in std::fs::read_to_string(&a[1])?.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 2 || f[1].trim().is_empty() {
            continue;
        }
        let seq = f[1].trim().as_bytes();
        if seq.len() > 32 {
            bail!("{}: {} bases; at most 32 supported", f[0], seq.len());
        }
        let mut table = Table::new();
        for s in expand(seq).with_context(|| f[0].to_string())? {
            add_variants(&mut table, &s);
        }
        oligos.push(Oligo {
            name: f[0].to_string(),
            len: seq.len(),
            table,
        });
    }
    let genomes: Vec<(String, String)> = std::fs::read_to_string(&a[2])?
        .lines()
        .filter_map(|l| {
            l.split_once('\t')
                .map(|(s, p)| (s.to_string(), p.to_string()))
        })
        .collect();
    eprintln!("{} oligos, {} genomes", oligos.len(), genomes.len());
    let tables = merged(&oligos);
    let rows: Vec<Result<String>> = genomes
        .par_iter()
        .map(|(s, p)| {
            let (best, sites) = scan(p, oligos.len(), &tables, all)?;
            let mut out = String::new();
            if all {
                for (oi, mm, mask, contig, pos, strand) in sites {
                    let o = &oligos[oi];
                    let n = TAIL.min(o.len);
                    let tail: String = (0..n)
                        .rev()
                        .map(|j| if mask & (1 << j) != 0 { 'x' } else { '.' })
                        .collect();
                    out += &format!(
                        "{s}\t{}\t{mm}\t{}\t{tail}\t{contig}\t{pos}\t{strand}\n",
                        o.name,
                        mask.count_ones()
                    );
                }
                return Ok(out);
            }
            for (o, h) in oligos.iter().zip(best) {
                match h {
                    Some(h) => {
                        out += &format!(
                            "{s}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                            o.name, h.mm, h.tail, h.n, h.contig, h.pos, h.strand
                        )
                    }
                    None => out += &format!("{s}\t{}\tnone\t\t0\t\t\t\n", o.name),
                }
            }
            Ok(out)
        })
        .collect();
    let mut w = std::io::BufWriter::new(std::fs::File::create(&a[3])?);
    if all {
        writeln!(
            w,
            "sample\toligo\tmismatches\tmm_3prime\ttail\tcontig\tpos\tstrand"
        )?;
    } else {
        writeln!(
            w,
            "sample\toligo\tmismatches\tmm_3prime\tn_sites\tcontig\tpos\tstrand"
        )?;
    }
    for r in rows {
        w.write_all(r?.as_bytes())?;
    }
    Ok(())
}
