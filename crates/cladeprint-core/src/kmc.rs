//! Reader for KMC k-mer databases (`.kmc_pre` + `.kmc_suf`), following `kmc_api/kmc_file.cpp`
//! (listing mode). Replaces `kmc_tools transform dump` + text parsing. Both on-disk versions:
//!   - 0x200 (KMC 2/3, written by `kmc`): "KMCP" | LUT (u64 per prefix, cumulative record index,
//!     one block of 4^prefix_len per signature bin, in a row) | guard u64 | signature map
//!     (u32 x (4^sig_len + 1)) | header | header offset | "KMCP". Header: k, mode, counter size,
//!     prefix length, signature length, min count, max count (u32 each), total (u64), strand flag.
//!   - 0 (KMC 1 layout, written by `kmc_tools`): "KMCP" | LUT (4^prefix_len + 1 u64) | header |
//!     header offset | "KMCP". Header: k, mode, counter size, prefix length, min count, max count
//!     low (u32 each), total (u64), strand flag (u8), max count high (u32). Records globally sorted.
//!
//! The version is the u32 12 bytes from the end; the header offset is read as one byte.
//! `.kmc_suf`: "KMCS" | records | "KMCS"; a record is the suffix (2 bits per base, big-endian,
//! (k - prefix length) / 4 bytes) and a little-endian counter of `counter_size` bytes.
//! A k-mer is prefix << (2 * suffix bases) | suffix: the same 2-bit encoding as `kmer` (A=0, C=1,
//! G=2, T=3, first base most significant). Format 0x200 records are sorted only within a signature
//! bin; `read_sorted` sorts.

use anyhow::{Context, Result, bail, ensure};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub struct Header {
    pub k: u32,
    pub counter_size: u32,
    pub lut_prefix_len: u32,
    pub signature_len: u32,
    pub min_count: u32,
    pub max_count: u32,
    pub total_kmers: u64,
}

fn with_ext(db: &Path, ext: &str) -> PathBuf {
    let mut s = db.as_os_str().to_owned();
    s.push(ext);
    PathBuf::from(s)
}

fn u32_at(b: &[u8], i: usize) -> u32 {
    u32::from_le_bytes(b[i..i + 4].try_into().unwrap())
}

/// Opens database `db` (path without extension). Returns the header and the LUT.
fn read_prefix(db: &Path) -> Result<(Header, Vec<u64>)> {
    let path = with_ext(db, ".kmc_pre");
    let buf = std::fs::read(&path).with_context(|| format!("reading {}", path.display()))?;
    let n = buf.len();
    ensure!(
        n >= 16 && &buf[..4] == b"KMCP" && &buf[n - 4..] == b"KMCP",
        "{}: not a KMC prefix file",
        path.display()
    );
    let version = u32_at(&buf, n - 12);
    if version != 0 && version != 0x200 {
        bail!(
            "{}: unsupported KMC database version {version:#x}",
            path.display()
        );
    }
    let header_offset = buf[n - 8] as usize;
    let h = n - 8 - header_offset;
    ensure!(
        u32_at(&buf, h + 4) == 0,
        "{}: Quake-mode counters are not supported",
        path.display()
    );
    let header = if version == 0x200 {
        Header {
            k: u32_at(&buf, h),
            counter_size: u32_at(&buf, h + 8),
            lut_prefix_len: u32_at(&buf, h + 12),
            signature_len: u32_at(&buf, h + 16),
            min_count: u32_at(&buf, h + 20),
            max_count: u32_at(&buf, h + 24),
            total_kmers: u64::from_le_bytes(buf[h + 28..h + 36].try_into().unwrap()),
        }
    } else {
        let hi = u32_at(&buf, h + 33);
        Header {
            k: u32_at(&buf, h),
            counter_size: u32_at(&buf, h + 8),
            lut_prefix_len: u32_at(&buf, h + 12),
            signature_len: 0,
            min_count: u32_at(&buf, h + 16),
            // counters are at most 4 bytes, so a max count above u32 only means "no maximum"
            max_count: if hi > 0 {
                u32::MAX
            } else {
                u32_at(&buf, h + 20)
            },
            total_kmers: u64::from_le_bytes(buf[h + 24..h + 32].try_into().unwrap()),
        }
    };
    ensure!(
        header.k <= 32,
        "{}: k = {} > 32 does not fit a u64",
        path.display(),
        header.k
    );
    ensure!(
        header.counter_size <= 4,
        "{}: counters over 4 bytes",
        path.display()
    );
    ensure!(
        (header.k - header.lut_prefix_len).is_multiple_of(4),
        "{}: suffix is not whole bytes",
        path.display()
    );
    let lut_len = if version == 0x200 {
        // size without markers and the 4-byte header offset, as in ReadParamsFrom_prefix_file_buf
        let size = n - 8 - 4;
        let sig_map = ((1usize << (2 * header.signature_len)) + 1) * 4;
        size.checked_sub(sig_map + header_offset + 8)
            .with_context(|| format!("{}: inconsistent sizes", path.display()))?
            / 8
    } else {
        1usize << (2 * header.lut_prefix_len)
    };
    ensure!(
        4 + 8 * lut_len <= h,
        "{}: LUT runs into the header",
        path.display()
    );
    let mut lut: Vec<u64> = buf[4..4 + 8 * lut_len]
        .as_chunks::<8>()
        .0
        .iter()
        .map(|c| u64::from_le_bytes(*c))
        .collect();
    lut.push(header.total_kmers); // guard: end of the last prefix
    ensure!(
        lut.windows(2).all(|w| w[0] <= w[1]),
        "{}: LUT not monotone",
        path.display()
    );
    ensure!(lut[0] == 0, "{}: LUT does not start at 0", path.display());
    Ok((header, lut))
}

pub fn header(db: &Path) -> Result<Header> {
    Ok(read_prefix(db)?.0)
}

/// Calls `f(kmer, count)` for every record, in file order (sorted within signature bins only).
/// Records outside [min_count, max_count] are skipped, as KMC's listing does.
pub fn for_each(db: &Path, mut f: impl FnMut(u64, u32)) -> Result<Header> {
    let (h, lut) = read_prefix(db)?;
    let path = with_ext(db, ".kmc_suf");
    let mut file = File::open(&path).with_context(|| format!("opening {}", path.display()))?;
    let len = file.seek(SeekFrom::End(0))?;
    let suf_bytes = ((h.k - h.lut_prefix_len) / 4) as usize;
    let rec = suf_bytes + h.counter_size as usize;
    ensure!(
        len == 8 + h.total_kmers * rec as u64,
        "{}: size does not match the header",
        path.display()
    );
    file.seek(SeekFrom::Start(0))?;
    let mut r = BufReader::with_capacity(1 << 20, file);
    let mut marker = [0u8; 4];
    r.read_exact(&mut marker)?;
    ensure!(
        &marker == b"KMCS",
        "{}: not a KMC suffix file",
        path.display()
    );

    let prefix_mask = (1u64 << (2 * h.lut_prefix_len)) - 1;
    let suf_bits = 2 * (h.k - h.lut_prefix_len);
    let mut buf = vec![0u8; rec];
    for p in 0..lut.len() - 1 {
        let prefix = (p as u64) & prefix_mask; // the LUT holds one block of 4^prefix_len per bin
        for _ in lut[p]..lut[p + 1] {
            r.read_exact(&mut buf)?;
            let suffix = buf[..suf_bytes]
                .iter()
                .fold(0u64, |v, &b| (v << 8) | b as u64);
            let count = buf[suf_bytes..]
                .iter()
                .rev()
                .fold(0u32, |c, &b| (c << 8) | b as u32);
            let count = if h.counter_size == 0 { 1 } else { count };
            if h.counter_size != 0 && (count < h.min_count || count > h.max_count) {
                continue;
            }
            f((prefix << suf_bits) | suffix, count);
        }
    }
    Ok(h)
}

/// All k-mers of a database, sorted.
pub fn read_sorted(db: &Path) -> Result<Vec<u64>> {
    let mut v = Vec::with_capacity(header(db)?.total_kmers as usize);
    for_each(db, |k, _| v.push(k))?;
    v.sort_unstable();
    Ok(v)
}

/// The k-mers of `db` that are in `wanted` (sorted), sorted.
pub fn read_present(db: &Path, wanted: &[u64]) -> Result<Vec<u64>> {
    let mut v = Vec::new();
    for_each(db, |k, _| {
        if wanted.binary_search(&k).is_ok() {
            v.push(k)
        }
    })?;
    v.sort_unstable();
    Ok(v)
}

/// One pass over `db`: its k-mers that are in `wanted` (sorted), and for the hash sample
/// (`mix(kmer) < threshold`) the number of sampled k-mers and how many of them are in `sketch`
/// (sorted). Their ratio estimates the share of the database's k-mers found in the sketched set.
pub fn read_present_sketch(
    db: &Path,
    wanted: &[u64],
    threshold: u64,
    sketch: &[u64],
) -> Result<(Vec<u64>, u64, u64)> {
    let (mut v, mut sampled, mut known) = (Vec::new(), 0u64, 0u64);
    for_each(db, |k, _| {
        if wanted.binary_search(&k).is_ok() {
            v.push(k);
        }
        if crate::kmer::mix(k) < threshold {
            sampled += 1;
            known += sketch.binary_search(&k).is_ok() as u64;
        }
    })?;
    v.sort_unstable();
    Ok((v, sampled, known))
}

/// Sequential reader of a globally sorted database (KMC 1 layout, as `kmc_tools` writes), in
/// steps of k-mer space: `read_until(end, out)` appends every k-mer whose first `bases` bases,
/// as a number, are below `end`. Many databases can be read in step, a batch at a time. The LUT
/// is streamed from the prefix file alongside the records, so memory per reader is two buffers
/// (a 4^11 LUT alone is 33.5 MB).
pub struct SortedReader {
    header: Header,
    lut: BufReader<File>,
    lut_len: usize,
    reader: BufReader<File>,
    /// next prefix (LUT index) to read, and the record range [start, end) of that prefix
    prefix: usize,
    start: u64,
    end: u64,
    suf_bytes: usize,
    rec: Vec<u8>,
}

impl SortedReader {
    pub fn open(db: &Path) -> Result<Self> {
        let path = with_ext(db, ".kmc_pre");
        let mut pre = File::open(&path).with_context(|| format!("opening {}", path.display()))?;
        let n = pre.seek(SeekFrom::End(0))? as usize;
        ensure!(n >= 16, "{}: not a KMC prefix file", path.display());
        // header: from the tail of the file (at most 255 + 8 bytes)
        let tail_len = n.min(263);
        let mut tail = vec![0u8; tail_len];
        pre.seek(SeekFrom::Start((n - tail_len) as u64))?;
        pre.read_exact(&mut tail)?;
        let tn = tail.len();
        ensure!(
            &tail[tn - 4..] == b"KMCP",
            "{}: not a KMC prefix file",
            path.display()
        );
        let version = u32_at(&tail, tn - 12);
        ensure!(
            version == 0,
            "{}: sequential reading needs the KMC 1 layout (kmc_tools output); this is {version:#x}",
            path.display()
        );
        let h = tn - 8 - tail[tn - 8] as usize;
        ensure!(
            u32_at(&tail, h + 4) == 0,
            "{}: Quake-mode counters are not supported",
            path.display()
        );
        let header = Header {
            k: u32_at(&tail, h),
            counter_size: u32_at(&tail, h + 8),
            lut_prefix_len: u32_at(&tail, h + 12),
            signature_len: 0,
            min_count: u32_at(&tail, h + 16),
            max_count: if u32_at(&tail, h + 33) > 0 {
                u32::MAX
            } else {
                u32_at(&tail, h + 20)
            },
            total_kmers: u64::from_le_bytes(tail[h + 24..h + 32].try_into().unwrap()),
        };
        ensure!(
            header.k <= 32 && header.counter_size <= 4,
            "{}: unsupported k or counter size",
            path.display()
        );
        ensure!(
            (header.k - header.lut_prefix_len).is_multiple_of(4),
            "{}: suffix is not whole bytes",
            path.display()
        );
        let lut_len = 1usize << (2 * header.lut_prefix_len);
        ensure!(
            4 + 8 * lut_len <= n - tail_len + h,
            "{}: LUT runs into the header",
            path.display()
        );
        pre.seek(SeekFrom::Start(4))?;
        let mut lut = BufReader::with_capacity(1 << 16, pre);
        let mut first = [0u8; 8];
        lut.read_exact(&mut first)?;
        ensure!(
            u64::from_le_bytes(first) == 0,
            "{}: LUT does not start at 0",
            path.display()
        );

        let suf = with_ext(db, ".kmc_suf");
        let mut file = File::open(&suf).with_context(|| format!("opening {}", suf.display()))?;
        let suf_bytes = ((header.k - header.lut_prefix_len) / 4) as usize;
        let rec = suf_bytes + header.counter_size as usize;
        ensure!(
            file.seek(SeekFrom::End(0))? == 8 + header.total_kmers * rec as u64,
            "{}: size does not match the header",
            suf.display()
        );
        file.seek(SeekFrom::Start(4))?;
        let mut r = SortedReader {
            header,
            lut,
            lut_len,
            reader: BufReader::with_capacity(1 << 16, file),
            prefix: 0,
            start: 0,
            end: 0,
            suf_bytes,
            rec: vec![0; rec],
        };
        r.end = r.next_lut()?;
        Ok(r)
    }

    /// LUT entry after the current prefix (the total after the last one).
    fn next_lut(&mut self) -> Result<u64> {
        if self.prefix + 1 >= self.lut_len {
            return Ok(self.header.total_kmers);
        }
        let mut b = [0u8; 8];
        self.lut.read_exact(&mut b)?;
        let v = u64::from_le_bytes(b);
        ensure!(
            v >= self.start && v <= self.header.total_kmers,
            "LUT not monotone"
        );
        Ok(v)
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Appends (sorted) the k-mers whose first `bases` bases are below `end`; `bases` must not
    /// exceed the database's prefix length.
    pub fn read_until(&mut self, bases: u32, end: u64, out: &mut Vec<u64>) -> Result<()> {
        let (lut_prefix_len, k) = (self.header.lut_prefix_len, self.header.k);
        ensure!(
            bases <= lut_prefix_len,
            "partition of {bases} bases is finer than the LUT ({lut_prefix_len})"
        );
        let stop = ((end << (2 * (lut_prefix_len - bases))) as usize).min(self.lut_len);
        let suf_bits = 2 * (k - lut_prefix_len);
        let (min, max, counter) = (
            self.header.min_count,
            self.header.max_count,
            self.header.counter_size,
        );
        while self.prefix < stop {
            let p = self.prefix as u64;
            for _ in self.start..self.end {
                self.reader.read_exact(&mut self.rec)?;
                let count = self.rec[self.suf_bytes..]
                    .iter()
                    .rev()
                    .fold(0u32, |c, &b| (c << 8) | b as u32);
                if counter != 0 && (count < min || count > max) {
                    continue;
                }
                let suffix = self.rec[..self.suf_bytes]
                    .iter()
                    .fold(0u64, |v, &b| (v << 8) | b as u64);
                out.push((p << suf_bits) | suffix);
            }
            self.start = self.end;
            self.prefix += 1;
            if self.prefix < self.lut_len {
                self.end = self.next_lut()?;
            }
        }
        Ok(())
    }
}

/// Writes a presence database in the KMC 1 layout that `kmc_tools` writes (counter size 1, all
/// counts 1): "KMCP" | LUT of 4^prefix_len u64 | 64-byte header, zero-padded | header length
/// (u32) | "KMCP"; and "KMCS" | records | "KMCS". `kmers`: sorted, unique, k <= 32.
pub fn write_kmc1(db: &Path, kmers: &[u64], k: u32, lut_prefix_len: u32) -> Result<()> {
    ensure!(
        k <= 32 && lut_prefix_len < k && (k - lut_prefix_len).is_multiple_of(4),
        "k - prefix length must be a multiple of 4"
    );
    ensure!(
        kmers.windows(2).all(|w| w[0] < w[1]),
        "k-mers must be sorted and unique"
    );
    let suf_bases = k - lut_prefix_len;
    let suf_bytes = (suf_bases / 4) as usize;
    let n_pre = 1usize << (2 * lut_prefix_len);
    let mut lut = vec![0u64; n_pre];
    for &x in kmers {
        let p = (x >> (2 * suf_bases)) as usize;
        ensure!(p < n_pre, "k-mer longer than k");
        if p + 1 < n_pre {
            lut[p + 1] += 1;
        }
    }
    for i in 1..n_pre {
        lut[i] += lut[i - 1];
    }
    let mut pre = BufWriter::new(File::create(with_ext(db, ".kmc_pre"))?);
    pre.write_all(b"KMCP")?;
    for v in &lut {
        pre.write_all(&v.to_le_bytes())?;
    }
    let mut header = Vec::with_capacity(64);
    for v in [k, 0, 1, lut_prefix_len, 1, u32::MAX] {
        header.extend_from_slice(&v.to_le_bytes()); // k, mode, counter size, prefix, min, max (low)
    }
    header.extend_from_slice(&(kmers.len() as u64).to_le_bytes());
    header.resize(64, 0); // strand flag, max count (high), padding
    pre.write_all(&header)?;
    pre.write_all(&64u32.to_le_bytes())?;
    pre.write_all(b"KMCP")?;
    pre.flush()?;

    let mut suf = BufWriter::new(File::create(with_ext(db, ".kmc_suf"))?);
    suf.write_all(b"KMCS")?;
    for &x in kmers {
        let s = x & ((1u64 << (2 * suf_bases)) - 1);
        suf.write_all(&s.to_be_bytes()[8 - suf_bytes..])?;
        suf.write_all(&[1])?;
    }
    suf.write_all(b"KMCS")?;
    suf.flush()?;
    Ok(())
}
