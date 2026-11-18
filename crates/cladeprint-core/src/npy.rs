//! Read and write 1-D little-endian uint64 `.npy` arrays, the interchange format of the Python
//! prototype (`groups/T.*.u64.npy`).

use anyhow::{Context, Result, bail};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

pub fn read_u64(path: &Path) -> Result<Vec<u64>> {
    let mut f =
        BufReader::new(File::open(path).with_context(|| format!("opening {}", path.display()))?);
    let mut magic = [0u8; 8];
    f.read_exact(&mut magic)?;
    if &magic[..6] != b"\x93NUMPY" {
        bail!("{}: not a .npy file", path.display());
    }
    let hlen = match magic[6] {
        1 => {
            let mut b = [0u8; 2];
            f.read_exact(&mut b)?;
            u16::from_le_bytes(b) as usize
        }
        2 | 3 => {
            let mut b = [0u8; 4];
            f.read_exact(&mut b)?;
            u32::from_le_bytes(b) as usize
        }
        v => bail!("{}: unsupported .npy version {v}", path.display()),
    };
    let mut header = vec![0u8; hlen];
    f.read_exact(&mut header)?;
    let header = String::from_utf8_lossy(&header);
    if !header.contains("'descr': '<u8'") || !header.contains("'fortran_order': False") {
        bail!(
            "{}: expected a little-endian uint64 array, header {header}",
            path.display()
        );
    }
    let shape = header
        .split("'shape': (")
        .nth(1)
        .and_then(|s| s.split(')').next())
        .with_context(|| format!("{}: no shape in header", path.display()))?;
    let dims: Vec<&str> = shape
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .collect();
    if dims.len() != 1 {
        bail!("{}: expected a 1-D array, shape ({shape})", path.display());
    }
    let n: usize = dims[0].parse()?;
    let mut bytes = vec![0u8; n * 8];
    f.read_exact(&mut bytes)?;
    Ok(bytes
        .as_chunks::<8>()
        .0
        .iter()
        .map(|c| u64::from_le_bytes(*c))
        .collect())
}

/// Writes a version 1.0 file with the same header layout as `numpy.save`.
pub fn write_u64(path: &Path, v: &[u64]) -> Result<()> {
    let mut f =
        BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    let mut header = format!(
        "{{'descr': '<u8', 'fortran_order': False, 'shape': ({},), }}",
        v.len()
    );
    // magic (6) + version (2) + length (2) + header, padded with spaces to a multiple of 64, then '\n'
    let total = 10 + header.len() + 1;
    header.push_str(&" ".repeat((64 - total % 64) % 64));
    header.push('\n');
    f.write_all(b"\x93NUMPY\x01\x00")?;
    f.write_all(&(header.len() as u16).to_le_bytes())?;
    f.write_all(header.as_bytes())?;
    for x in v {
        f.write_all(&x.to_le_bytes())?;
    }
    f.flush()?;
    Ok(())
}

/// Streams a 1-D uint64 array to disk; the length goes into the header at `finish`. The header
/// is a fixed 128 bytes (numpy reads any header length), so it can be rewritten in place.
pub struct Writer {
    f: BufWriter<File>,
    n: u64,
}

const STREAM_HEADER: usize = 128;

fn stream_header(n: u64) -> Vec<u8> {
    let mut h = format!("{{'descr': '<u8', 'fortran_order': False, 'shape': ({n},), }}");
    h.push_str(&" ".repeat(STREAM_HEADER - 10 - h.len() - 1));
    h.push('\n');
    let mut out = b"\x93NUMPY\x01\x00".to_vec();
    out.extend_from_slice(&(h.len() as u16).to_le_bytes());
    out.extend_from_slice(h.as_bytes());
    out
}

impl Writer {
    pub fn create(path: &Path) -> Result<Self> {
        let mut f = BufWriter::new(
            File::create(path).with_context(|| format!("creating {}", path.display()))?,
        );
        f.write_all(&stream_header(0))?;
        Ok(Writer { f, n: 0 })
    }

    pub fn write(&mut self, v: &[u64]) -> Result<()> {
        for x in v {
            self.f.write_all(&x.to_le_bytes())?;
        }
        self.n += v.len() as u64;
        Ok(())
    }

    pub fn finish(self) -> Result<u64> {
        use std::io::{Seek, SeekFrom};
        let mut f = self.f.into_inner().map_err(|e| e.into_error())?;
        f.seek(SeekFrom::Start(0))?;
        f.write_all(&stream_header(self.n))?;
        Ok(self.n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn streamed_roundtrip() {
        let p =
            std::env::temp_dir().join(format!("cladeprint-npy-stream-{}.npy", std::process::id()));
        let mut w = Writer::create(&p).unwrap();
        w.write(&[1, 2, 3]).unwrap();
        w.write(&[]).unwrap();
        w.write(&[u64::MAX]).unwrap();
        assert_eq!(w.finish().unwrap(), 4);
        assert_eq!(read_u64(&p).unwrap(), vec![1, 2, 3, u64::MAX]);
        std::fs::remove_file(&p).unwrap();
    }

    #[test]
    fn roundtrip() {
        let dir = std::env::temp_dir().join(format!("cladeprint-npy-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        for n in [0usize, 1, 5, 1000] {
            let v: Vec<u64> = (0..n as u64)
                .map(|i| i.wrapping_mul(0x9E37_79B9_7F4A_7C15))
                .collect();
            let p = dir.join(format!("{n}.npy"));
            write_u64(&p, &v).unwrap();
            assert_eq!(
                std::fs::metadata(&p).unwrap().len() as usize % 64,
                (n * 8) % 64
            );
            assert_eq!(read_u64(&p).unwrap(), v);
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
