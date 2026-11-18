//! Parity with the Python prototype: the CLI must reproduce its outputs k-mer for k-mer.
//!
//! Needs the prototype's data (kmers/groups, kmers/tmp/present, kmers/samples.tsv, its outputs), so
//! the tests are ignored by default. Run with:
//!     CLADEPRINT_KMERS=$PWD/kmers cargo test --release -- --ignored
//! VNII (the largest target, ~7 min) only with CLADEPRINT_PARITY_VNII=1.

use cladeprint_core::{kmer::revcomp_str, npy};
use std::path::{Path, PathBuf};
use std::process::Command;

fn kmers_dir() -> PathBuf {
    PathBuf::from(
        std::env::var("CLADEPRINT_KMERS")
            .expect("set CLADEPRINT_KMERS to the prototype's kmers/ directory"),
    )
}

fn tmp(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("cladeprint-parity-{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d.join(name)
}

fn run(args: &[&str]) {
    let st = Command::new(env!("CARGO_BIN_EXE_cladeprint"))
        .args(args)
        .status()
        .unwrap();
    assert!(st.success(), "cladeprint {args:?} failed");
}

fn p(x: &Path) -> &str {
    x.to_str().unwrap()
}

fn filters(target: &str) {
    let g = kmers_dir().join("groups");
    let out = g.join(format!("{target}.out.u64.npy"));
    for (cmd, input, want) in [("near-miss", "strict", "nm2"), ("edit-miss", "nm2", "ed2")] {
        let got = tmp(&format!("{target}.{want}.npy"));
        let input = g.join(format!("{target}.{input}.u64.npy"));
        run(&[
            cmd,
            "--candidates",
            p(&input),
            "--outgroup",
            p(&out),
            "-o",
            p(&got),
        ]);
        let want_v = npy::read_u64(&g.join(format!("{target}.{want}.u64.npy"))).unwrap();
        assert_eq!(npy::read_u64(&got).unwrap(), want_v, "{cmd} {target}");
    }
}

#[test]
#[ignore]
fn filters_small_targets() {
    for t in ["VNIc", "VNIa", "VNIb", "VNBII"] {
        filters(t);
    }
}

#[test]
#[ignore]
fn filters_vnii() {
    if std::env::var("CLADEPRINT_PARITY_VNII").is_ok() {
        filters("VNII");
    }
}

fn unitig_set(fa: &Path) -> Vec<String> {
    let mut v: Vec<String> = std::fs::read_to_string(fa)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('>'))
        .map(|s| s.to_string().min(revcomp_str(s)))
        .collect();
    v.sort();
    v
}

#[test]
#[ignore]
fn unitigs_match() {
    let g = kmers_dir().join("groups");
    for t in ["VNII", "VNI", "VNBI", "VNBII", "VNIa", "VNIb", "VNIc"] {
        let got = tmp(&format!("{t}.unitigs.fa"));
        run(&[
            "unitigs",
            p(&g.join(format!("{t}.ed2.u64.npy"))),
            "-o",
            p(&got),
        ]);
        // isolated cycles could be cut at a different k-mer; there are none in these sets
        assert_eq!(
            unitig_set(&got),
            unitig_set(&g.join(format!("{t}.ed2.unitigs.fa"))),
            "{t}"
        );
    }
}

fn table(path: &Path) -> Vec<Vec<String>> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|l| l.split('\t').map(String::from).collect())
        .collect()
}

#[test]
#[ignore]
fn scores_match() {
    let k = kmers_dir();
    let prefix = tmp("rs_");
    run(&[
        "score",
        "--groups",
        p(&k.join("groups")),
        "--samples",
        p(&k.join("samples.tsv")),
        "--present",
        p(&k.join("tmp/present")),
        "--prefix",
        p(&prefix),
    ]);
    for f in ["scores.tsv", "loo_sensitivity.tsv"] {
        let got = table(Path::new(&format!("{}{f}", p(&prefix))));
        let want = table(&k.join(f));
        assert_eq!(got.len(), want.len(), "{f} rows");
        for (r, (a, b)) in got.iter().zip(&want).enumerate() {
            assert_eq!(a.len(), b.len(), "{f} row {r}");
            for (x, y) in a.iter().zip(b) {
                match (x.parse::<f64>(), y.parse::<f64>()) {
                    (Ok(u), Ok(v)) => assert!((u - v).abs() < 1e-9, "{f} row {r}: {x} vs {y}"),
                    _ => assert_eq!(x, y, "{f} row {r}"),
                }
            }
        }
    }
}

#[test]
#[ignore]
fn kmc_dump_matches_python_npy() {
    let g = kmers_dir().join("groups");
    for db in ["VNIc.strict", "VNIa.strict", "VNBII.strict", "VNIc.out"] {
        let got = tmp(&format!("{db}.dump.npy"));
        run(&["kmc-dump", p(&g.join(db)), "-o", p(&got)]);
        let want = npy::read_u64(&g.join(format!("{db}.u64.npy"))).unwrap();
        assert_eq!(npy::read_u64(&got).unwrap(), want, "{db}");
    }
}

/// Scoring straight from the 123 per-run KMC databases (~9 min) must give the same tables as the
/// prototype's kmc_tools intersect + dump route. Only with CLADEPRINT_PARITY_KMCDB=1.
#[test]
#[ignore]
fn scores_from_kmc_databases() {
    if std::env::var("CLADEPRINT_PARITY_KMCDB").is_err() {
        return;
    }
    let k = kmers_dir();
    let prefix = tmp("kmcdb_");
    run(&[
        "score",
        "--groups",
        p(&k.join("groups")),
        "--samples",
        p(&k.join("samples.tsv")),
        "--kmc-db",
        p(&k.join("db")),
        "--prefix",
        p(&prefix),
    ]);
    for f in ["scores.tsv", "loo_sensitivity.tsv"] {
        let got = std::fs::read_to_string(format!("{}{f}", p(&prefix))).unwrap();
        assert_eq!(got, std::fs::read_to_string(k.join(f)).unwrap(), "{f}");
    }
}

/// Group counting over the 108 training databases (~10 min) must reproduce group_kmers.py and the
/// loo_extra.sh candidates: the summary table and every strict and extra k-mer set.
/// Only with CLADEPRINT_PARITY_GROUPS=1.
#[test]
#[ignore]
fn groups_match() {
    if std::env::var("CLADEPRINT_PARITY_GROUPS").is_err() {
        return;
    }
    let k = kmers_dir();
    let out = tmp("groups");
    run(&[
        "groups",
        "--samples",
        p(&k.join("samples.tsv")),
        "--kmc-db",
        p(&k.join("db")),
        "-o",
        p(&out),
    ]);
    let got = table(&out.join("groups.tsv"));
    let want = table(&k.join("group_kmers_fin0.9.tsv"));
    assert_eq!(got.len(), want.len());
    for (a, b) in got.iter().zip(&want) {
        assert_eq!(a[..7], b[..7], "summary row {}", b[0]);
    }
    let g = k.join("groups");
    for row in &want[1..] {
        let t = &row[0];
        for w in ["strict", "x.strict"] {
            let mine = npy::read_u64(&out.join(format!("{t}.{w}.u64.npy"))).unwrap();
            let npy_path = g.join(format!("{t}.{w}.u64.npy"));
            let theirs = if npy_path.exists() {
                npy::read_u64(&npy_path).unwrap()
            } else {
                cladeprint_core::kmc::read_sorted(&g.join(format!("{t}.{w}"))).unwrap()
            };
            assert_eq!(mine, theirs, "{t}.{w}");
        }
    }
}
