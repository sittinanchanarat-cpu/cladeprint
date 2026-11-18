//! KMC reader against small databases written by KMC 3.2.4, with `kmc_tools transform dump` as
//! the expected output (tests/data). v2: `kmc -k31 -ci1 -cs65535` (format 0x200, signature bins);
//! v1: `kmc_tools transform v2 reduce` (KMC 1 layout); v1_ge3: `kmc_tools transform v2 -ci3
//! set_counts 1`, the recipe of the presence databases. Input: overlapping reads of a random 3 kb
//! sequence, part of it at higher depth, one read with an N.

use cladeprint_core::{kmc, kmer};
use std::collections::BTreeMap;
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn dump(name: &str) -> BTreeMap<u64, u32> {
    std::fs::read_to_string(data(&format!("{name}.txt")))
        .unwrap()
        .lines()
        .map(|l| {
            let (k, c) = l.split_once('\t').unwrap();
            (kmer::encode(k.as_bytes()).unwrap(), c.parse().unwrap())
        })
        .collect()
}

fn read(name: &str) -> (kmc::Header, BTreeMap<u64, u32>) {
    let mut got = BTreeMap::new();
    let h = kmc::for_each(&data(name), |k, c| {
        assert!(got.insert(k, c).is_none(), "{name}: k-mer listed twice");
    })
    .unwrap();
    (h, got)
}

#[test]
fn kmc2_format_with_counts() {
    let (h, got) = read("v2");
    assert_eq!((h.k, h.signature_len > 0), (31, true));
    let want = dump("v2");
    assert_eq!(got, want);
    assert!(want.values().any(|&c| c > 1) && want.values().any(|&c| c == 1));
}

#[test]
fn kmc1_format_from_kmc_tools() {
    let (h, got) = read("v1");
    assert_eq!(h.signature_len, 0);
    assert_eq!(got, dump("v1"));
    assert_eq!(got, dump("v2")); // same k-mers and counts, other layout
}

#[test]
fn presence_database_recipe() {
    let (_, got) = read("v1_ge3");
    let want = dump("v1_ge3");
    assert_eq!(got, want);
    assert!(got.values().all(|&c| c == 1));
    let v2 = dump("v2");
    let ge3: Vec<u64> = v2
        .iter()
        .filter(|&(_, &c)| c >= 3)
        .map(|(&k, _)| k)
        .collect();
    assert_eq!(got.keys().copied().collect::<Vec<_>>(), ge3);
}

#[test]
fn sorted_and_present_helpers() {
    let all = kmc::read_sorted(&data("v2")).unwrap();
    assert!(all.windows(2).all(|w| w[0] < w[1]));
    assert_eq!(all, dump("v2").keys().copied().collect::<Vec<_>>());
    // canonical k-mers: every k-mer is <= its reverse complement
    assert!(all.iter().all(|&k| k <= kmer::revcomp(k)));
    let wanted: Vec<u64> = all
        .iter()
        .step_by(3)
        .copied()
        .chain([0, u64::MAX >> 2])
        .collect();
    let mut wanted_sorted = wanted.clone();
    wanted_sorted.sort_unstable();
    let present = kmc::read_present(&data("v1"), &wanted_sorted).unwrap();
    assert_eq!(present, all.iter().step_by(3).copied().collect::<Vec<_>>());
}

/// `groups::run` against a brute-force count on six small isolates (tests/data/groups, built with
/// kmc + `kmc_tools transform -ci2 set_counts 1` from reads of a shared 2 kb genome with private
/// SNPs). Planted 150 bp segments: one in all of group A (1-3); one in all of group B (4-6); one in
/// B isolates 4-5 only (in = need - 1: extra); one in all of A plus isolate 4 (out = 1).
#[test]
fn group_sketch_matches_brute_force() {
    use cladeprint_core::groups::{Target, run_sketch};
    use cladeprint_core::kmer::{mix, sketch_threshold};
    use std::collections::HashMap;
    let dbs: Vec<PathBuf> = (1..=6).map(|i| data(&format!("groups/iso{i}"))).collect();
    let targets = vec![Target::new("A", vec![0, 1, 2], 0.9)];
    let thr = sketch_threshold(4);
    let mut sketch = Vec::new();
    let (_, sk) = run_sketch(&dbs, &targets, false, 7, Some(thr), |b| {
        sketch.extend(b.sketch);
        Ok(())
    })
    .unwrap();
    let mut holders: HashMap<u64, Vec<usize>> = HashMap::new();
    for (i, d) in dbs.iter().enumerate() {
        for k in kmc::read_sorted(d).unwrap() {
            holders.entry(k).or_default().push(i);
        }
    }
    let mut want: Vec<u64> = holders.keys().copied().filter(|&k| mix(k) < thr).collect();
    want.sort_unstable();
    assert_eq!(sketch, want, "sampled k-mers, in order");
    for i in 0..dbs.len() {
        let mine: Vec<&Vec<usize>> = holders
            .iter()
            .filter(|(k, h)| mix(**k) < thr && h.contains(&i))
            .map(|(_, h)| h)
            .collect();
        assert_eq!(
            sk.total[i],
            mine.len() as u64,
            "sampled k-mers of isolate {i}"
        );
        assert_eq!(
            sk.unique[i],
            mine.iter().filter(|h| h.len() == 1).count() as u64,
            "unique k-mers of isolate {i}"
        );
    }
    assert!(
        sk.unique.iter().all(|&u| u > 0),
        "private SNPs give every isolate some unique k-mers"
    );
}

#[test]
fn group_counting_matches_brute_force() {
    use cladeprint_core::groups::{Target, cluster_weights};
    group_counting_case(
        &[
            Target::new("A", vec![0, 1, 2], 0.9),
            Target::new("B", vec![3, 4, 5], 0.9),
        ],
        |strict, extra, stats| {
            assert!(
                extra[1].len() >= 100 && extra[0].is_empty(),
                "extras A {} B {}",
                extra[0].len(),
                extra[1].len()
            );
            assert!(
                stats[0].out_le1 - stats[0].out0 >= 100,
                "A: k-mers shared with one B isolate"
            );
            assert_eq!(
                stats[0].strict, stats[0].out0,
                "no tolerance: strict = out0"
            );
            assert!(strict[0].len() >= 100);
        },
    );
    // clade weights (isolates 1-2 one cluster in A, 5-6 in B) and a tolerance of one out-group
    // isolate for A: the segment in all of A plus isolate 4 becomes an A marker
    let wa = cluster_weights(&[Some("x"), Some("x"), None]);
    let wb = cluster_weights(&[None, Some("y"), Some("y")]);
    group_counting_case(
        &[
            Target::weighted("A", vec![0, 1, 2], wa, 0.9, 1),
            Target::weighted("B", vec![3, 4, 5], wb, 0.9, 0),
        ],
        |_, extra, stats| {
            assert!(
                stats[0].strict - stats[0].out0 >= 100,
                "A: the k-mers shared with isolate 4 are tolerated"
            );
            // B's 4-5 segment: in = 1 + 0.5 = 1.5 of 2, below need 1.8, above extra_lo 0.9
            assert!(extra[1].len() >= 100, "extras B {}", extra[1].len());
        },
    );
}

/// Ignored isolates: A = isolates 1-2 with isolate 3 unlabelled at A's level (ignored). The A
/// segment (in 1-3) stays a marker because 3 counts neither in nor out; without the ignore it
/// would be held by an out-group isolate. B's out-group no longer includes 3 either.
#[test]
fn group_counting_ignores_unlabelled() {
    use cladeprint_core::groups::Target;
    let plain = || Target::new("A", vec![0, 1], 0.9);
    // without the ignore, isolate 3 is out-group and holds the A segment: few strict k-mers
    group_counting_case(&[plain()], |_, _, stats| {
        assert!(stats[0].strict < 50, "A strict {}", stats[0].strict)
    });
    group_counting_case(
        &[
            plain().with_ignore(vec![2]),
            Target::new("B", vec![3, 4, 5], 0.9).with_ignore(vec![2]),
        ],
        |strict, _, stats| {
            assert!(
                stats[0].strict >= 100,
                "A segment kept with isolate 3 ignored: {}",
                stats[0].strict
            );
            assert!(strict[1].len() >= 100);
        },
    );
}

/// `groups::run` for `targets` against a brute-force count, at several batch sizes; `check` gets
/// strict and extra sets and the stats for planted-segment assertions.
fn group_counting_case(
    targets: &[cladeprint_core::groups::Target],
    check: impl Fn(&[Vec<u64>], &[Vec<u64>], &[cladeprint_core::groups::Stats]),
) {
    use cladeprint_core::groups::{Stats, run};
    use std::collections::HashMap;
    let dbs: Vec<PathBuf> = (1..=6).map(|i| data(&format!("groups/iso{i}"))).collect();
    let nt = targets.len();
    let mut holders: HashMap<u64, Vec<usize>> = HashMap::new();
    for (i, d) in dbs.iter().enumerate() {
        for k in kmc::read_sorted(d).unwrap() {
            holders.entry(k).or_default().push(i);
        }
    }
    for batches in [1, 7, 64] {
        let (mut strict, mut extra, mut out) =
            (vec![vec![]; nt], vec![vec![]; nt], vec![vec![]; nt]);
        let (mut counts, mut xcounts): (Vec<Vec<f64>>, Vec<Vec<f64>>) =
            (vec![vec![]; nt], vec![vec![]; nt]);
        let stats = run(&dbs, targets, true, batches, |b| {
            for t in 0..nt {
                strict[t].extend(&b.strict[t]);
                counts[t].extend(&b.strict_count[t]);
                extra[t].extend(&b.extra[t]);
                xcounts[t].extend(&b.extra_count[t]);
                out[t].extend(&b.out[t]);
            }
            Ok(())
        })
        .unwrap();
        for (t, tg) in targets.iter().enumerate() {
            let in_w = |h: &[usize]| -> f64 {
                tg.members
                    .iter()
                    .zip(&tg.weights)
                    .filter(|(m, _)| h.contains(m))
                    .map(|(_, w)| w)
                    .sum()
            };
            let (mut ws, mut wx, mut wo) = (vec![], vec![], vec![]);
            let mut want = Stats::default();
            for (&k, h) in &holders {
                let i = in_w(h);
                let o = h
                    .iter()
                    .filter(|x| !tg.members.contains(x) && !tg.ignore.contains(x))
                    .count() as u32;
                if i >= tg.need {
                    want.candidates += 1;
                    want.out0 += (o == 0) as u64;
                    want.out_le1 += (o <= 1) as u64;
                    if o <= tg.max_out {
                        ws.push(k);
                    }
                } else if i >= tg.extra_lo && o <= tg.max_out {
                    wx.push(k);
                }
                if o > tg.max_out {
                    wo.push(k);
                }
            }
            ws.sort_unstable();
            wx.sort_unstable();
            wo.sort_unstable();
            (want.strict, want.extra, want.out_kmers) =
                (ws.len() as u64, wx.len() as u64, wo.len() as u64);
            assert_eq!(stats[t], want, "target {} batches {batches}", tg.name);
            assert_eq!(
                counts[t].len(),
                strict[t].len(),
                "one in-group weight per strict k-mer"
            );
            assert_eq!(
                xcounts[t].len(),
                extra[t].len(),
                "one in-group weight per extra"
            );
            for (k, c) in strict[t]
                .iter()
                .zip(&counts[t])
                .chain(extra[t].iter().zip(&xcounts[t]))
            {
                assert!(
                    (c - in_w(&holders[k])).abs() < 1e-9,
                    "in-group weight of {k}"
                );
            }
            assert_eq!(
                (&strict[t], &extra[t], &out[t]),
                (&ws, &wx, &wo),
                "target {}",
                tg.name
            );
            // the planted segments: ~120 strict k-mers per group (when the rule lets them through)
            if tg.members.len() == 3 || !tg.ignore.is_empty() {
                assert!(ws.len() >= 100, "{} strict {}", tg.name, ws.len());
            }
        }
        check(&strict, &extra, &stats);
    }
}

/// `clusters::run` against a brute-force genotyping of the six fixture isolates (shared genome,
/// private SNPs): sites, per-sample stats and every pair's exact distance, keeping all alleles and
/// ignoring run-private ones.
#[test]
fn clusters_match_brute_force() {
    use cladeprint_core::clusters::{Options, run, split};
    use cladeprint_core::kmer::{mix, sketch_threshold};
    use std::collections::{BTreeMap, HashMap};
    let dbs: Vec<PathBuf> = (1..=6).map(|i| data(&format!("groups/iso{i}"))).collect();
    let scale = 4;
    // site key -> per-sample allele mask
    let mut raw: BTreeMap<u64, Vec<u8>> = BTreeMap::new();
    for (i, d) in dbs.iter().enumerate() {
        for k in kmc::read_sorted(d).unwrap() {
            if let Some(v) = split(k) {
                raw.entry(v >> 2).or_insert_with(|| vec![0; 6])[i] |= 1 << (v & 3);
            }
        }
    }
    let thr = sketch_threshold(scale);
    let single = |m: u8| m.count_ones() == 1;
    for min_allele in [1u32, 2] {
        let mut want_pos = [0f64; 6];
        let mut geno: Vec<Vec<Option<u8>>> = Vec::new(); // per kept site
        let (mut want_sites, mut want_mixed) = ([0u64; 6], [0u64; 6]);
        for (&key, m) in &raw {
            // drop alleles held by fewer than min_allele samples
            let keep: Vec<bool> = (0..4)
                .map(|a| m.iter().filter(|&&x| x >> a & 1 == 1).count() as u32 >= min_allele)
                .collect();
            let m: Vec<u8> = m
                .iter()
                .map(|&x| {
                    (0..4)
                        .filter(|&a| keep[a] && x >> a & 1 == 1)
                        .map(|a| 1u8 << a)
                        .sum()
                })
                .collect();
            if mix(key) < thr {
                for i in 0..6 {
                    want_pos[i] += single(m[i]) as u8 as f64 * scale as f64;
                }
            }
            let alleles: std::collections::HashSet<u32> = m
                .iter()
                .filter(|&&x| single(x))
                .map(|x| x.trailing_zeros())
                .collect();
            if alleles.len() >= 2 {
                geno.push(
                    m.iter()
                        .map(|&x| single(x).then(|| x.trailing_zeros() as u8))
                        .collect(),
                );
                for i in 0..6 {
                    want_sites[i] += single(m[i]) as u64;
                    want_mixed[i] += (m[i].count_ones() > 1) as u64;
                }
            }
        }
        for batches in [0, 7] {
            let opt = Options {
                scale,
                overlap_scale: scale,
                min_positions: 0.0,
                exact_max: f64::INFINITY,
                min_allele,
                batches,
            };
            let r = run(&dbs, &opt).unwrap();
            let what = format!("min_allele {min_allele} batches {batches}");
            assert_eq!(r.sites, geno.len(), "{what}");
            for i in 0..6 {
                let s = &r.samples[i];
                assert_eq!(
                    (s.sites, s.mixed, s.positions),
                    (want_sites[i], want_mixed[i], want_pos[i]),
                    "{what}, sample {i}"
                );
            }
            assert_eq!(r.pairs.len(), 15, "every pair exact");
            let got: HashMap<(u32, u32), (u64, u64)> = r
                .pairs
                .iter()
                .map(|p| ((p.i, p.j), (p.sites, p.diffs)))
                .collect();
            for i in 0..6 {
                for j in i + 1..6 {
                    let (mut sh, mut df) = (0, 0);
                    for g in &geno {
                        if let (Some(a), Some(b)) = (g[i], g[j]) {
                            sh += 1;
                            df += (a != b) as u64;
                        }
                    }
                    assert_eq!(got[&(i as u32, j as u32)], (sh, df), "{what}, pair {i} {j}");
                }
            }
            // the fixture's SNPs are all private: sites only when private alleles count
            if min_allele == 1 {
                assert!(r.sites >= 10, "{what}: {} sites", r.sites);
            } else {
                assert_eq!(r.sites, 0, "{what}");
            }
        }
    }
}
