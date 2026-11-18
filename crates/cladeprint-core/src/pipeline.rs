//! `cladeprint run`: per-sample counting, training-set QC, group k-mers for every group at every
//! level, near-miss filters, unitigs, leave-one-out sensitivity and scores, and the marker bundle.
//!
//!
//! Output directory:
//!   db/<sample>.kmc_pre|suf + .stats.toml   per-sample presence databases (reused when present)
//!   qc.tsv        per sample: input, counting stats, whether it trained and why
//!   scores.tsv    per sample x group: share of the group's markers present (leave-one-out
//!                 markers for a training sample's own groups)
//!   loo.tsv       per training sample x own group: leave-one-out sensitivity
//!   bundle/       what `classify` needs: bundle.toml, groups.tsv, markers/<level>/<group>.kmc_*|.fa
//!   tmp/          group-counting intermediates, removed as they are used

use crate::config::{Config, Group, Input, Sample, Train, hierarchy, load_samples};
use crate::count::{self, CountStats};
use crate::groups::{self, Target};
use crate::score;
use crate::{edit, hamming, kmc, npy, unitig};
use anyhow::{Context, Result, ensure};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

pub const BUNDLE_FORMAT: u32 = 4;
/// the training k-mer set is kept as a 1/SKETCH_SCALE hash sample (bundle/sketch.u64.npy)
pub const SKETCH_SCALE: u64 = 64;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BundleGroup {
    pub id: String,
    pub level: usize,
    pub name: String,
    /// marker file stem under bundle/, e.g. "markers/lineage/Cn__VNI"
    pub file: String,
    pub n_train: usize,
    /// summed clade weight of the training members (= n_train without clade weighting)
    #[serde(default)]
    pub n_eff: f64,
    /// in-group weight a marker needed (a member count without clade weighting)
    pub need: f64,
    pub markers: usize,
    /// leave-one-out sensitivity of the training members: min, median, max
    pub loo: [f64; 3],
    /// 5th percentile of the leave-one-out sensitivity (nearest rank; the minimum below 11
    /// members): classify's bar for a high-confidence call. None in bundles built before it
    /// (classify then falls back to the minimum; census/add_loo_p05.py fills it from loo.tsv).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub loo_p05: Option<f64>,
    /// leave-one-cluster-out sensitivity (the honest figure for clonal data): min, median, max
    #[serde(default)]
    pub loco: [f64; 3],
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Bundle {
    pub format: u32,
    pub cladeprint_version: String,
    pub k: u32,
    pub levels: Vec<String>,
    pub min_peak_coverage: u32,
    pub groups: Vec<BundleGroup>,
    /// groups not built (too few training samples), as "id (n training)"
    pub skipped: Vec<String>,
    /// hash sample of every training k-mer: 1/sketch_scale of them, in `sketch_file`
    pub sketch_scale: u64,
    pub sketch_file: String,
    /// leave-one-out share of a training sample's k-mers found in no other training sample:
    /// min, median, 95th, 99th percentile, max
    pub unexplained_loo: [f64; 5],
    /// classify flags UNEXPLAINED above this share
    pub unexplained_flag: f64,
    /// median solid k-mers of good-coverage samples per top-level group and overall ("*")
    pub solid_reference: std::collections::BTreeMap<String, f64>,
    /// classify flags EXCESS_KMERS above this ratio to the reference (= training.max_kmer_excess)
    pub max_kmer_excess: f64,
    pub config: Config,
}

/// File-name-safe form of a sample name or group id.
/// Peak resident memory so far (VmHWM, Linux), for the stage log.
fn peak_rss() -> String {
    fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            let l = s.lines().find(|l| l.starts_with("VmHWM:"))?;
            l.split_whitespace().nth(1)?.parse::<f64>().ok()
        })
        .map_or("n/a".into(), |kb| format!("{:.2} GB", kb / 1048576.0))
}

pub fn sanitize(s: &str) -> String {
    s.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || "._-".contains(c) {
                c
            } else {
                '_'
            }
        })
        .collect()
}

fn db_path(cfg: &Config, s: &Sample) -> PathBuf {
    cfg.output.join("db").join(sanitize(&s.name))
}

fn with_ext(p: &Path, ext: &str) -> PathBuf {
    let mut s = p.as_os_str().to_owned();
    s.push(ext);
    PathBuf::from(s)
}

/// Presence database of every sample; an existing database with its stats file is reused.
fn count_all(cfg: &Config, samples: &[Sample]) -> Result<Vec<CountStats>> {
    fs::create_dir_all(cfg.output.join("db"))?;
    let tmp = cfg.tmp().join("count");
    let mut out = Vec::new();
    for (i, s) in samples.iter().enumerate() {
        let db = db_path(cfg, s);
        let stats_file = with_ext(&db, ".stats.toml");
        if stats_file.exists() && with_ext(&db, ".kmc_suf").exists() {
            out.push(toml::from_str(&fs::read_to_string(&stats_file)?)?);
            continue;
        }
        eprintln!("[count {}/{}] {}", i + 1, samples.len(), s.name);
        let st = match &s.input {
            Input::Reads(files) => count::count_reads(
                files,
                &db,
                &tmp,
                cfg.resources.threads,
                cfg.resources.memory_gb,
                cfg.kmers.min_peak_coverage,
                0, // training: samples short of the reference are dropped by min_kmer_excess
            )?,
            Input::Assembly(fa) => count::count_assembly(fa, &db)?,
        };
        fs::write(&stats_file, toml::to_string(&st)?)?; // written last: marks the database complete
        out.push(st);
    }
    Ok(out)
}

/// Training decision and its reason.
fn training(
    s: &Sample,
    st: &CountStats,
    excess: f64,
    cfg: &crate::config::Training,
) -> (bool, &'static str) {
    let high = excess > cfg.max_kmer_excess;
    let low = excess < cfg.min_kmer_excess;
    match (s.train, st.low_coverage, high, low) {
        (Train::No, ..) => (false, "train=no in sample sheet"),
        (Train::Yes, true, ..) => (true, "train=yes in sample sheet (low coverage)"),
        (Train::Yes, false, true, _) => (true, "train=yes in sample sheet (excess k-mers)"),
        (Train::Yes, false, _, true) => (true, "train=yes in sample sheet (too few k-mers)"),
        (Train::Yes, ..) => (true, "train=yes in sample sheet"),
        (Train::Auto, true, ..) => (false, "low coverage"),
        (Train::Auto, false, true, _) => (false, "excess k-mers: mixture, hybrid or contamination"),
        (Train::Auto, false, _, true) => (
            false,
            "too few k-mers: foreign organism, partial genome or degraded run",
        ),
        (Train::Auto, ..) => (true, "passes QC"),
    }
}

fn median_of(mut v: Vec<f64>) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    Some(if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    })
}

/// Reference solid k-mer counts: the median over good-coverage samples, per top-level group with
/// at least 3 such samples, and over all of them ("*"). Two genomes in one sample (a mixture, a
/// heterozygous hybrid) or contamination add distinct k-mers: haploid isolates sit at <= 1.03 x
/// their reference, mixtures and hybrids at >= 1.19 x.
pub fn solid_reference(
    samples: &[Sample],
    stats: &[CountStats],
) -> std::collections::BTreeMap<String, f64> {
    let good = |i: usize| !stats[i].low_coverage;
    let mut r = std::collections::BTreeMap::new();
    if let Some(m) = median_of(
        (0..samples.len())
            .filter(|&i| good(i))
            .map(|i| stats[i].solid_kmers as f64)
            .collect(),
    ) {
        r.insert("*".to_string(), m);
    }
    let mut tops: Vec<&str> = samples
        .iter()
        .filter_map(|s| s.labels.first()?.as_deref())
        .collect();
    tops.sort();
    tops.dedup();
    for top in tops {
        let v: Vec<f64> = (0..samples.len())
            .filter(|&i| {
                good(i) && samples[i].labels.first().and_then(|l| l.as_deref()) == Some(top)
            })
            .map(|i| stats[i].solid_kmers as f64)
            .collect();
        if v.len() >= 3 {
            r.insert(top.to_string(), median_of(v).unwrap());
        }
    }
    r
}

/// Solid k-mers relative to the reference of the top-level group (or "*").
pub fn kmer_excess(
    reference: &std::collections::BTreeMap<String, f64>,
    top: Option<&str>,
    solid: u64,
) -> f64 {
    let r = top
        .and_then(|t| reference.get(t))
        .or_else(|| reference.get("*"));
    r.map_or(f64::NAN, |&r| solid as f64 / r)
}

/// Markers of one group after both near-miss filters, with their in-group weights, and the
/// filtered leave-one-out extras with theirs.
struct GroupMarkers {
    markers: Vec<u64>,
    /// `stride` weights per marker: the group's in-weight then each child's
    counts: Vec<f64>,
    extras: Vec<u64>,
    extras_count: Vec<f64>,
    /// 1 + number of child constraints
    stride: usize,
    stats: groups::Stats,
    hamming_survive: usize,
    extras_raw: usize,
    /// share of the strict candidates kept by the cap (1 = all)
    kept_share: f64,
}

fn filter(cand: &[u64], out: &[u64]) -> (Vec<usize>, usize) {
    let d1 = hamming::candidate_dist(cand, out);
    let h: Vec<usize> = (0..cand.len()).filter(|&i| d1[i] == hamming::FAR).collect();
    let after_h: Vec<u64> = h.iter().map(|&i| cand[i]).collect();
    let d2 = edit::candidate_dist(&after_h, out);
    let keep = h
        .iter()
        .zip(&d2)
        .filter(|&(_, &d)| d == edit::FAR)
        .map(|(&i, _)| i)
        .collect();
    (keep, after_h.len())
}

/// Leave-one-cluster-out result of one training member: (cluster size, markers, present, score).
type LocoResult = (usize, usize, usize, f64);

/// Leave-one-cluster-out scores of one group's members that sit in multi-member clusters inside
/// the current scoring batch (`slot`: sample index -> position in `present`). Clusters of one are
/// left out: holding one out is exactly leave-one-out. Per cluster, the weight it contributes to
/// every marker and extra and the thresholds without it are built, its members scored, and the
/// arrays dropped, so memory holds one cluster's arrays per thread, not every cluster's.
fn loco_scores(
    ti: usize,
    targets: &[Target],
    m: &GroupMarkers,
    samples: &[Sample],
    train_idx: &[usize],
    present: &[Vec<u64>],
    slot: &std::collections::HashMap<usize, usize>,
) -> Result<std::collections::HashMap<usize, LocoResult>> {
    use std::collections::HashMap;
    let t = &targets[ti];
    let stride = m.stride;
    // training members of this group, by cluster name
    let mut by_cluster: HashMap<&str, Vec<usize>> = HashMap::new();
    for (j, &p) in t.members.iter().enumerate() {
        let name = samples[train_idx[p]].cluster.as_deref().unwrap_or("");
        if !name.is_empty() {
            by_cluster.entry(name).or_default().push(j);
        }
    }
    by_cluster.retain(|_, v| v.len() > 1);
    // a batch holds whole clusters (score_batches): a cluster is in it entirely or not at all
    let in_batch = |j: usize| slot.contains_key(&train_idx[t.members[j]]);
    by_cluster.retain(|_, js| in_batch(js[0]));
    for (name, js) in &by_cluster {
        ensure!(
            js.iter().all(|&j| in_batch(j)),
            "cluster {name} is split across scoring batches"
        );
    }
    // weight of training sample `p` in level `c` of this target (0 = the group itself)
    let weight_at = |c: usize, p: usize| -> f64 {
        if c == 0 {
            t.members
                .iter()
                .position(|&q| q == p)
                .map_or(0.0, |k| t.weights[k])
        } else {
            let ct = &targets[t.children[c - 1]];
            ct.members
                .iter()
                .position(|&q| q == p)
                .map_or(0.0, |k| ct.weights[k])
        }
    };
    let held = |p: usize| &present[slot[&train_idx[p]]];
    Ok(by_cluster
        .into_par_iter()
        .flat_map_iter(|(_, js)| {
            // weight the cluster takes out of the group and of each child, and the thresholds left
            let needs: Vec<f64> = (0..stride)
                .map(|c| {
                    let lost: f64 = js.iter().map(|&j| weight_at(c, t.members[j])).sum();
                    let g = if c == 0 {
                        t
                    } else {
                        &targets[t.children[c - 1]]
                    };
                    groups::threshold(g.f_in, g.total - lost, g.weighted)
                })
                .collect();
            // per k-mer, the weight the cluster contributes at every level
            let fill = |ks: &[u64]| {
                let mut acc = vec![0f64; ks.len() * stride];
                for &j in &js {
                    let p = t.members[j];
                    let w: Vec<f64> = (0..stride).map(|c| weight_at(c, p)).collect();
                    if w.iter().all(|&x| x == 0.0) {
                        continue;
                    }
                    let h = held(p);
                    for (idx, k) in ks.iter().enumerate() {
                        if h.binary_search(k).is_ok() {
                            for c in 0..stride {
                                acc[idx * stride + c] += w[c];
                            }
                        }
                    }
                }
                acc
            };
            let (removed, extras_removed) = (fill(&m.markers), fill(&m.extras));
            js.iter()
                .map(|&j| {
                    let p = t.members[j];
                    let sc = score::loco_score(
                        &m.markers,
                        &m.counts,
                        &removed,
                        &m.extras,
                        &m.extras_count,
                        &extras_removed,
                        &needs,
                        held(p),
                    );
                    (train_idx[p], (js.len(), sc.markers, sc.present, sc.score))
                })
                .collect::<Vec<_>>()
        })
        .collect())
}

/// Sample indices in scoring batches of about `n`. Each batch holds whole genotype clusters (a
/// named cluster's samples go together, at the position of its first one), so leave-one-cluster-
/// out can be computed inside a batch; a batch exceeds `n` only by the rest of its last cluster.
fn score_batches(samples: &[Sample], n: usize) -> Vec<Vec<usize>> {
    let mut units: Vec<Vec<usize>> = Vec::new();
    let mut unit_of: std::collections::HashMap<&str, usize> = Default::default();
    for (i, s) in samples.iter().enumerate() {
        match s.cluster.as_deref().filter(|c| !c.is_empty()) {
            Some(c) => {
                let u = *unit_of.entry(c).or_insert_with(|| {
                    units.push(Vec::new());
                    units.len() - 1
                });
                units[u].push(i);
            }
            None => units.push(vec![i]),
        }
    }
    let mut out: Vec<Vec<usize>> = Vec::new();
    for u in units {
        match out.last_mut() {
            Some(b) if b.len() < n => b.extend(u),
            _ => out.push(u),
        }
    }
    out
}

/// Chunks of at most 32 targets, as indices into `targets`, keeping every parent together with its
/// children (their weights are counted in the same pass).
fn target_chunks(targets: &[Target]) -> Result<Vec<Vec<usize>>> {
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    let n = targets.len();
    let mut parent: Vec<usize> = (0..n).collect();
    for (t, tg) in targets.iter().enumerate() {
        for &c in &tg.children {
            let (a, b) = (find(&mut parent, t), find(&mut parent, c));
            if a != b {
                parent[a.max(b)] = a.min(b);
            }
        }
    }
    let mut comps: std::collections::BTreeMap<usize, Vec<usize>> = Default::default();
    for i in 0..n {
        let r = find(&mut parent, i);
        comps.entry(r).or_default().push(i);
    }
    let mut out: Vec<Vec<usize>> = Vec::new();
    for (_, c) in comps {
        ensure!(
            c.len() <= 32,
            "{} groups are linked by child constraints; at most 32 fit in one pass",
            c.len()
        );
        match out.last_mut() {
            Some(last) if last.len() + c.len() <= 32 => last.extend(c),
            _ => out.push(c),
        }
    }
    Ok(out)
}

/// Group counting for every target (in passes of up to 32), then the filters, one group at a time.
/// Also returns the hash-sampled training k-mer set and its per-training-sample counts.
fn build_markers(
    cfg: &Config,
    dbs: &[PathBuf],
    targets: &[Target],
) -> Result<(Vec<GroupMarkers>, Vec<u64>, groups::SketchCounts)> {
    let gdir = cfg.tmp().join("groups");
    fs::create_dir_all(&gdir)?;
    let mut done: Vec<Option<GroupMarkers>> = (0..targets.len()).map(|_| None).collect();
    let (mut sketch, mut sketch_counts) = (Vec::new(), groups::SketchCounts::default());
    for (pass, idxs) in target_chunks(targets)?.iter().enumerate() {
        // children are indices into the full slice: remap them to positions inside this chunk
        let pos: std::collections::HashMap<usize, usize> =
            idxs.iter().enumerate().map(|(k, &i)| (i, k)).collect();
        let chunk: Vec<Target> = idxs
            .iter()
            .map(|&i| {
                let mut t = targets[i].clone();
                t.children = t.children.iter().map(|c| pos[c]).collect();
                t
            })
            .collect();
        let chunk = &chunk[..];
        let files = |t: &Target, what: &str| gdir.join(format!("{}.{what}.npy", sanitize(&t.name)));
        // in-group weights are stored as f64 bit patterns in the u64 files
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<u64>>();
        let mut w: Vec<[npy::Writer; 5]> = chunk
            .iter()
            .map(|t| {
                Ok([
                    npy::Writer::create(&files(t, "strict"))?,
                    npy::Writer::create(&files(t, "count"))?,
                    npy::Writer::create(&files(t, "extra"))?,
                    npy::Writer::create(&files(t, "xcount"))?,
                    npy::Writer::create(&files(t, "out"))?,
                ])
            })
            .collect::<Result<_>>()?;
        eprintln!(
            "[groups] {} groups over {} training samples",
            chunk.len(),
            dbs.len()
        );
        let thr = (pass == 0).then(|| crate::kmer::sketch_threshold(SKETCH_SCALE));
        let (stats, sc) = groups::run_sketch(dbs, chunk, true, 0, thr, |b| {
            sketch.extend(&b.sketch);
            for (t, ws) in w.iter_mut().enumerate() {
                ws[0].write(&b.strict[t])?;
                ws[1].write(&bits(&b.strict_count[t]))?;
                ws[2].write(&b.extra[t])?;
                ws[3].write(&bits(&b.extra_count[t]))?;
                ws[4].write(&b.out[t])?;
            }
            Ok(())
        })?;
        for ws in w {
            for x in ws {
                x.finish()?;
            }
        }
        if pass == 0 {
            sketch_counts = sc;
        }
        eprintln!(
            "[mem] peak after group counting pass {}: {}",
            pass + 1,
            peak_rss()
        );
        for ((t, st), &orig) in chunk.iter().zip(stats).zip(idxs) {
            let stride = 1 + t.children.len();
            eprintln!(
                "[filter] {}: {} candidates, {} strict",
                t.name, st.candidates, st.strict
            );
            let floats = |v: Vec<u64>| v.into_iter().map(f64::from_bits).collect::<Vec<f64>>();
            let mut strict = npy::read_u64(&files(t, "strict"))?;
            let mut count = floats(npy::read_u64(&files(t, "count"))?);
            let mut extra = npy::read_u64(&files(t, "extra"))?;
            let mut xcount = floats(npy::read_u64(&files(t, "xcount"))?);
            // cap: keep a hash sample of the candidates (and extras, same hash), so scores and
            // leave-one-out shares stay unbiased estimates
            let take = |v: &[f64], idx: &[usize]| -> Vec<f64> {
                idx.iter()
                    .flat_map(|&i| v[i * stride..(i + 1) * stride].to_vec())
                    .collect()
            };
            let cap = cfg.markers.max_candidates;
            let kept_share = if cap > 0 && strict.len() > cap {
                let share = cap as f64 / strict.len() as f64;
                let thr = (share * u64::MAX as f64) as u64;
                let keep = |k: u64| crate::kmer::mix(k) < thr;
                let idx: Vec<usize> = (0..strict.len()).filter(|&i| keep(strict[i])).collect();
                count = take(&count, &idx);
                strict = idx.iter().map(|&i| strict[i]).collect();
                let idx: Vec<usize> = (0..extra.len()).filter(|&i| keep(extra[i])).collect();
                xcount = take(&xcount, &idx);
                extra = idx.iter().map(|&i| extra[i]).collect();
                eprintln!(
                    "[filter] {}: capped to a {share:.3} hash sample ({} candidates)",
                    t.name,
                    strict.len()
                );
                share
            } else {
                1.0
            };
            let out = npy::read_u64(&files(t, "out"))?;
            let (keep, hamming_survive, keep_x) = if cfg.markers.near_miss {
                let (keep, survive) = filter(&strict, &out);
                (keep, survive, filter(&extra, &out).0)
            } else {
                (
                    (0..strict.len()).collect(),
                    strict.len(),
                    (0..extra.len()).collect(),
                )
            };
            drop(out);
            for what in ["strict", "count", "extra", "xcount", "out"] {
                fs::remove_file(files(t, what)).ok();
            }
            done[orig] = Some(GroupMarkers {
                markers: keep.iter().map(|&i| strict[i]).collect(),
                counts: take(&count, &keep),
                extras: keep_x.iter().map(|&i| extra[i]).collect(),
                extras_count: take(&xcount, &keep_x),
                stride,
                stats: st,
                hamming_survive,
                extras_raw: extra.len(),
                kept_share,
            });
        }
    }
    let all = done
        .into_iter()
        .map(|g| g.context("a target was not counted"))
        .collect::<Result<Vec<_>>>()?;
    Ok((all, sketch, sketch_counts))
}

/// UNEXPLAINED threshold from the training samples' leave-one-out shares: twice the 99th
/// percentile, at least 0.05 and at most 0.5.
///
/// It used to be twice the largest, which one contaminated training sample can push above 1, so
/// the flag never fires (census 17 Oct: 4 foreign runs took it to 1.96). A percentile cannot be
/// moved by a few samples, and the ceiling keeps the flag alive whatever the training set holds:
/// pure isolates sit far below (census p99 0.025, pilot max 0.053) and a foreign genome far above
/// (44-96%), so any threshold in between separates them.
pub fn unexplained_threshold(loo: &[f64; 5]) -> f64 {
    (2.0 * loo[3]).clamp(0.05, 0.5)
}

/// Nearest-rank quantile (index round((n - 1) f) of the sorted values), as for unexplained_loo.
fn quantile(v: &[f64], f: f64) -> Option<f64> {
    let mut v = v.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (!v.is_empty()).then(|| v[((v.len() - 1) as f64 * f).round() as usize])
}

fn median(mut v: Vec<f64>) -> [f64; 3] {
    if v.is_empty() {
        return [f64::NAN; 3];
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let m = if v.len() % 2 == 1 {
        v[v.len() / 2]
    } else {
        (v[v.len() / 2 - 1] + v[v.len() / 2]) / 2.0
    };
    [v[0], m, v[v.len() - 1]]
}

/// An integer as an integer, anything else to 4 decimals.
fn fnum(x: f64) -> String {
    if x.fract() == 0.0 {
        format!("{x:.0}")
    } else {
        format!("{x:.4}")
    }
}

fn f4(x: f64) -> String {
    if x.is_nan() {
        "NA".into()
    } else {
        format!("{x:.4}")
    }
}

pub fn run(cfg: &Config) -> Result<Bundle> {
    let samples = load_samples(&cfg.samples, &cfg.levels)?;
    let hier = hierarchy(&samples, &cfg.levels)?;
    fs::create_dir_all(&cfg.output)?;
    eprintln!(
        "{} samples, {} groups on {} levels",
        samples.len(),
        hier.len(),
        cfg.levels.len()
    );

    // 1-2: counting and training set
    let stats = count_all(cfg, &samples)?;
    let solid_ref = solid_reference(&samples, &stats);
    let excess: Vec<f64> = samples
        .iter()
        .zip(&stats)
        .map(|(s, st)| {
            kmer_excess(
                &solid_ref,
                s.labels.first().and_then(|l| l.as_deref()),
                st.solid_kmers,
            )
        })
        .collect();
    let train: Vec<(bool, &str)> = samples
        .iter()
        .zip(&stats)
        .zip(&excess)
        .map(|((s, st), &x)| training(s, st, x, &cfg.training))
        .collect();
    let train_idx: Vec<usize> = (0..samples.len()).filter(|&i| train[i].0).collect();
    ensure!(train_idx.len() >= 2, "fewer than 2 training samples");
    let pos_in_train = |i: usize| train_idx.iter().position(|&j| j == i);
    let dbs: Vec<PathBuf> = train_idx
        .iter()
        .map(|&i| db_path(cfg, &samples[i]))
        .collect();
    let mut qc = BufWriter::new(fs::File::create(cfg.output.join("qc.tsv"))?);
    writeln!(
        qc,
        "sample\tinput\tcutoff\tpeak\tsolid_kmers\terror_kmers\tlow_coverage\tkmer_excess\ttrain\treason\tcluster"
    )?;
    for (((s, st), (t, why)), x) in samples.iter().zip(&stats).zip(&train).zip(&excess) {
        let input = if matches!(s.input, Input::Reads(_)) {
            "reads"
        } else {
            "assembly"
        };
        let peak = st.peak.map_or("NA".into(), |p| p.to_string());
        writeln!(
            qc,
            "{}\t{input}\t{}\t{peak}\t{}\t{}\t{}\t{}\t{}\t{why}\t{}",
            s.name,
            st.cutoff,
            st.solid_kmers,
            st.error_kmers,
            st.low_coverage,
            f4(*x),
            t,
            s.cluster.as_deref().unwrap_or("")
        )?;
    }
    qc.flush()?;

    eprintln!("[mem] peak after counting and QC: {}", peak_rss());
    // 3-4: groups with enough training samples, then markers
    let mut built: Vec<&Group> = Vec::new();
    let mut targets = Vec::new();
    let mut skipped = Vec::new();
    for g in &hier {
        let tm: Vec<usize> = g.members.iter().filter_map(|&i| pos_in_train(i)).collect();
        if tm.len() < cfg.markers.min_group_size {
            skipped.push(format!("{} ({} training)", g.id(), tm.len()));
            continue;
        }
        // clade weights within the group: training members sharing a cluster split one weight
        let clusters: Vec<Option<&str>> = tm
            .iter()
            .map(|&p| samples[train_idx[p]].cluster.as_deref())
            .collect();
        let w = groups::cluster_weights(&clusters);
        // training samples whose label is unknown at this level but whose known labels match the
        // group's ancestors may belong to it: neither in nor out
        let ignore: Vec<usize> = (0..train_idx.len())
            .filter(|&p| {
                let labels = &samples[train_idx[p]].labels;
                let known = labels.iter().take_while(|l| l.is_some()).count();
                known <= g.level
                    && labels[..known]
                        .iter()
                        .zip(&g.path)
                        .all(|(l, q)| l.as_deref() == Some(q.as_str()))
            })
            .collect();
        targets.push(
            Target::weighted(&g.id(), tm, w, cfg.markers.f_in, cfg.markers.max_out)
                .with_ignore(ignore),
        );
        built.push(g);
    }
    // hierarchy-aware prevalence: each built group's children are the built groups one level down
    // whose path extends it (children of skipped groups simply constrain nothing)
    if cfg.markers.balance_children {
        for t in 0..built.len() {
            let g = built[t];
            let kids: Vec<usize> = built
                .iter()
                .enumerate()
                .filter(|(_, c)| {
                    c.level == g.level + 1
                        && c.path.len() == g.path.len() + 1
                        && c.path[..g.path.len()] == g.path[..]
                })
                .map(|(i, _)| i)
                .collect();
            if !kids.is_empty() {
                eprintln!(
                    "[groups] {} must also reach f_in in: {}",
                    g.id(),
                    kids.iter()
                        .map(|&k| built[k].id())
                        .collect::<Vec<_>>()
                        .join(", ")
                );
                let t0 = std::mem::replace(&mut targets[t], Target::new("", vec![], 1.0));
                targets[t] = t0.with_children(kids);
            }
        }
    }
    for s in &skipped {
        eprintln!("skipped group {s}");
    }
    ensure!(
        !targets.is_empty(),
        "no group has {} training samples",
        cfg.markers.min_group_size
    );
    let (gm, sketch, sketch_counts) = build_markers(cfg, &dbs, &targets)?;

    eprintln!("[mem] peak after markers: {}", peak_rss());
    // 5: bundle files
    let bdir = cfg.output.join("bundle");
    let mut bgroups = Vec::new();
    let mut gtsv = BufWriter::new(fs::File::create({
        fs::create_dir_all(&bdir)?;
        bdir.join("groups.tsv")
    })?);
    writeln!(
        gtsv,
        "group\tlevel\tn_train\tneed\tcandidates\tstrict\tkept_share\thamming_survive\tmarkers\tunitigs\textras_raw\textras\tn_eff\tout0\tn_ignored\tn_children"
    )?;
    for ((g, t), m) in built.iter().zip(&targets).zip(&gm) {
        let level = &cfg.levels[g.level];
        let file = format!(
            "markers/{}/{}",
            sanitize(level),
            sanitize(&g.id().replace('/', "__"))
        );
        let stem = bdir.join(&file);
        fs::create_dir_all(stem.parent().unwrap())?;
        kmc::write_kmc1(&stem, &m.markers, 31, 7)?;
        let u = unitig::unitigs(&m.markers);
        let mut fa = BufWriter::new(fs::File::create(with_ext(&stem, ".fa"))?);
        for (i, x) in u.iter().enumerate() {
            writeln!(fa, ">{}_u{i}_n{}\n{}", sanitize(g.name()), x.n_kmers, x.seq)?;
        }
        fa.flush()?;
        writeln!(
            gtsv,
            "{}\t{level}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            g.id(),
            t.members.len(),
            fnum(t.need),
            m.stats.candidates,
            m.stats.strict,
            f4(m.kept_share),
            m.hamming_survive,
            m.markers.len(),
            u.len(),
            m.extras_raw,
            m.extras.len(),
            fnum(t.total),
            m.stats.out0,
            t.ignore.len(),
            t.children.len()
        )?;
        bgroups.push(BundleGroup {
            id: g.id(),
            level: g.level,
            name: g.name().to_string(),
            file,
            n_train: t.members.len(),
            n_eff: t.total,
            need: t.need,
            markers: m.markers.len(),
            loo: [f64::NAN; 3],
            loo_p05: None,
            loco: [f64::NAN; 3],
        });
    }
    gtsv.flush()?;

    eprintln!("[mem] peak after bundle files: {}", peak_rss());
    // 6: scores for every sample; leave-one-out for training members of their own groups
    let mut union: Vec<u64> = gm
        .iter()
        .flat_map(|m| m.markers.iter().chain(&m.extras))
        .copied()
        .collect();
    union.sort_unstable();
    union.dedup();
    eprintln!(
        "[score] {} samples against {} marker k-mers",
        samples.len(),
        union.len()
    );
    let thr = crate::kmer::sketch_threshold(SKETCH_SCALE);
    let mut sc = BufWriter::new(fs::File::create(cfg.output.join("scores.tsv"))?);
    let mut loo = BufWriter::new(fs::File::create(cfg.output.join("loo.tsv"))?);
    writeln!(
        sc,
        "sample\ttrain\tunexplained\t{}",
        built.iter().map(|g| g.id()).collect::<Vec<_>>().join("\t")
    )?;
    writeln!(
        loo,
        "sample\tgroup\tmarkers\tpresent\tsensitivity\tcluster\tcluster_size\tloco_markers\tloco_present\tloco_sensitivity"
    )?;
    let mut loo_by_group: Vec<Vec<f64>> = vec![Vec::new(); built.len()];
    let mut loco_by_group: Vec<Vec<f64>> = vec![Vec::new(); built.len()];
    // Scoring streams over the samples a batch at a time: memory holds one batch's present marker
    // k-mers, not every sample's (the census: 3,736 samples x ~310k present = 9 GB). A batch holds
    // whole genotype clusters, since leave-one-cluster-out needs all of a cluster's training
    // members at once. Rows are kept per sample and written in sample-sheet order.
    let mut unexplained = vec![0f64; samples.len()];
    let mut sc_rows = vec![String::new(); samples.len()];
    let mut loo_rows = vec![String::new(); samples.len()];
    for batch in score_batches(&samples, cfg.resources.score_batch.max(1)) {
        let read: Vec<(Vec<u64>, u64, u64)> = batch
            .par_iter()
            .map(|&i| kmc::read_present_sketch(&db_path(cfg, &samples[i]), &union, thr, &sketch))
            .collect::<Result<_>>()?;
        let slot: std::collections::HashMap<usize, usize> =
            batch.iter().enumerate().map(|(b, &i)| (i, b)).collect();
        let (present, sampled): (Vec<Vec<u64>>, Vec<(u64, u64)>) =
            read.into_iter().map(|(p, n, k)| (p, (n, k))).unzip();
        // leave-one-cluster-out: per group, the multi-member clusters of this batch (a cluster of
        // one gives exactly the leave-one-out answer, so it is reused)
        let loco: Vec<std::collections::HashMap<usize, LocoResult>> = (0..targets.len())
            .map(|ti| loco_scores(ti, &targets, &gm[ti], &samples, &train_idx, &present, &slot))
            .collect::<Result<_>>()?;
        type Row = (f64, String, String, Vec<(usize, f64, f64)>);
        let rows: Vec<Row> = batch
            .par_iter()
            .zip(present.par_iter().zip(&sampled))
            .map(|(&i, (present, &(n_sampled, n_known)))| {
                let s = &samples[i];
                // unexplained share: training samples leave-one-out (k-mers in no other training
                // sample), the rest against all training samples
                let u = match pos_in_train(i) {
                    Some(p) => {
                        sketch_counts.unique[p] as f64 / sketch_counts.total[p].max(1) as f64
                    }
                    None => 1.0 - n_known as f64 / n_sampled.max(1) as f64,
                };
                let mut row = vec![s.name.clone(), train[i].0.to_string(), f4(u)];
                let mut lines = String::new();
                let mut own = Vec::new(); // (group, leave-one-out, leave-one-cluster-out)
                for (gi, ((g, t), m)) in built.iter().zip(&targets).zip(&gm).enumerate() {
                    if m.markers.is_empty() {
                        row.push("NA".into()); // no markers: the score is undefined, not 0
                        continue;
                    }
                    let j = pos_in_train(i).and_then(|p| t.members.iter().position(|&q| q == p));
                    let r = if let Some(j) = j {
                        // thresholds of each child group with this sample held out
                        let rules: Vec<score::ChildRule> = t
                            .children
                            .iter()
                            .map(|&c| {
                                let ct = &targets[c];
                                let has = pos_in_train(i)
                                    .and_then(|p| ct.members.iter().position(|&q| q == p));
                                score::ChildRule {
                                    need: ct.need,
                                    need_without_x: has.map_or(ct.need, |k| ct.need_without(k)),
                                    has_x: has.is_some(),
                                }
                            })
                            .collect();
                        debug_assert_eq!(
                            m.stride,
                            1 + rules.len(),
                            "one weight per child, plus the group"
                        );
                        let r = score::loo_score(
                            &m.markers,
                            &m.counts,
                            &m.extras,
                            &m.extras_count,
                            &rules,
                            present,
                            t.weights[j],
                            t.need_without(j),
                        );
                        // leave-one-cluster-out: the cluster's own answer, or this run's when alone
                        let key = s.cluster.as_deref().unwrap_or("");
                        let cl = loco[gi]
                            .get(&i)
                            .copied()
                            .unwrap_or((1, r.markers, r.present, r.score));
                        lines.push_str(&format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                            s.name,
                            g.id(),
                            r.markers,
                            r.present,
                            f4(r.score),
                            if key.is_empty() { "-" } else { key },
                            cl.0,
                            cl.1,
                            cl.2,
                            f4(cl.3)
                        ));
                        own.push((gi, r.score, cl.3));
                        r
                    } else {
                        score::score(&m.markers, &[], present, false)
                    };
                    row.push(f4(r.score));
                }
                (u, row.join("\t"), lines, own)
            })
            .collect();
        for (&i, (u, row, lines, own)) in batch.iter().zip(rows) {
            unexplained[i] = u;
            sc_rows[i] = row;
            loo_rows[i] = lines;
            for (gi, l, c) in own {
                loo_by_group[gi].push(l);
                loco_by_group[gi].push(c);
            }
        }
    }
    for (row, lines) in sc_rows.iter().zip(&loo_rows) {
        writeln!(sc, "{row}")?;
        loo.write_all(lines.as_bytes())?;
    }
    sc.flush()?;
    loo.flush()?;
    eprintln!("[mem] peak after scoring: {}", peak_rss());
    for ((bg, v), w) in bgroups.iter_mut().zip(loo_by_group).zip(loco_by_group) {
        bg.loo_p05 = quantile(&v, 0.05);
        bg.loo = median(v);
        bg.loco = median(w);
    }

    npy::write_u64(&bdir.join("sketch.u64.npy"), &sketch)?;
    let mut u: Vec<f64> = train_idx.iter().map(|&i| unexplained[i]).collect();
    u.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |f: f64| u[((u.len() - 1) as f64 * f).round() as usize];
    let unexplained_loo = [u[0], q(0.5), q(0.95), q(0.99), u[u.len() - 1]];
    let bundle = Bundle {
        format: BUNDLE_FORMAT,
        cladeprint_version: env!("CARGO_PKG_VERSION").to_string(),
        k: cfg.kmers.k,
        levels: cfg.levels.clone(),
        min_peak_coverage: cfg.kmers.min_peak_coverage,
        groups: bgroups,
        skipped,
        sketch_scale: SKETCH_SCALE,
        sketch_file: "sketch.u64.npy".into(),
        unexplained_loo,
        unexplained_flag: unexplained_threshold(&unexplained_loo),
        solid_reference: solid_ref,
        max_kmer_excess: cfg.training.max_kmer_excess,
        config: cfg.clone(),
    };
    fs::write(
        bdir.join("bundle.toml"),
        toml::to_string(&bundle).context("writing bundle.toml")?,
    )?;
    fs::remove_dir_all(cfg.tmp().join("groups")).ok();
    Ok(bundle)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn score_batches_keep_clusters_whole() {
        let s = |name: &str, cluster: &str| Sample {
            name: name.into(),
            input: Input::Assembly(PathBuf::from("x.fa")),
            labels: vec![],
            train: Train::Auto,
            cluster: (!cluster.is_empty()).then(|| cluster.to_string()),
        };
        // a, c, e share cluster x; b, d share y; f is alone
        let samples = [
            s("a", "x"),
            s("b", "y"),
            s("c", "x"),
            s("d", "y"),
            s("e", "x"),
            s("f", ""),
        ];
        assert_eq!(
            score_batches(&samples, 1),
            vec![vec![0, 2, 4], vec![1, 3], vec![5]]
        );
        assert_eq!(
            score_batches(&samples, 4),
            vec![vec![0, 2, 4, 1, 3], vec![5]]
        );
        let mut all: Vec<usize> = score_batches(&samples, 100).concat();
        all.sort_unstable();
        assert_eq!(all, (0..6).collect::<Vec<_>>());
    }
}
