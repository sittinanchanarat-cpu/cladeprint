//! `cladeprint classify`: score a sample against a marker bundle and call it level by level.
//! Thresholds are first values, to be calibrated on the external test
//! set (pure isolates: off-target scores <= 0.12; hybrids: both parents >= 0.38).
//!
//! Per group: score = share of the group's markers present in the sample. Calls go top-down:
//! level 0 takes the best-scoring group, each lower level only the children of the call above.
//! A group is called when its score reaches `CALL_MIN`; no call at a level stops the descent.
//! Confidence: high when the score reaches the group's 5th-percentile leave-one-out score (the
//! minimum in bundles without it) and nothing else at the level is MIXED; medium when it reaches
//! half the leave-one-out median; else low. A low-confidence call also stops the descent: the
//! R2 holdouts showed an unsampled clade lands on its nearest sibling, and calls below that
//! were confidently wrong (31 VNIa-high under a held-out VNIa sublineage). Not at LOW_COVERAGE:
//! there a low score means missing k-mers, not a novel clade (on the external test set the stop
//! would drop 21 right low-coverage calls and fix none). Also at LOW_COVERAGE, a call is no more
//! confident than the call above it: groups with few markers in few unitigs (VNIa 303 in 21, VNIc
//! 146 in 14) gain or lose markers a unitig at a time and passed their high threshold at 0.4-0.8x
//! while species and lineage were low or medium (coverage calibration, 22 Oct: 6 rows, all right).
//! Flags: MIXED (another group at the level, not an ancestor or descendant of the call, scores at
//! least MIXED_MIN), LOW_COVERAGE (from counting), NO_CALL (nothing reached CALL_MIN),
//! NOT_DESCENDED (the level has groups under the call above, which was low confidence).
//! UNEXPLAINED (sample level): the share of the sample's k-mers found in no training sample,
//! estimated on a hash sample, is above the bundle's threshold (set from the training samples'
//! leave-one-out shares). Something in the sample is outside the trained groups: an untrained
//! lineage or species, a hybrid with one, or contamination. At LOW_COVERAGE the flag is
//! UNEXPLAINED?: sequencing errors that pass the cutoff also count as unexplained there.
//! EXCESS_KMERS (sample level, in the CLI): more solid k-mers than `max_kmer_excess` x the
//! called top-level group's reference: two genomes (mixture, heterozygous hybrid) or
//! contamination; EXCESS_KMERS? at low coverage.

use crate::kmc;
use crate::pipeline::{Bundle, BundleGroup};
use crate::score::n_shared;
use anyhow::{Context, Result};
use std::path::Path;

pub const CALL_MIN: f64 = 0.3;
pub const MIXED_MIN: f64 = 0.3;
/// reads classification: fewer solid k-mers than this share of the smallest reference genome
/// means the histogram valley was in the repeat tail, not between errors and the genome (below
/// 1x; low-coverage samples with a cutoff of 2 keep 0.36-0.48, the misread valleys 0.002-0.004)
pub const MIN_GENOME_SHARE: f64 = 0.1;

pub struct LoadedBundle {
    pub bundle: Bundle,
    /// markers of each group (sorted), in bundle.groups order
    pub markers: Vec<Vec<u64>>,
    /// union of all markers (sorted)
    pub union: Vec<u64>,
    /// hash sample of the training k-mers (sorted) and its threshold
    pub sketch: Vec<u64>,
    pub sketch_threshold: u64,
}

impl LoadedBundle {
    /// The sample's marker k-mers and its unexplained share (share of its hash-sampled k-mers in
    /// no training sample), from one pass over its presence database.
    pub fn read(&self, db: &Path) -> Result<(Vec<u64>, f64)> {
        let (present, sampled, known) =
            kmc::read_present_sketch(db, &self.union, self.sketch_threshold, &self.sketch)?;
        Ok((present, 1.0 - known as f64 / sampled.max(1) as f64))
    }
}

pub fn load(dir: &Path) -> Result<LoadedBundle> {
    let text = std::fs::read_to_string(dir.join("bundle.toml"))
        .with_context(|| format!("reading {}/bundle.toml", dir.display()))?;
    // check the format first, so an old bundle gives a clear message, not a missing-field error
    let format = text
        .parse::<toml::Table>()
        .ok()
        .and_then(|t| t.get("format")?.as_integer());
    anyhow::ensure!(
        format == Some(crate::pipeline::BUNDLE_FORMAT as i64),
        "bundle format {:?}, this cladeprint reads format {}: rebuild the bundle with `cladeprint run`",
        format,
        crate::pipeline::BUNDLE_FORMAT
    );
    let bundle: Bundle = toml::from_str(&text).context("parsing bundle.toml")?;
    anyhow::ensure!(
        bundle.format == crate::pipeline::BUNDLE_FORMAT,
        "bundle format {} not supported",
        bundle.format
    );
    let markers: Vec<Vec<u64>> = bundle
        .groups
        .iter()
        .map(|g| kmc::read_sorted(&dir.join(&g.file)))
        .collect::<Result<_>>()?;
    let mut union: Vec<u64> = markers.iter().flatten().copied().collect();
    union.sort_unstable();
    union.dedup();
    let sketch = crate::npy::read_u64(&dir.join(&bundle.sketch_file))?;
    let sketch_threshold = crate::kmer::sketch_threshold(bundle.sketch_scale);
    Ok(LoadedBundle {
        bundle,
        markers,
        union,
        sketch,
        sketch_threshold,
    })
}

#[derive(Debug, Clone)]
pub struct GroupScore {
    pub id: String,
    pub level: usize,
    /// None when the group has no markers
    pub score: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LevelCall {
    pub level: String,
    pub call: Option<String>,
    pub score: Option<f64>,
    pub confidence: &'static str,
    pub flags: Vec<String>,
}

fn related(a: &str, b: &str) -> bool {
    a == b || a.starts_with(&format!("{b}/")) || b.starts_with(&format!("{a}/"))
}

/// Scores of every group from the sample's marker k-mers (sorted).
pub fn scores(b: &LoadedBundle, present: &[u64]) -> Vec<GroupScore> {
    b.bundle
        .groups
        .iter()
        .zip(&b.markers)
        .map(|(g, m)| GroupScore {
            id: g.id.clone(),
            level: g.level,
            score: (!m.is_empty()).then(|| n_shared(m, present) as f64 / m.len() as f64),
        })
        .collect()
}

/// Level-by-level calls from group scores.
pub fn call(
    groups: &[BundleGroup],
    levels: &[String],
    sc: &[GroupScore],
    low_coverage: bool,
) -> Vec<LevelCall> {
    let mut out = Vec::new();
    let mut parent: Option<String> = None;
    // confidence of the call above (caps this level's at LOW_COVERAGE)
    let mut parent_conf: Option<&'static str> = None;
    // the call above was low confidence: report this level, but do not call it
    let mut stopped = false;
    for (l, name) in levels.iter().enumerate() {
        let at_level: Vec<(usize, f64)> = sc
            .iter()
            .enumerate()
            .filter(|(_, s)| s.level == l)
            .filter_map(|(i, s)| s.score.map(|x| (i, x)))
            .collect();
        let candidates: Vec<(usize, f64)> = at_level
            .iter()
            .copied()
            .filter(|&(i, _)| match &parent {
                None => l == 0,
                Some(p) => sc[i].id.starts_with(&format!("{p}/")),
            })
            .collect();
        let mut flags = Vec::new();
        if low_coverage {
            flags.push("LOW_COVERAGE".to_string());
        }
        let best = candidates
            .iter()
            .copied()
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        if stopped {
            if !candidates.is_empty() {
                flags.push("NOT_DESCENDED".to_string());
            }
            out.push(LevelCall {
                level: name.clone(),
                call: None,
                score: best.map(|b| b.1),
                confidence: "none",
                flags,
            });
            stopped = false;
            parent = None;
            parent_conf = None;
            continue;
        }
        let Some((bi, bs)) = best.filter(|&(_, s)| s >= CALL_MIN) else {
            // NO_CALL only when there was something to call: a called parent without child groups
            // in the bundle (e.g. no sublineages defined under VNBI) leaves the level blank
            if !candidates.is_empty() {
                flags.push("NO_CALL".to_string());
            }
            out.push(LevelCall {
                level: name.clone(),
                call: None,
                score: best.map(|b| b.1),
                confidence: "none",
                flags,
            });
            parent = None;
            parent_conf = None;
            continue;
        };
        let id = sc[bi].id.clone();
        let others: Vec<String> = at_level
            .iter()
            .filter(|&&(i, s)| i != bi && s >= MIXED_MIN && !related(&sc[i].id, &id))
            .map(|&(i, s)| format!("MIXED:{}={s:.2}", sc[i].id))
            .collect();
        let g = &groups[bi];
        let mut confidence = if bs >= g.loo_p05.unwrap_or(g.loo[0]) && others.is_empty() {
            "high"
        } else if bs >= 0.5 * g.loo[1] {
            "medium"
        } else {
            "low"
        };
        let rank = |c: &str| ["low", "medium", "high"].iter().position(|&x| x == c);
        if let Some(pc) = parent_conf.filter(|&pc| low_coverage && rank(confidence) > rank(pc)) {
            confidence = pc;
        }
        stopped = confidence == "low" && !low_coverage;
        flags.extend(others);
        out.push(LevelCall {
            level: name.clone(),
            call: Some(id.clone()),
            score: Some(bs),
            confidence,
            flags,
        });
        parent = Some(id);
        parent_conf = Some(confidence);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bg(id: &str, level: usize, loo: [f64; 3]) -> BundleGroup {
        BundleGroup {
            id: id.into(),
            level,
            name: id.rsplit('/').next().unwrap().into(),
            file: String::new(),
            n_train: 10,
            n_eff: 10.0,
            need: 9.0,
            markers: 100,
            loo,
            loo_p05: None,
            loco: loo,
        }
    }

    fn gs(id: &str, level: usize, score: f64) -> GroupScore {
        GroupScore {
            id: id.into(),
            level,
            score: Some(score),
        }
    }

    #[test]
    fn calls_descend_and_flag_mixtures() {
        let levels = vec!["lineage".to_string(), "sub".to_string()];
        let groups = vec![
            bg("VNI", 0, [0.9, 0.97, 1.0]),
            bg("VNII", 0, [0.8, 0.98, 1.0]),
            bg("VNI/VNIa", 1, [0.7, 0.98, 1.0]),
            bg("VNI/VNIb", 1, [0.6, 0.91, 1.0]),
        ];
        // a clean VNIa
        let sc = vec![
            gs("VNI", 0, 0.96),
            gs("VNII", 0, 0.0),
            gs("VNI/VNIa", 1, 0.95),
            gs("VNI/VNIb", 1, 0.01),
        ];
        let c = call(&groups, &levels, &sc, false);
        assert_eq!(
            (c[0].call.as_deref(), c[0].confidence),
            (Some("VNI"), "high")
        );
        assert_eq!(
            (c[1].call.as_deref(), c[1].confidence),
            (Some("VNI/VNIa"), "high")
        );
        assert!(c.iter().all(|x| x.flags.is_empty()));
        // a VNI x VNII hybrid: called VNI, flagged MIXED, so not high confidence
        let sc = vec![
            gs("VNI", 0, 0.8),
            gs("VNII", 0, 0.7),
            gs("VNI/VNIa", 1, 0.2),
            gs("VNI/VNIb", 1, 0.1),
        ];
        let c = call(&groups, &levels, &sc, true);
        assert_eq!(c[0].call.as_deref(), Some("VNI"));
        assert_ne!(c[0].confidence, "high");
        assert!(c[0].flags.contains(&"LOW_COVERAGE".to_string()));
        assert!(c[0].flags.iter().any(|f| f.starts_with("MIXED:VNII=")));
        // no sublineage reaches CALL_MIN
        assert_eq!(c[1].call, None);
        assert!(c[1].flags.contains(&"NO_CALL".to_string()));
        // nothing at all: no call at the top stops the descent without a second NO_CALL
        let sc = vec![
            gs("VNI", 0, 0.01),
            gs("VNII", 0, 0.02),
            gs("VNI/VNIa", 1, 0.0),
            gs("VNI/VNIb", 1, 0.0),
        ];
        let c = call(&groups, &levels, &sc, false);
        assert_eq!((c[0].call.clone(), c[1].call.clone()), (None, None));
        assert!(
            c[0].flags.contains(&"NO_CALL".to_string())
                && !c[1].flags.contains(&"NO_CALL".to_string())
        );
        // a lineage without sublineages in the bundle: the lower level is blank, not NO_CALL
        let sc = vec![
            gs("VNI", 0, 0.0),
            gs("VNII", 0, 0.97),
            gs("VNI/VNIa", 1, 0.0),
            gs("VNI/VNIb", 1, 0.0),
        ];
        let c = call(&groups, &levels, &sc, false);
        assert_eq!(c[0].call.as_deref(), Some("VNII"));
        assert_eq!(c[1].call, None);
        assert!(c[1].flags.is_empty(), "{:?}", c[1].flags);
    }

    #[test]
    fn low_confidence_stops_the_descent() {
        let levels = vec!["species".into(), "lineage".into(), "sub".into()];
        let groups = vec![
            bg("Cn", 0, [0.9, 0.97, 1.0]),
            bg("Cn/VNI", 1, [0.9, 0.97, 1.0]),
            bg("Cn/VNI/VNIa", 2, [0.7, 0.98, 1.0]),
        ];
        // VNI at 0.4 < 0.5 x 0.97: low, so VNIa is not called although it scores 0.99
        let sc = vec![
            gs("Cn", 0, 0.96),
            gs("Cn/VNI", 1, 0.4),
            gs("Cn/VNI/VNIa", 2, 0.99),
        ];
        let c = call(&groups, &levels, &sc, false);
        assert_eq!(
            (c[1].call.as_deref(), c[1].confidence),
            (Some("Cn/VNI"), "low")
        );
        assert_eq!((c[2].call.clone(), c[2].confidence), (None, "none"));
        assert_eq!(c[2].flags, vec!["NOT_DESCENDED".to_string()]);
        // at low coverage a low call still descends: low there means missing k-mers
        let c = call(&groups, &levels, &sc, true);
        assert_eq!(c[1].confidence, "low");
        assert_eq!(c[2].call.as_deref(), Some("Cn/VNI/VNIa"));
        // ... and is capped at its parent's confidence: 0.99 alone would be high
        assert_eq!(c[2].confidence, "low");
        // medium still descends
        let sc = vec![
            gs("Cn", 0, 0.96),
            gs("Cn/VNI", 1, 0.6),
            gs("Cn/VNI/VNIa", 2, 0.99),
        ];
        let c = call(&groups, &levels, &sc, false);
        assert_eq!(c[1].confidence, "medium");
        assert_eq!(c[2].call.as_deref(), Some("Cn/VNI/VNIa"));
        assert_eq!(c[2].confidence, "high", "no cap with coverage");
        let c = call(&groups, &levels, &sc, true);
        assert_eq!(c[2].confidence, "medium", "capped at low coverage");
        // the cap only lowers: a medium child under a high parent stays medium
        let sc = vec![
            gs("Cn", 0, 0.96),
            gs("Cn/VNI", 1, 0.95),
            gs("Cn/VNI/VNIa", 2, 0.6),
        ];
        let c = call(&groups, &levels, &sc, true);
        assert_eq!((c[1].confidence, c[2].confidence), ("high", "medium"));
    }

    #[test]
    fn high_uses_the_5th_percentile_when_present() {
        let levels = vec!["species".to_string()];
        let mut g = bg("Cn", 0, [0.25, 0.97, 1.0]);
        let sc = vec![gs("Cn", 0, 0.5)];
        // one odd member puts the minimum at 0.25: 0.5 would count as high
        assert_eq!(
            call(&[g.clone()], &levels, &sc, false)[0].confidence,
            "high"
        );
        g.loo_p05 = Some(0.9);
        assert_eq!(
            call(&[g.clone()], &levels, &sc, false)[0].confidence,
            "medium"
        );
        let sc = vec![gs("Cn", 0, 0.92)];
        assert_eq!(call(&[g], &levels, &sc, false)[0].confidence, "high");
    }
}
