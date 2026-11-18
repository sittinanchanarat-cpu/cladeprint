//! Lineage scores and leave-one-out sensitivity (port of `kmers/score_isolates.py`).
//!
//! Markers per target T: `m` (current markers) and `mx` (k-mers in exactly need-1 target isolates
//! that pass the same filters). Holding out a target isolate X leaves the out-group unchanged and
//! lowers the threshold by one, so the markers built without X are m + (mx - X), and X's
//! leave-one-out sensitivity is |m & X| / (|m| + |mx - X|). Exact, no rebuild.
//! Every other isolate is scored as |m & X| / |m|. All sets are sorted, deduplicated u64.

/// |a & b| for sorted, deduplicated slices.
pub fn n_shared(a: &[u64], b: &[u64]) -> usize {
    let (mut i, mut j, mut n) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                n += 1;
                i += 1;
                j += 1;
            }
        }
    }
    n
}

pub struct Score {
    pub score: f64,
    /// markers in the denominator (|m|, or |m| + |mx - X| for a held-out member)
    pub markers: usize,
    pub present: usize,
}

/// Score of one isolate (`present`: its marker k-mers) for one target. `held_out_member`: the
/// isolate is a training member of the target, so the leave-one-out markers apply.
pub fn score(m: &[u64], mx: &[u64], present: &[u64], held_out_member: bool) -> Score {
    let m_in = n_shared(m, present);
    let markers = if held_out_member {
        m.len() + mx.len() - n_shared(mx, present)
    } else {
        m.len()
    };
    Score {
        score: if markers > 0 {
            m_in as f64 / markers as f64
        } else {
            0.0
        },
        markers,
        present: m_in,
    }
}

/// A child group's thresholds for one held-out sample X (hierarchy-aware prevalence).
#[derive(Debug, Clone, Copy)]
pub struct ChildRule {
    /// the child's normal threshold
    pub need: f64,
    /// its threshold with X removed (= `need` when X is not a member)
    pub need_without_x: f64,
    /// X is a member of this child
    pub has_x: bool,
}

/// Leave-one-out score of a training member X of a group, in general (weighted or not, with or
/// without child constraints).
///
/// `m_in` / `extras_in` hold `1 + children.len()` weights per k-mer: the group's in-weight then
/// each child's, in the order of `children` (as `groups::Batch`). Holding X out lowers every
/// threshold and removes X's weight from the k-mers X holds, so the markers without X are exactly
/// the k-mers (markers and extras alike) that satisfy every threshold after that:
///   group:  in - [w_x if X holds it] >= need1
///   child:  in_c - [w_x if X holds it and X is in c] >= that child's threshold without X
/// A k-mer X lacks keeps its weights while the thresholds only fall, so every marker X lacks stays
/// one; an extra X holds can never become one (its weight falls by w_x, more than the threshold
/// falls), which the same predicate expresses. All slices sorted.
#[allow(clippy::too_many_arguments)]
pub fn loo_score(
    m: &[u64],
    m_in: &[f64],
    extras: &[u64],
    extras_in: &[f64],
    children: &[ChildRule],
    present: &[u64],
    w_x: f64,
    need1: f64,
) -> Score {
    let stride = 1 + children.len();
    debug_assert!(m.len() * stride == m_in.len() && extras.len() * stride == extras_in.len());
    let has = |k: &u64| present.binary_search(k).is_ok();
    let keeps = |w: &[f64], idx: usize, held: bool| -> bool {
        let base = idx * stride;
        let drop = if held { w_x } else { 0.0 };
        if w[base] - drop < need1 {
            return false;
        }
        children.iter().enumerate().all(|(c, r)| {
            let d = if held && r.has_x { w_x } else { 0.0 };
            w[base + 1 + c] - d >= if r.has_x { r.need_without_x } else { r.need }
        })
    };
    let (mut kept, mut held_n) = (0usize, 0usize);
    for (idx, k) in m.iter().enumerate() {
        let held = has(k);
        if keeps(m_in, idx, held) {
            kept += 1;
            held_n += held as usize;
        }
    }
    for (idx, k) in extras.iter().enumerate() {
        if !has(k) && keeps(extras_in, idx, false) {
            kept += 1;
        }
    }
    Score {
        score: if kept > 0 {
            held_n as f64 / kept as f64
        } else {
            0.0
        },
        markers: kept,
        present: held_n,
    }
}

/// Leave-one-cluster-out: the sensitivity of a whole genotype cluster (clonal group), not of one
/// run.
///
/// Leave-one-out holds out a single run while its clone-mates stay in training, so a member of a
/// 50-run clonal cluster is nearly always "rediscovered" by markers its own clones kept alive. On
/// the census, where 3,893 runs fall into 2,666 clusters, that inflates sensitivity. Holding out
/// the whole cluster asks the question an assay user actually cares about: would this lineage's
/// markers detect an isolate from a clonal group the training set never saw?
///
/// `weights` holds, per k-mer, `needs.len()` in-group weights (the group's, then each child's, as
/// `groups::Batch`), and `removed` the weight the held-out cluster contributes to each of them.
/// `needs` are the thresholds recomputed without the cluster. A k-mer is a marker without the
/// cluster when every reduced weight still clears its reduced threshold; the score is the share of
/// those that `present` (a member of the cluster) carries. Markers and extras are pooled, because
/// an extra can rise above a lowered threshold exactly as in `loo_score`.
///
/// A cluster of one gives the same answer as `loo_score`, so only multi-member clusters need this.
#[allow(clippy::too_many_arguments)]
pub fn loco_score(
    m: &[u64],
    m_in: &[f64],
    m_removed: &[f64],
    extras: &[u64],
    extras_in: &[f64],
    extras_removed: &[f64],
    needs: &[f64],
    present: &[u64],
) -> Score {
    let stride = needs.len();
    debug_assert!(m.len() * stride == m_in.len() && m_in.len() == m_removed.len());
    debug_assert!(
        extras.len() * stride == extras_in.len() && extras_in.len() == extras_removed.len()
    );
    let (mut kept, mut held) = (0usize, 0usize);
    let mut count = |ks: &[u64], w: &[f64], rem: &[f64]| {
        for (idx, k) in ks.iter().enumerate() {
            let base = idx * stride;
            if (0..stride).all(|c| w[base + c] - rem[base + c] >= needs[c]) {
                kept += 1;
                held += present.binary_search(k).is_ok() as usize;
            }
        }
    };
    count(m, m_in, m_removed);
    count(extras, extras_in, extras_removed);
    Score {
        score: if kept > 0 {
            held as f64 / kept as f64
        } else {
            0.0
        },
        markers: kept,
        present: held,
    }
}

/// ceil(f * n) as the threshold everywhere (same float arithmetic as `groups::Target::new`).
pub fn need_for(f_in: f64, n: usize) -> u32 {
    (f_in * n as f64).ceil() as u32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loo_general_cases() {
        // n = 20 at f 0.9: need 18, need1 = ceil(17.1) = 18 -> markers X was needed for are lost
        assert_eq!((need_for(0.9, 20), need_for(0.9, 19)), (18, 18));
        let m = [1, 2, 3, 4];
        let c = [18.0, 19.0, 18.0, 20.0];
        let x = [1, 2, 9];
        let s = loo_score(&m, &c, &[], &[], &[], &x, 1.0, 18.0);
        // k-mer 1: X has it and it is in exactly 18 -> lost; 2: kept, present; 3, 4: X lacks, kept
        assert_eq!((s.present, s.markers), (1, 3));
        // n = 10 at f 0.9: need 9, need1 = ceil(8.1) = 9 too; n = 11: need 10, need1 = 9
        assert_eq!((need_for(0.9, 10), need_for(0.9, 9)), (9, 9));
        assert_eq!((need_for(0.9, 11), need_for(0.9, 10)), (10, 9));
        let c = [10.0, 10.0, 11.0, 11.0];
        let s = loo_score(&m, &c, &[10, 11], &[9.0, 9.0], &[], &[2, 10], 1.0, 9.0);
        assert_eq!((s.present, s.markers), (1, 5)); // the extras-based case, as score()
        assert_eq!(s.markers, score(&m, &[10, 11], &[2, 10], true).markers);
    }

    #[test]
    fn loo_weighted() {
        // X weighs 0.25 (a clone of 4); need 3.6, need1 = 0.9 * (4 - 0.25) = 3.375
        let m = [1, 2, 3];
        let m_in = [3.6, 3.75, 4.0];
        let ex = [7, 8];
        let ex_in = [3.4, 3.0];
        let x = [1, 2, 8];
        let s = loo_score(&m, &m_in, &ex, &ex_in, &[], &x, 0.25, 3.375);
        // 1: 3.6 - 0.25 < 3.375 lost; 2: 3.5 kept, held; 3: lacked, kept; 7: lacked, 3.4 kept;
        // 8: held by X (an extra X holds is never a marker without X)
        assert_eq!((s.present, s.markers), (1, 3));
    }

    #[test]
    fn loo_with_child_constraints() {
        // group of 10 (need 9, need1 = 9), one child of 5 holding X (need 4.5, without X 3.6)
        // stride 2: [group weight, child weight] per k-mer
        let m = [1, 2, 3];
        let m_in = [10.0, 5.0, 9.5, 4.5, 9.0, 5.0];
        let x = [1, 2, 3];
        let ch = |has_x| ChildRule {
            need: 4.5,
            need_without_x: 3.6,
            has_x,
        };
        // X holds all three. 1: 9 >= 9 and child 4 >= 3.6 -> kept. 2: 8.5 < 9 -> lost.
        // 3: 8 < 9 -> lost.
        let s = loo_score(&m, &m_in, &[], &[], &[ch(true)], &x, 1.0, 9.0);
        assert_eq!((s.present, s.markers), (1, 1));
        // the same k-mers with X in no child: the child keeps its full threshold 4.5
        // 1: group 9 ok, child 5 >= 4.5 -> kept; 2: group 8.5 < 9 lost; 3: 8 < 9 lost
        let s = loo_score(&m, &m_in, &[], &[], &[ch(false)], &x, 1.0, 9.0);
        assert_eq!((s.present, s.markers), (1, 1));
        // an extra X lacks is promoted only if it also clears the child threshold
        let ex = [20, 21];
        let ex_in = [9.0, 4.6, 9.0, 3.0];
        let s = loo_score(&[], &[], &ex, &ex_in, &[ch(false)], &[], 1.0, 9.0);
        assert_eq!(s.markers, 1, "21 fails the child threshold 4.5");
    }

    #[test]
    fn loco_matches_loo_for_a_cluster_of_one() {
        // group of 10 at f 0.9: need 9; one member of weight 1 held out -> need1 9
        let m = [1, 2, 3, 4];
        let m_in = [9.0, 10.0, 9.0, 10.0];
        let x = [1, 2];
        let a = loo_score(&m, &m_in, &[], &[], &[], &x, 1.0, 9.0);
        // the same as a cluster whose only member is X: it removes its own weight where X holds
        let removed: Vec<f64> = m.iter().map(|k| x.contains(k) as u8 as f64).collect();
        let b = loco_score(&m, &m_in, &removed, &[], &[], &[], &[9.0], &x);
        assert_eq!((a.present, a.markers), (b.present, b.markers));
        assert_eq!((b.present, b.markers), (1, 3)); // k-mer 1 falls to 8 < 9 and is lost
    }

    #[test]
    fn loco_removes_a_whole_clonal_cluster() {
        // group of 10 runs, one cluster of 4 clones sharing weight 1 (0.25 each): total 7,
        // need 6.3; without the cluster total 6, need 5.4
        let m = [1, 2, 3];
        // k-mer 1 is held only by the cluster, 2 by the cluster and others, 3 by others only
        let m_in = [6.5, 6.5, 7.0];
        let removed = [1.0, 1.0, 0.0];
        let x = [1, 2];
        let s = loco_score(&m, &m_in, &removed, &[], &[], &[], &[5.4], &x);
        // 1: 6.5 - 1 = 5.5 >= 5.4 kept and held; 2: same kept and held; 3: 7 kept, not held
        assert_eq!((s.present, s.markers), (2, 3));
        // the point of the whole exercise: a marker propped up by one clonal cluster
        // leave-one-out keeps it, because the held-out clone only takes 0.25 with it
        let one = loo_score(&[9], &[5.7], &[], &[], &[], &[9], 0.25, 5.4);
        assert_eq!(one.markers, 1, "5.7 - 0.25 = 5.45 >= 5.4");
        assert_eq!(one.score, 1.0, "the clone looks perfectly detected");
        // leave-one-cluster-out drops it: without the cluster's whole weight, 4.7 < 5.4
        let all = loco_score(&[9], &[5.7], &[1.0], &[], &[], &[], &[5.4], &[9]);
        assert_eq!((all.present, all.markers, all.score), (0, 0, 0.0));
    }

    #[test]
    fn loo_formula() {
        let m = [1, 2, 3, 4];
        let mx = [10, 11];
        let x = [2, 3, 10, 99];
        assert_eq!(n_shared(&m, &x), 2);
        let s = score(&m, &mx, &x, true);
        assert_eq!((s.present, s.markers), (2, 5)); // 4 markers + mx without the isolate's own (11)
        assert!((s.score - 0.4).abs() < 1e-12);
        let s = score(&m, &mx, &x, false);
        assert_eq!(s.markers, 4);
        assert!((s.score - 0.5).abs() < 1e-12);
    }
}
