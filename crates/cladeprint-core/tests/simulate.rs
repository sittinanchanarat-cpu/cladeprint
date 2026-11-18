//! End-to-end on simulated genomes with planted markers (design phase 2 acceptance test):
//! presence databases (KMC 1 layout) -> group counting -> Hamming near-miss -> edit near-miss.
//!
//! Three clades of 6 isolates share a random 30 kb ancestor; each isolate has private SNPs. Planted
//! on the branch to clade A, in all A isolates:
//!   - a 300 bp insertion of new sequence: true markers, must be recovered;
//!   - 15 SNPs: their k-mers are A-specific but 1 mismatch from the ancestral allele that the other
//!     clades carry, so the Hamming filter must remove them (the prototype's main finding: most
//!     exact group-specific k-mers are single-SNP near misses);
//!   - 15 single-base deletions: their k-mers shift the sequence, so they are far in Hamming terms
//!     but 1 edit away; the edit filter must remove the ones the Hamming filter keeps.
//!
//! Also in A: a 300 bp insertion in 5 of 6 isolates (need = 6, so in = need - 1: extras only).
//! Clade B gets its own 300 bp insertion, as a second target.

use cladeprint_core::groups::{Target, run};
use cladeprint_core::kmer::{K, canonical, encode};
use cladeprint_core::testutil::Rng;
use cladeprint_core::{edit, hamming, kmc};
use std::collections::HashSet;

fn random_seq(rng: &mut Rng, n: usize) -> Vec<u8> {
    (0..n).map(|_| b"ACGT"[rng.below(4) as usize]).collect()
}

fn kmers(seq: &[u8]) -> HashSet<u64> {
    seq.windows(K)
        .map(|w| canonical(encode(w).unwrap()))
        .collect()
}

fn snp(rng: &mut Rng, s: &mut [u8], p: usize) {
    let old = s[p];
    s[p] = loop {
        let b = b"ACGT"[rng.below(4) as usize];
        if b != old {
            break b;
        }
    };
}

#[test]
fn planted_markers_recovered_near_misses_removed() {
    let mut rng = Rng::new(2026);
    let anc = random_seq(&mut rng, 30_000);
    let ins_a = random_seq(&mut rng, 300);
    let ins_a5 = random_seq(&mut rng, 300);
    let ins_b = random_seq(&mut rng, 300);

    // clade A branch: 15 SNPs and 15 deletions, at well-separated sites (no k-mer spans two)
    let mut a = anc.clone();
    let sites: Vec<usize> = (0..30).map(|i| 1_000 + i * 900).collect();
    for &p in &sites[..15] {
        snp(&mut rng, &mut a, p);
    }
    let dels: Vec<usize> = sites[15..].to_vec();
    for &p in dels.iter().rev() {
        a.remove(p);
    }
    let a_ins_at = 500;
    let mut genomes = Vec::new();
    for clade in 0..3 {
        for i in 0..6 {
            let mut g = if clade == 0 { a.clone() } else { anc.clone() };
            for _ in 0..20 {
                let p = rng.below(g.len() as u64) as usize; // private SNPs
                snp(&mut rng, &mut g, p);
            }
            let mut g2 = g[..a_ins_at].to_vec();
            if clade == 0 {
                g2.extend(&ins_a);
                if i < 5 {
                    g2.extend(&ins_a5);
                }
            }
            if clade == 1 {
                g2.extend(&ins_b);
            }
            g2.extend(&g[a_ins_at..]);
            genomes.push(g2);
        }
    }

    let dir = std::env::temp_dir().join(format!("cladeprint-sim-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let dbs: Vec<_> = genomes
        .iter()
        .enumerate()
        .map(|(i, g)| {
            let mut v: Vec<u64> = kmers(g).into_iter().collect();
            v.sort_unstable();
            let db = dir.join(format!("iso{i}"));
            kmc::write_kmc1(&db, &v, 31, 7).unwrap();
            db
        })
        .collect();
    let targets = vec![
        Target::new("A", (0..6).collect(), 0.9),
        Target::new("B", (6..12).collect(), 0.9),
    ];
    let (mut strict, mut extra, mut out) = (vec![vec![]; 2], vec![vec![]; 2], vec![vec![]; 2]);
    run(&dbs, &targets, true, 16, |b| {
        for t in 0..2 {
            strict[t].extend(&b.strict[t]);
            extra[t].extend(&b.extra[t]);
            out[t].extend(&b.out[t]);
        }
        Ok(())
    })
    .unwrap();
    std::fs::remove_dir_all(&dir).unwrap();

    // what each planted feature contributes to clade A
    let core_a = kmers(&ins_a); // k-mers entirely inside the insertion
    let core_a5 = kmers(&ins_a5);
    let window = |s: &[u8], p: usize| kmers(&s[p.saturating_sub(K - 1)..(p + K).min(s.len())]);
    let snp_kmers: HashSet<u64> = sites[..15].iter().flat_map(|&p| window(&a, p)).collect();
    // deletion k-mers: those spanning the joined bases in a (positions shift by earlier deletions)
    let del_kmers: HashSet<u64> = dels
        .iter()
        .enumerate()
        .flat_map(|(j, &p)| window(&a, p - j))
        .filter(|k| !snp_kmers.contains(k))
        .collect();

    let s_a: HashSet<u64> = strict[0].iter().copied().collect();
    assert!(
        core_a.is_subset(&s_a),
        "insertion k-mers must all be strict candidates"
    );
    assert!(
        snp_kmers.iter().filter(|k| s_a.contains(k)).count() > 300,
        "SNP k-mers should start as candidates"
    );
    assert!(
        core_a5.is_disjoint(&s_a),
        "the 5-of-6 insertion is not a full marker"
    );
    let x_a: HashSet<u64> = extra[0].iter().copied().collect();
    assert!(
        core_a5.is_subset(&x_a),
        "the 5-of-6 insertion goes to the leave-one-out extras"
    );

    // near-miss filters, as in the pipeline: Hamming, then edit distance on the survivors
    let d1 = hamming::candidate_dist(&strict[0], &out[0]);
    let after_h: Vec<u64> = strict[0]
        .iter()
        .zip(&d1)
        .filter(|&(_, &d)| d == hamming::FAR)
        .map(|(&k, _)| k)
        .collect();
    let d2 = edit::candidate_dist(&after_h, &out[0]);
    let markers: HashSet<u64> = after_h
        .iter()
        .zip(&d2)
        .filter(|&(_, &d)| d == edit::FAR)
        .map(|(&k, _)| k)
        .collect();

    assert!(
        core_a.is_subset(&markers),
        "planted insertion markers must survive both filters"
    );
    let h_set: HashSet<u64> = after_h.iter().copied().collect();
    assert!(
        snp_kmers.is_disjoint(&h_set),
        "every SNP k-mer is a 1-mismatch near miss"
    );
    let del_after_h = del_kmers.iter().filter(|k| h_set.contains(k)).count();
    assert!(
        del_after_h > 100,
        "many deletion k-mers pass the Hamming filter ({del_after_h})"
    );
    assert!(
        del_kmers.is_disjoint(&markers),
        "the edit filter removes every deletion k-mer"
    );
    // every surviving marker comes from the A insertion (core, or junctions with >= 3 new bases)
    let a0 = &genomes[0];
    let ins_region: HashSet<u64> = kmers(&a0[a_ins_at - (K - 1)..a_ins_at + 300 + K - 1]);
    assert!(
        markers.is_subset(&ins_region),
        "{} markers outside the insertion",
        markers.difference(&ins_region).count()
    );

    // clade B: its insertion, and nothing from A's features
    let s_b: HashSet<u64> = strict[1].iter().copied().collect();
    assert!(kmers(&ins_b).is_subset(&s_b));
    assert!(core_a.is_disjoint(&s_b) && snp_kmers.is_disjoint(&s_b));
}

/// Split-k-mer clusters on simulated genomes (200 kb): two founders 1,000 SNPs from an ancestor,
/// each with clones carrying 10 private SNPs (clone pairs differ by 20 SNPs = 100 per Mb), and 3
/// unrelated isolates. One clone keeps only 60% of its k-mers (low coverage): it must stay in its
/// clone group at about the same distance, since missing k-mers remove sites, not add differences.
#[test]
fn split_kmer_clusters_find_clones_despite_low_coverage() {
    use cladeprint_core::clusters::{Options, run, single_linkage};
    let mut rng = Rng::new(21);
    let len = 200_000;
    let anc = random_seq(&mut rng, len);
    let derive = |rng: &mut Rng, base: &[u8], n: usize| {
        let mut s = base.to_vec();
        for _ in 0..n {
            let p = rng.below(len as u64) as usize;
            snp(rng, &mut s, p);
        }
        s
    };
    let f1 = derive(&mut rng, &anc, 1000);
    let f2 = derive(&mut rng, &anc, 1000);
    let mut genomes: Vec<(String, Vec<u8>)> = Vec::new();
    for i in 0..4 {
        genomes.push((format!("f1_{i}"), derive(&mut rng, &f1, 10)));
    }
    for i in 0..3 {
        genomes.push((format!("f2_{i}"), derive(&mut rng, &f2, 10)));
    }
    for i in 0..3 {
        genomes.push((format!("u{i}"), derive(&mut rng, &anc, 1000)));
    }
    let dir = std::env::temp_dir().join(format!("cladeprint-clu-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let mut dbs = Vec::new();
    for (name, g) in &genomes {
        let mut v: Vec<u64> = kmers(g).into_iter().collect();
        if name == "f1_3" {
            v.retain(|_| rng.below(10) < 6);
        }
        v.sort_unstable();
        let db = dir.join(name);
        kmc::write_kmc1(&db, &v, 31, 7).unwrap();
        dbs.push(db);
    }
    let opt = Options {
        scale: 4,
        overlap_scale: 16,
        min_positions: 50_000.0,
        exact_max: 1000.0,
        min_allele: 1,
        batches: 0,
    };
    let r = run(&dbs, &opt).unwrap();
    let d = |a: usize, b: usize| {
        r.pairs
            .iter()
            .find(|p| (p.i as usize, p.j as usize) == (a.min(b), a.max(b)))
            .map(|p| p.per_mb)
    };
    // clone pairs ~100 per Mb (loose: Poisson noise on 20 SNPs, the covered-genome estimate)
    for (a, b) in [(0, 1), (0, 2), (1, 2), (4, 5), (4, 6)] {
        let x = d(a, b).unwrap();
        assert!((60.0..140.0).contains(&x), "clones {a} {b}: {x}");
    }
    // the low-coverage clone: covers ~60% of the genome, same distance scale
    let low = &r.samples[3];
    assert!(
        (0.5..0.7).contains(&(low.positions / r.samples[0].positions)),
        "{low:?}"
    );
    for b in 0..3 {
        let x = d(b, 3).unwrap();
        assert!((50.0..160.0).contains(&x), "low-coverage clone vs {b}: {x}");
    }
    // unrelated pairs are far (~5,000-10,000 per Mb, fewer seen where SNPs crowd) and never exact
    assert!(d(0, 4).is_none() && d(7, 8).is_none() && d(0, 7).is_none());
    assert_eq!(
        single_linkage(genomes.len(), &r.pairs, 500.0),
        vec![0, 0, 0, 0, 1, 1, 1, 2, 3, 4]
    );
    // every pair is in the screen histogram; nothing near 1,000 per Mb besides the clones
    assert_eq!(r.hist.iter().map(|h| h.1).sum::<u64>(), 45);
    assert_eq!(r.too_little_overlap, 0);
    let far: u64 = r.hist.iter().filter(|h| h.0 >= 2000.0).map(|h| h.1).sum();
    assert_eq!(far, 45 - 6 - 3, "all but the clone pairs");
    assert_eq!(r.samples.iter().map(|s| s.mixed).sum::<u64>(), 0);
    std::fs::remove_dir_all(&dir).unwrap();
}
