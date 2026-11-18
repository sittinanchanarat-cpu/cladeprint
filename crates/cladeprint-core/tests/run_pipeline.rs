//! `pipeline::run` end to end on simulated assemblies (no KMC needed): config + sample sheet ->
//! bundle. Two levels: clades A, B, C (6 genomes each) and, within A, subclades A1, A2 (3 each).
//! Planted: a 300 bp insertion in all of A, one in A1 only, one in A2 only. Every genome has
//! private SNPs. One extra A1 genome is marked train=no, and one B genome has no clade label.

use cladeprint_core::config::Config;
use cladeprint_core::kmer::{K, canonical, encode};
use cladeprint_core::testutil::Rng;
use cladeprint_core::{classify, count, kmc, pipeline};
use std::collections::{HashMap, HashSet};
use std::path::Path;

fn random_seq(rng: &mut Rng, n: usize) -> Vec<u8> {
    (0..n).map(|_| b"ACGT"[rng.below(4) as usize]).collect()
}

fn kmers(seq: &[u8]) -> HashSet<u64> {
    seq.windows(K)
        .map(|w| canonical(encode(w).unwrap()))
        .collect()
}

fn table(path: &Path) -> Vec<HashMap<String, String>> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let head: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(String::from)
        .collect();
    lines
        .map(|l| {
            head.iter()
                .cloned()
                .zip(l.split('\t').map(String::from))
                .collect()
        })
        .collect()
}

#[test]
fn run_builds_a_hierarchical_bundle() {
    let mut rng = Rng::new(7);
    let anc = random_seq(&mut rng, 20_000);
    let (ins_a, ins_a1, ins_a2) = (
        random_seq(&mut rng, 300),
        random_seq(&mut rng, 300),
        random_seq(&mut rng, 300),
    );
    let dir = std::env::temp_dir().join(format!("cladeprint-run-{}", std::process::id()));
    std::fs::create_dir_all(dir.join("asm")).unwrap();
    let mut sheet = String::from("sample\tassembly\tclade\tsub\ttrain\n");
    let mut genomes = HashMap::new();
    let plan: Vec<(String, &str, &str, &str)> = (0..3)
        .map(|i| (format!("a1_{i}"), "A", "A1", ""))
        .chain((0..3).map(|i| (format!("a2_{i}"), "A", "A2", "")))
        .chain((0..6).map(|i| (format!("b{i}"), "B", "", "")))
        .chain((0..6).map(|i| (format!("c{i}"), "C", "", "")))
        .chain([
            ("a1_extra".to_string(), "A", "A1", "no"),
            ("unlabelled".to_string(), "", "", ""),
        ])
        .collect();
    for (name, clade, sub, train) in &plan {
        let mut g = anc.clone();
        for _ in 0..15 {
            let p = rng.below(g.len() as u64) as usize;
            g[p] = if g[p] == b'A' { b'C' } else { b'A' };
        }
        let mut ins = Vec::new();
        if *clade == "A" {
            ins.extend(&ins_a);
        }
        if *sub == "A1" {
            ins.extend(&ins_a1);
        }
        if *sub == "A2" {
            ins.extend(&ins_a2);
        }
        g.splice(5_000..5_000, ins);
        std::fs::write(
            dir.join(format!("asm/{name}.fa")),
            format!(">{name}\n{}\n", String::from_utf8_lossy(&g)),
        )
        .unwrap();
        sheet += &format!("{name}\tasm/{name}.fa\t{clade}\t{sub}\t{train}\n");
        genomes.insert(name.clone(), g);
    }
    // a contaminated A1 sample: its genome plus 10 kb of foreign sequence (train=auto)
    let mut contam = genomes["a1_0"].clone();
    contam.extend(random_seq(&mut rng, 10_000));
    std::fs::write(
        dir.join("asm/a1_contam.fa"),
        format!(">c\n{}\n", String::from_utf8_lossy(&contam)),
    )
    .unwrap();
    sheet += "a1_contam\tasm/a1_contam.fa\tA\tA1\t\n";
    std::fs::write(dir.join("samples.tsv"), sheet).unwrap();
    std::fs::write(
        dir.join("c.toml"),
        "samples = \"samples.tsv\"\noutput = \"out\"\nlevels = [\"clade\", \"sub\"]\n",
    )
    .unwrap();

    let cfg = Config::load(&dir.join("c.toml")).unwrap();
    let b = pipeline::run(&cfg).unwrap();
    let ids: Vec<&str> = b.groups.iter().map(|g| g.id.as_str()).collect();
    assert_eq!(ids, ["A", "B", "C", "A/A1", "A/A2"]);
    assert!(b.skipped.is_empty());

    let bdir = cfg.output.join("bundle");
    let markers = |id: &str| -> HashSet<u64> {
        let g = b.groups.iter().find(|g| g.id == id).unwrap();
        kmc::read_sorted(&bdir.join(&g.file))
            .unwrap()
            .into_iter()
            .collect()
    };
    let (ma, ma1, ma2) = (markers("A"), markers("A/A1"), markers("A/A2"));
    assert!(kmers(&ins_a).is_subset(&ma), "clade A insertion recovered");
    assert!(
        kmers(&ins_a1).is_subset(&ma1) && kmers(&ins_a2).is_subset(&ma2),
        "subclade insertions recovered"
    );
    // out-group = everyone else: no A1 marker in an A2 genome, none in B or C
    for (name, g) in &genomes {
        let k = kmers(g);
        if !name.starts_with("a1") {
            assert!(k.is_disjoint(&ma1), "A1 marker found in {name}");
        }
        if !name.starts_with('a') {
            assert!(k.is_disjoint(&ma), "A marker found in {name}");
        }
    }
    // complete genomes: every training member carries every marker, leave-one-out included;
    // B and C differ only by private SNPs, so they have no markers and undefined sensitivity
    for g in &b.groups {
        if g.id == "B" || g.id == "C" {
            assert_eq!(g.markers, 0, "{}", g.id);
            assert!(g.loo.iter().all(|x| x.is_nan()), "{} {:?}", g.id, g.loo);
        } else {
            assert_eq!(g.loo, [1.0, 1.0, 1.0], "{} {:?}", g.id, g.loo);
        }
    }
    let qc = table(&cfg.output.join("qc.tsv"));
    let trains = |n: &str| qc.iter().find(|r| r["sample"] == n).unwrap()["train"].clone();
    assert_eq!(
        (trains("a1_extra").as_str(), trains("a1_0").as_str()),
        ("false", "true")
    );
    let row = |n: &str| qc.iter().find(|r| r["sample"] == n).unwrap().clone();
    let c = row("a1_contam");
    let x: f64 = c["kmer_excess"].parse().unwrap();
    assert_eq!(c["train"], "false", "contaminated sample must not train");
    assert!(
        c["reason"].starts_with("excess k-mers") && (1.4..1.6).contains(&x),
        "{c:?}"
    );
    for r in &qc {
        if r["sample"] != "a1_contam" {
            let x: f64 = r["kmer_excess"].parse().unwrap();
            assert!((0.97..1.03).contains(&x), "{} {x}", r["sample"]);
        }
    }
    // scores: the held-out A1 genome and the unlabelled one are scored too
    let sc = table(&cfg.output.join("scores.tsv"));
    let score = |n: &str, g: &str| {
        sc.iter().find(|r| r["sample"] == n).unwrap()[g]
            .parse::<f64>()
            .unwrap()
    };
    assert_eq!(
        (
            score("a1_extra", "A"),
            score("a1_extra", "A/A1"),
            score("a1_extra", "A/A2")
        ),
        (1.0, 1.0, 0.0)
    );
    assert_eq!((score("b0", "A"), score("unlabelled", "A")), (0.0, 0.0));
    let raw = |n: &str, g: &str| sc.iter().find(|r| r["sample"] == n).unwrap()[g].clone();
    assert_eq!(raw("b0", "B"), "NA");
    assert!(bdir.join("groups.tsv").exists() && bdir.join("bundle.toml").exists());
    // classify against the bundle: the held-out A1 genome is called A / A1 with high confidence;
    // a B genome gets no call (B has no markers, and it carries none of A's)
    let lb = classify::load(&bdir).unwrap();
    let classify_one = |name: &str| {
        let db = dir.join(format!("cls_{name}"));
        count::count_assembly(&dir.join(format!("asm/{name}.fa")), &db).unwrap();
        let present = kmc::read_present(&db, &lb.union).unwrap();
        let sc = classify::scores(&lb, &present);
        classify::call(&lb.bundle.groups, &lb.bundle.levels, &sc, false)
    };
    let c = classify_one("a1_extra");
    assert_eq!((c[0].call.as_deref(), c[0].confidence), (Some("A"), "high"));
    assert_eq!(
        (c[1].call.as_deref(), c[1].confidence),
        (Some("A/A1"), "high")
    );
    let c = classify_one("b0");
    assert_eq!((c[0].call.clone(), c[1].call.clone()), (None, None));
    assert!(c[0].flags.contains(&"NO_CALL".to_string()));
    // UNEXPLAINED: share of a sample's k-mers in no training sample (hash-sampled)
    assert!(
        lb.bundle.unexplained_loo[3] < 0.1,
        "training samples: private SNPs only, {:?}",
        lb.bundle.unexplained_loo
    );
    let novel = random_seq(&mut rng, 20_000);
    let a1 = &genomes["a1_0"];
    std::fs::write(
        dir.join("asm/hybrid.fa"),
        format!(
            ">a\n{}\n>n\n{}\n",
            String::from_utf8_lossy(a1),
            String::from_utf8_lossy(&novel)
        ),
    )
    .unwrap();
    std::fs::write(
        dir.join("asm/untrained.fa"),
        format!(
            ">n\n{}\n",
            String::from_utf8_lossy(&random_seq(&mut rng, 20_000))
        ),
    )
    .unwrap();
    let unexplained = |name: &str| {
        let db = dir.join(format!("u_{name}"));
        count::count_assembly(&dir.join(format!("asm/{name}.fa")), &db).unwrap();
        lb.read(&db).unwrap().1
    };
    let (held, hyb, unt) = (
        unexplained("a1_extra"),
        unexplained("hybrid"),
        unexplained("untrained"),
    );
    assert!(
        held < lb.bundle.unexplained_flag,
        "held-out A1 genome {held} vs flag {}",
        lb.bundle.unexplained_flag
    );
    assert!(
        (0.35..0.65).contains(&hyb) && hyb > lb.bundle.unexplained_flag,
        "hybrid {hyb}"
    );
    assert!(unt > 0.95, "untrained genome {unt}");
    // the candidate cap: rebuild with max_candidates 200 -> clade A keeps a hash sample (its
    // markers a subset of the uncapped ones), and every score stays near the uncapped one
    std::fs::write(
        dir.join("c2.toml"),
        "samples = \"samples.tsv\"\noutput = \"out2\"\nlevels = [\"clade\", \"sub\"]\n[markers]\nmax_candidates = 200\n",
    )
    .unwrap();
    std::fs::create_dir_all(dir.join("out2")).unwrap();
    for e in std::fs::read_dir(cfg.output.join("db")).unwrap() {
        let e = e.unwrap();
        std::fs::create_dir_all(dir.join("out2/db")).unwrap();
        std::fs::copy(e.path(), dir.join("out2/db").join(e.file_name())).unwrap();
    }
    let cfg2 = Config::load(&dir.join("c2.toml")).unwrap();
    let b2 = pipeline::run(&cfg2).unwrap();
    let ga2 = b2.groups.iter().find(|g| g.id == "A").unwrap();
    let ma2: HashSet<u64> = kmc::read_sorted(&cfg2.output.join("bundle").join(&ga2.file))
        .unwrap()
        .into_iter()
        .collect();
    assert!(
        ma2.is_subset(&ma) && ma2.len() < ma.len() && ma2.len() > 50,
        "{} of {}",
        ma2.len(),
        ma.len()
    );
    assert_eq!(ga2.loo, [1.0, 1.0, 1.0]);
    std::fs::remove_dir_all(&dir).unwrap();
}

/// Clade weighting and out-group tolerance through `run`. Clade A: a clonal cluster of 4 (k0-k3)
/// and 3 unrelated isolates (s0-s2); clade B: 6 isolates. Planted insertions:
///   ins_w in all of A but clone k0: 6/7 = 0.86 unweighted (fails f_in 0.9); weighted, the clones
///         share one weight: 3.75/4 = 0.94 (a marker)
///   ins_u in all of A but s2: 6/7 unweighted, 3/4 weighted (fails both ways)
///   ins_t in all of A and in b0 (a mislabelled or contaminated B run): a marker only with
///         max_out >= 1
#[test]
fn run_weights_clones_and_tolerates_one_outgroup_run() {
    let mut rng = Rng::new(11);
    let anc = random_seq(&mut rng, 20_000);
    let (ins_w, ins_u, ins_t) = (
        random_seq(&mut rng, 300),
        random_seq(&mut rng, 300),
        random_seq(&mut rng, 300),
    );
    let dir = std::env::temp_dir().join(format!("cladeprint-wt-{}", std::process::id()));
    std::fs::create_dir_all(dir.join("asm")).unwrap();
    let names: Vec<(String, &str, &str)> = (0..4)
        .map(|i| (format!("k{i}"), "A", "clone"))
        .chain((0..3).map(|i| (format!("s{i}"), "A", "")))
        .chain((0..6).map(|i| (format!("b{i}"), "B", "")))
        .collect();
    let mut plain = String::from("sample\tassembly\tclade\n");
    let mut clustered = String::from("sample\tassembly\tclade\tcluster\n");
    let clone_base = {
        let mut g = anc.clone();
        for _ in 0..30 {
            let p = rng.below(g.len() as u64) as usize;
            g[p] = if g[p] == b'A' { b'C' } else { b'A' };
        }
        g
    };
    for (name, clade, cluster) in &names {
        let mut g = if *cluster == "clone" {
            clone_base.clone()
        } else {
            anc.clone()
        };
        for _ in 0..if *cluster == "clone" { 1 } else { 15 } {
            let p = rng.below(g.len() as u64) as usize;
            g[p] = if g[p] == b'A' { b'C' } else { b'A' };
        }
        let mut ins = Vec::new();
        if *clade == "A" && name != "k0" {
            ins.extend(&ins_w);
        }
        if *clade == "A" && name != "s2" {
            ins.extend(&ins_u);
        }
        if *clade == "A" || name == "b0" {
            ins.extend(&ins_t);
        }
        g.splice(5_000..5_000, ins);
        std::fs::write(
            dir.join(format!("asm/{name}.fa")),
            format!(">{name}\n{}\n", String::from_utf8_lossy(&g)),
        )
        .unwrap();
        plain += &format!("{name}\tasm/{name}.fa\t{clade}\n");
        clustered += &format!("{name}\tasm/{name}.fa\t{clade}\t{cluster}\n");
    }
    std::fs::write(dir.join("plain.tsv"), plain).unwrap();
    std::fs::write(dir.join("clustered.tsv"), clustered).unwrap();

    let build_with = |tag: &str, sheet: &str, max_out: u32, extra: &str| {
        let c = dir.join(format!("{tag}.toml"));
        std::fs::write(
            &c,
            format!(
                "samples = \"{sheet}\"\noutput = \"{tag}\"\nlevels = [\"clade\"]\n[markers]\nmax_out = {max_out}\n{extra}"
            ),
        )
        .unwrap();
        let cfg = Config::load(&c).unwrap();
        let b = pipeline::run(&cfg).unwrap();
        let g = b.groups.iter().find(|g| g.id == "A").unwrap().clone();
        let m: HashSet<u64> = kmc::read_sorted(&cfg.output.join("bundle").join(&g.file))
            .unwrap()
            .into_iter()
            .collect();
        (g, m, cfg)
    };
    let build = |tag: &str, sheet: &str, max_out: u32| build_with(tag, sheet, max_out, "");
    let has = |m: &HashSet<u64>, ins: &[u8]| {
        let k = kmers(ins);
        (k.is_subset(m), k.is_disjoint(m))
    };

    // unweighted, strict: none of the three
    let (g, m, _) = build("plain", "plain.tsv", 0);
    assert_eq!((g.n_eff, g.need), (7.0, 7.0));
    assert_eq!(has(&m, &ins_w), (false, true), "6 of 7 is below 0.9");
    assert_eq!(has(&m, &ins_u), (false, true));
    assert_eq!(has(&m, &ins_t), (false, true), "b0 holds it");

    // clade weights + tolerance of one out-group run
    let (g, m, cfg) = build("weighted", "clustered.tsv", 1);
    assert!(
        (g.n_eff - 4.0).abs() < 1e-12 && (g.need - 3.6).abs() < 1e-6,
        "{g:?}"
    );
    assert!(has(&m, &ins_w).0, "the clones count as one isolate");
    assert_eq!(
        has(&m, &ins_u),
        (false, true),
        "an unrelated isolate still counts fully"
    );
    assert!(has(&m, &ins_t).0, "one out-group run tolerated");
    // leave-one-out: k0 lacks ins_w (a marker); without s2, ins_u is in all remaining members, so
    // it is a marker s2 lacks (about 1 in 3 of its markers); k1 holds everything
    let loo = table(&cfg.output.join("loo.tsv"));
    let sens = |n: &str| -> f64 {
        loo.iter()
            .find(|r| r["sample"] == n && r["group"] == "A")
            .unwrap()["sensitivity"]
            .parse()
            .unwrap()
    };
    assert!(sens("k0") < 1.0 && sens("k0") > 0.3, "k0 {}", sens("k0"));
    assert!((0.5..0.8).contains(&sens("s2")), "s2 {}", sens("s2"));
    assert_eq!(sens("k1"), 1.0);
    assert_eq!(loo.iter().filter(|r| r["group"] == "A").count(), 7);
    // the tolerated out-group run carries A markers: visible in its A score
    let sc = table(&cfg.output.join("scores.tsv"));
    let score = |n: &str| -> f64 {
        sc.iter().find(|r| r["sample"] == n).unwrap()["A"]
            .parse()
            .unwrap()
    };
    assert!(
        score("b0") > 0.2 && score("b1") == 0.0,
        "b0 {}",
        score("b0")
    );
    let gt = table(&cfg.output.join("bundle/groups.tsv"));
    let a = gt.iter().find(|r| r["group"] == "A").unwrap();
    assert_eq!(a["n_eff"], "4");
    // scoring in batches of one sample (each clonal cluster still kept whole) changes nothing
    let (_, _, c1) = build_with(
        "weighted_b1",
        "clustered.tsv",
        1,
        "[resources]\nscore_batch = 1\n",
    );
    for f in ["scores.tsv", "loo.tsv"] {
        let read = |c: &Config| std::fs::read_to_string(c.output.join(f)).unwrap();
        assert_eq!(read(&cfg), read(&c1), "{f} with score_batch = 1");
    }
    assert!(a["strict"].parse::<u64>().unwrap() > a["out0"].parse::<u64>().unwrap());
    std::fs::remove_dir_all(&dir).unwrap();
}

/// Unknown labels are not out-group evidence (the census case: most C. neoformans runs have no
/// lineage). Clade A: subclade A1 (3 labelled) and A2 (3); 3 more A runs carry the A1 insertion but
/// have no subclade label; 1 run with no clade label carries the A insertion. Clade B: 5 runs, and
/// b_a1 (labelled B) carries the A1 insertion too. With max_out = 1:
///   - the A1 insertion stays an A1 marker: the 3 unlabelled A runs are ignored for A1, and b_a1
///     alone is within the tolerance
///   - a second insertion in A1 and in 2 B runs is not a marker (2 out-group runs > max_out)
///   - the A insertion stays an A marker (the clade-unlabelled run is ignored for A)
#[test]
fn run_ignores_unknown_labels_for_outgroups() {
    let mut rng = Rng::new(31);
    let anc = random_seq(&mut rng, 20_000);
    let (ins_a, ins_a1, ins_x) = (
        random_seq(&mut rng, 300),
        random_seq(&mut rng, 300),
        random_seq(&mut rng, 300),
    );
    let dir = std::env::temp_dir().join(format!("cladeprint-ign-{}", std::process::id()));
    std::fs::create_dir_all(dir.join("asm")).unwrap();
    // name, clade, sub, has A insertion, has A1 insertion, has X insertion
    let mut plan: Vec<(String, &str, &str, bool, bool, bool)> = Vec::new();
    for i in 0..3 {
        plan.push((format!("a1_{i}"), "A", "A1", true, true, true));
        plan.push((format!("a2_{i}"), "A", "A2", true, false, false));
        plan.push((format!("a_nosub{i}"), "A", "", true, true, false));
    }
    plan.push(("noclade".into(), "", "", true, false, false));
    for i in 0..5 {
        plan.push((format!("b{i}"), "B", "", false, false, i < 2));
    }
    plan.push(("b_a1".into(), "B", "", false, true, false));
    let mut sheet = String::from("sample\tassembly\tclade\tsub\n");
    for (name, clade, sub, a, a1, x) in &plan {
        let mut g = anc.clone();
        for _ in 0..15 {
            let p = rng.below(g.len() as u64) as usize;
            g[p] = if g[p] == b'A' { b'C' } else { b'A' };
        }
        let mut ins = Vec::new();
        for (on, seq) in [(a, &ins_a), (a1, &ins_a1), (x, &ins_x)] {
            if *on {
                ins.extend(seq.iter());
            }
        }
        g.splice(5_000..5_000, ins);
        std::fs::write(
            dir.join(format!("asm/{name}.fa")),
            format!(">{name}\n{}\n", String::from_utf8_lossy(&g)),
        )
        .unwrap();
        sheet += &format!("{name}\tasm/{name}.fa\t{clade}\t{sub}\n");
    }
    std::fs::write(dir.join("samples.tsv"), sheet).unwrap();
    std::fs::write(
        dir.join("c.toml"),
        "samples = \"samples.tsv\"\noutput = \"out\"\nlevels = [\"clade\", \"sub\"]\n[markers]\nmax_out = 1\n",
    )
    .unwrap();
    let cfg = Config::load(&dir.join("c.toml")).unwrap();
    let b = pipeline::run(&cfg).unwrap();
    let markers = |id: &str| -> HashSet<u64> {
        let g = b.groups.iter().find(|g| g.id == id).unwrap();
        kmc::read_sorted(&cfg.output.join("bundle").join(&g.file))
            .unwrap()
            .into_iter()
            .collect()
    };
    let (ma, ma1) = (markers("A"), markers("A/A1"));
    assert!(kmers(&ins_a1).is_subset(&ma1), "A1 insertion recovered");
    assert!(
        kmers(&ins_x).is_disjoint(&ma1),
        "two B runs hold X: not a marker"
    );
    assert!(kmers(&ins_a).is_subset(&ma), "A insertion recovered");
    let gt = table(&cfg.output.join("bundle/groups.tsv"));
    let ign = |g: &str| gt.iter().find(|r| r["group"] == g).unwrap()["n_ignored"].clone();
    // A1/A2: the 3 A runs without a subclade + the clade-unlabelled run; A, B: the latter only
    assert_eq!(
        (ign("A/A1"), ign("A/A2"), ign("A"), ign("B")),
        ("4".into(), "4".into(), "1".into(), "1".into())
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

/// Hierarchy-aware prevalence (the census VNII case): clade A has a big child A1 (27 runs) and a
/// small one A2 (3 runs). An insertion carried by all of A1 but no A2 run reaches 27/30 = 0.9 of
/// the parent, so flat prevalence accepts it as an "A" marker although it detects no A2 isolate.
/// With balance_children it must be rejected, while an insertion in all of A survives both ways.
#[test]
fn run_balances_children_against_a_dominant_subgroup() {
    let mut rng = Rng::new(41);
    let anc = random_seq(&mut rng, 20_000);
    let (ins_a, ins_a1) = (random_seq(&mut rng, 300), random_seq(&mut rng, 300));
    let dir = std::env::temp_dir().join(format!("cladeprint-bal-{}", std::process::id()));
    std::fs::create_dir_all(dir.join("asm")).unwrap();
    let mut sheet = String::from("sample\tassembly\tclade\tsub\n");
    let plan: Vec<(String, &str, &str)> = (0..27)
        .map(|i| (format!("a1_{i}"), "A", "A1"))
        .chain((0..3).map(|i| (format!("a2_{i}"), "A", "A2")))
        .chain((0..6).map(|i| (format!("b{i}"), "B", "")))
        .collect();
    for (name, clade, sub) in &plan {
        let mut g = anc.clone();
        for _ in 0..15 {
            let p = rng.below(g.len() as u64) as usize;
            g[p] = if g[p] == b'A' { b'C' } else { b'A' };
        }
        let mut ins = Vec::new();
        if *clade == "A" {
            ins.extend(&ins_a);
        }
        if *sub == "A1" {
            ins.extend(&ins_a1);
        }
        g.splice(5_000..5_000, ins);
        std::fs::write(
            dir.join(format!("asm/{name}.fa")),
            format!(">{name}\n{}\n", String::from_utf8_lossy(&g)),
        )
        .unwrap();
        sheet += &format!("{name}\tasm/{name}.fa\t{clade}\t{sub}\n");
    }
    std::fs::write(dir.join("samples.tsv"), sheet).unwrap();
    let build = |tag: &str, balance: bool| {
        let c = dir.join(format!("{tag}.toml"));
        std::fs::write(
            &c,
            format!(
                "samples = \"samples.tsv\"\noutput = \"{tag}\"\nlevels = [\"clade\", \"sub\"]\n[markers]\nbalance_children = {balance}\n"
            ),
        )
        .unwrap();
        let cfg = Config::load(&c).unwrap();
        let b = pipeline::run(&cfg).unwrap();
        let g = b.groups.iter().find(|g| g.id == "A").unwrap().clone();
        let m: HashSet<u64> = kmc::read_sorted(&cfg.output.join("bundle").join(&g.file))
            .unwrap()
            .into_iter()
            .collect();
        (m, cfg)
    };
    let (flat, _) = build("flat", false);
    assert!(
        kmers(&ins_a1).is_subset(&flat),
        "flat prevalence accepts the A1-only insertion as an A marker"
    );
    assert!(kmers(&ins_a).is_subset(&flat), "the real A insertion too");

    let (bal, cfg) = build("bal", true);
    assert!(
        kmers(&ins_a1).is_disjoint(&bal),
        "balanced prevalence rejects an A marker absent from every A2 run"
    );
    assert!(
        kmers(&ins_a).is_subset(&bal),
        "the A insertion is in both children and survives"
    );
    let gt = table(&cfg.output.join("bundle/groups.tsv"));
    let row = |g: &str| gt.iter().find(|r| r["group"] == g).unwrap().clone();
    assert_eq!(row("A")["n_children"], "2");
    assert_eq!(row("A/A1")["n_children"], "0");
    // A2 members now carry the A markers they used to miss
    let sc = table(&cfg.output.join("scores.tsv"));
    let a2: Vec<f64> = sc
        .iter()
        .filter(|r| r["sample"].starts_with("a2_"))
        .map(|r| r["A"].parse().unwrap())
        .collect();
    assert!(
        a2.iter().all(|&x| x > 0.99),
        "A2 scores for A markers: {a2:?}"
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

/// Near-miss exclusion vs the exact-specificity baseline (`near_miss = false`). Every A isolate
/// carries one fixed SNP: the k-mers spanning it are exactly A-specific but one mismatch from B's
/// allele, so they are markers only in the baseline. A planted insertion is a marker both ways.
#[test]
fn exact_baseline_keeps_kmers_one_mismatch_from_the_outgroup() {
    let mut rng = Rng::new(53);
    let anc = random_seq(&mut rng, 20_000);
    let ins = random_seq(&mut rng, 300);
    let dir = std::env::temp_dir().join(format!("cladeprint-nm-{}", std::process::id()));
    std::fs::create_dir_all(dir.join("asm")).unwrap();
    let snp = 12_000;
    let mut sheet = String::from("sample\tassembly\tclade\n");
    let mut a_seq = Vec::new();
    for (name, clade) in (0..5)
        .map(|i| (format!("a{i}"), "A"))
        .chain((0..5).map(|i| (format!("b{i}"), "B")))
    {
        let mut g = anc.clone();
        for _ in 0..15 {
            let p = rng.below(g.len() as u64) as usize;
            if p.abs_diff(snp) > 40 {
                g[p] = if g[p] == b'A' { b'C' } else { b'A' };
            }
        }
        if clade == "A" {
            g[snp] = if g[snp] == b'G' { b'T' } else { b'G' };
            a_seq = g[snp - 30..=snp + 30].to_vec();
            g.splice(5_000..5_000, ins.clone());
        }
        std::fs::write(
            dir.join(format!("asm/{name}.fa")),
            format!(">{name}\n{}\n", String::from_utf8_lossy(&g)),
        )
        .unwrap();
        sheet += &format!("{name}\tasm/{name}.fa\t{clade}\n");
    }
    std::fs::write(dir.join("samples.tsv"), sheet).unwrap();
    let build = |tag: &str, near_miss: bool| {
        let c = dir.join(format!("{tag}.toml"));
        std::fs::write(
            &c,
            format!(
                "samples = \"samples.tsv\"\noutput = \"{tag}\"\nlevels = [\"clade\"]\n[markers]\nnear_miss = {near_miss}\n"
            ),
        )
        .unwrap();
        let cfg = Config::load(&c).unwrap();
        let b = pipeline::run(&cfg).unwrap();
        let g = b.groups.iter().find(|g| g.id == "A").unwrap().clone();
        let m: HashSet<u64> = kmc::read_sorted(&cfg.output.join("bundle").join(&g.file))
            .unwrap()
            .into_iter()
            .collect();
        (m, cfg)
    };
    let snp_kmers = kmers(&a_seq);
    assert_eq!(snp_kmers.len(), 31);
    let (filtered, cfg) = build("nm", true);
    let (exact, cfg_x) = build("exact", false);
    assert!(kmers(&ins).is_subset(&filtered) && kmers(&ins).is_subset(&exact));
    assert!(
        snp_kmers.is_disjoint(&filtered),
        "SNP k-mers are near misses"
    );
    assert!(snp_kmers.is_subset(&exact), "the baseline keeps them");
    assert!(filtered.is_subset(&exact) && exact.len() > filtered.len());
    let row = |cfg: &Config| {
        table(&cfg.output.join("bundle/groups.tsv"))
            .into_iter()
            .find(|r| r["group"] == "A")
            .unwrap()
    };
    let (r, rx) = (row(&cfg), row(&cfg_x));
    assert_eq!(
        rx["markers"], rx["strict"],
        "baseline: every strict candidate"
    );
    assert_eq!(r["strict"], rx["strict"]);
    std::fs::remove_dir_all(&dir).unwrap();
}

/// Leave-one-cluster-out vs leave-one-out on clonal data. Clade A is 6 clones (one cluster, so
/// one isolate's worth of weight between them) plus 3 unrelated isolates. An insertion carried by
/// all 6 clones and by exactly one unrelated isolate is a marker while the cluster is in training;
/// holding out one clone leaves it a marker (its 5 clone-mates keep it), but holding out the whole
/// cluster does not. That gap is what leave-one-out hides on clonal data.
#[test]
fn leave_one_cluster_out_is_stricter_than_leave_one_out() {
    let mut rng = Rng::new(53);
    let anc = random_seq(&mut rng, 20_000);
    let (ins_a, ins_shared) = (random_seq(&mut rng, 300), random_seq(&mut rng, 300));
    let dir = std::env::temp_dir().join(format!("cladeprint-loco-{}", std::process::id()));
    std::fs::create_dir_all(dir.join("asm")).unwrap();
    let mut sheet = String::from("sample\tassembly\tclade\tcluster\n");
    // name, clade, cluster, is a clone, carries ins_shared
    let plan: Vec<(String, &str, &str, bool, bool)> = (0..6)
        .map(|i| (format!("k{i}"), "A", "cl1", true, true))
        .chain((0..3).map(|i| (format!("a{i}"), "A", "", false, i == 0)))
        .chain((0..6).map(|i| (format!("b{i}"), "B", "", false, false)))
        .collect();
    for (name, clade, cluster, clone, shared) in &plan {
        let mut g = anc.clone();
        for _ in 0..if *clone { 1 } else { 15 } {
            let p = rng.below(g.len() as u64) as usize;
            g[p] = if g[p] == b'A' { b'C' } else { b'A' };
        }
        let mut ins = Vec::new();
        if *clade == "A" {
            ins.extend(&ins_a);
        }
        if *shared {
            ins.extend(&ins_shared);
        }
        g.splice(5_000..5_000, ins);
        std::fs::write(
            dir.join(format!("asm/{name}.fa")),
            format!(">{name}\n{}\n", String::from_utf8_lossy(&g)),
        )
        .unwrap();
        sheet += &format!("{name}\tasm/{name}.fa\t{clade}\t{cluster}\n");
    }
    std::fs::write(dir.join("samples.tsv"), sheet).unwrap();
    // f_in 0.4 over an effective n of 4 (the cluster counts once, plus 3 singletons): ins_shared
    // weighs 1 (cluster) + 1 (a0) = 2 >= 1.6, so it is a marker; without the cluster it weighs 1
    // against a threshold of 1.2, so it is not.
    std::fs::write(
        dir.join("c.toml"),
        "samples = \"samples.tsv\"\noutput = \"out\"\nlevels = [\"clade\"]\n[markers]\nf_in = 0.4\n",
    )
    .unwrap();
    let cfg = Config::load(&dir.join("c.toml")).unwrap();
    let b = pipeline::run(&cfg).unwrap();
    let ga = b.groups.iter().find(|g| g.id == "A").unwrap();
    let markers: HashSet<u64> = kmc::read_sorted(&cfg.output.join("bundle").join(&ga.file))
        .unwrap()
        .into_iter()
        .collect();
    assert!(
        kmers(&ins_a).is_subset(&markers),
        "the clade A insertion is a marker"
    );
    assert!(
        kmers(&ins_shared).is_subset(&markers),
        "so is the cluster-propped one"
    );

    let t = table(&cfg.output.join("loo.tsv"));
    let row = |n: &str| {
        t.iter()
            .find(|r| r["sample"] == n && r["group"] == "A")
            .unwrap()
            .clone()
    };
    let clone = row("k0");
    assert_eq!(
        (clone["cluster"].as_str(), clone["cluster_size"].as_str()),
        ("cl1", "6")
    );
    assert_eq!(
        row("a0")["cluster_size"],
        "1",
        "an unclustered run is its own cluster"
    );
    assert_eq!(
        row("a0")["sensitivity"],
        row("a0")["loco_sensitivity"],
        "a cluster of one gives the same answer both ways"
    );
    let n = |r: &std::collections::HashMap<String, String>, k: &str| r[k].parse::<usize>().unwrap();
    assert!(
        n(&clone, "loco_markers") < n(&clone, "markers"),
        "holding out the whole cluster leaves fewer markers: {clone:?}"
    );
    // the ~270 k-mers of ins_shared are exactly what the cluster was propping up
    assert!(
        n(&clone, "markers") - n(&clone, "loco_markers") > 200,
        "{clone:?}"
    );
    // and the clade A insertion still is a marker without the cluster, so the clone is still found
    assert!(
        clone["loco_sensitivity"].parse::<f64>().unwrap() > 0.9,
        "{clone:?}"
    );
    // group summary carries both figures
    assert!(
        ga.loco[1] <= ga.loo[1],
        "median loco {:?} loo {:?}",
        ga.loco,
        ga.loo
    );
    std::fs::remove_dir_all(&dir).unwrap();
}
