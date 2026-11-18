//! `cladeprint run` configuration (TOML) and sample sheet (TSV), and the group hierarchy built
//! from them.
//!
//! Paths in the config are relative to the config file; paths in the sample sheet are relative to
//! the sheet. Only settings the pipeline implements are accepted (unknown keys are errors), so a
//! config never silently asks for something that is not done.

use anyhow::{Context, Result, bail, ensure};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub samples: PathBuf,
    pub output: PathBuf,
    /// sample-sheet columns that hold the group labels, top level first
    pub levels: Vec<String>,
    #[serde(default)]
    pub kmers: Kmers,
    #[serde(default)]
    pub markers: Markers,
    #[serde(default)]
    pub training: Training,
    #[serde(default)]
    pub resources: Resources,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct Training {
    /// a sample with more solid k-mers than this many times the median of its top-level group
    /// (good-coverage samples) does not train: two genomes (mixture, heterozygous hybrid) or
    /// contamination. Haploid isolates sit at <= 1.03, mixtures and hybrids at >= 1.19.
    pub max_kmer_excess: f64,
    /// a sample with fewer solid k-mers than this many times its top-level group's median does not
    /// train: a foreign organism, a partial genome (RNA-seq) or a badly degraded run. The census
    /// found 4 such runs at good coverage that the species screen had passed (excess 0.23-0.68,
    /// group median 1.00) and whose k-mers matched no training sample.
    pub min_kmer_excess: f64,
}

impl Default for Training {
    fn default() -> Self {
        Training {
            max_kmer_excess: 1.10,
            min_kmer_excess: 0.80,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct Kmers {
    pub k: u32,
    /// samples whose genome peak is below this are LOW_COVERAGE and do not train (4, as the
    /// prototype: 4x isolates still hold ~17.5M of ~18M solid k-mers)
    pub min_peak_coverage: u32,
}

impl Default for Kmers {
    fn default() -> Self {
        Kmers {
            k: 31,
            min_peak_coverage: 4,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct Markers {
    /// share of a group's training samples a marker must be in
    pub f_in: f64,
    /// minimum training samples for a group to get markers
    pub min_group_size: usize,
    /// a group with more strict candidates keeps a hash sample of about this many (and of its
    /// leave-one-out extras, with the same hash); 0 = no cap. Species-level groups otherwise reach
    /// millions of markers, all slow to filter and none needed for identification.
    pub max_candidates: usize,
    /// out-group training samples allowed to hold a marker (or an out-group k-mer within the
    /// near-miss distance): tolerates a few mislabelled or contaminated runs. A count, not a share,
    /// so a large out-group cannot hide a whole small species. 0 = the prototype's strict rule.
    pub max_out: u32,
    /// hierarchy-aware prevalence: a marker must also reach `f_in` inside every child group that
    /// was built (a species marker in >= f_in of each of its lineages). Without it a large child
    /// outvotes a small one, and a marker absent from a whole small lineage still passes.
    pub balance_children: bool,
    /// near-miss exclusion: drop a candidate within Hamming or edit distance 2 of any out-group
    /// k-mer (probes and primers bind with a few mismatches). false = the exact-specificity
    /// baseline (every strict candidate is a marker), for benchmarking only.
    pub near_miss: bool,
}

impl Default for Markers {
    fn default() -> Self {
        Markers {
            f_in: 0.9,
            min_group_size: 3,
            max_candidates: 1_000_000,
            max_out: 0,
            balance_children: true,
            near_miss: true,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct Resources {
    /// 0 = all cores
    pub threads: usize,
    /// memory for KMC counting, GB
    pub memory_gb: u32,
    /// temporary files (on disk); default <output>/tmp
    pub tmp: Option<PathBuf>,
    /// samples scored at once (step 6); memory holds one batch's present marker k-mers
    pub score_batch: usize,
}

impl Default for Resources {
    fn default() -> Self {
        Resources {
            threads: 0,
            memory_gb: 4,
            tmp: None,
            score_batch: 256,
        }
    }
}

impl Config {
    /// Reads a config; relative paths become relative to the config file's directory.
    pub fn load(path: &Path) -> Result<Config> {
        let text =
            std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
        let mut c: Config =
            toml::from_str(&text).with_context(|| format!("parsing {}", path.display()))?;
        let base = path.parent().unwrap_or(Path::new("."));
        let abs = |p: &Path| {
            if p.is_absolute() {
                p.to_path_buf()
            } else {
                base.join(p)
            }
        };
        c.samples = abs(&c.samples);
        c.output = abs(&c.output);
        c.resources.tmp = c.resources.tmp.as_deref().map(abs);
        c.validate()?;
        Ok(c)
    }

    pub fn validate(&self) -> Result<()> {
        ensure!(
            self.kmers.k == 31,
            "k = {}: only k = 31 is supported",
            self.kmers.k
        );
        ensure!(
            !self.levels.is_empty(),
            "levels: at least one sample-sheet column"
        );
        ensure!(
            self.markers.f_in > 0.0 && self.markers.f_in <= 1.0,
            "markers.f_in must be in (0, 1]"
        );
        ensure!(
            self.markers.min_group_size >= 2,
            "markers.min_group_size must be at least 2"
        );
        ensure!(
            self.training.max_kmer_excess > 1.0,
            "training.max_kmer_excess must be above 1"
        );
        ensure!(
            self.training.min_kmer_excess > 0.0 && self.training.min_kmer_excess < 1.0,
            "training.min_kmer_excess must be in (0, 1)"
        );
        let mut seen = std::collections::HashSet::new();
        for l in &self.levels {
            ensure!(seen.insert(l), "level {l} listed twice");
            ensure!(
                !RESERVED.contains(&l.as_str()),
                "level {l} uses a reserved column name"
            );
        }
        Ok(())
    }

    pub fn tmp(&self) -> PathBuf {
        self.resources
            .tmp
            .clone()
            .unwrap_or_else(|| self.output.join("tmp"))
    }
}

const RESERVED: [&str; 5] = ["sample", "reads", "assembly", "train", "cluster"];

#[derive(Debug, Clone, PartialEq)]
pub enum Input {
    /// one (single-end or interleaved) or two (paired) FASTQ files
    Reads(Vec<PathBuf>),
    Assembly(PathBuf),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Train {
    Auto,
    Yes,
    No,
}

#[derive(Debug, Clone)]
pub struct Sample {
    pub name: String,
    pub input: Input,
    /// label per level (None = unknown at that level)
    pub labels: Vec<Option<String>>,
    pub train: Train,
    /// genotype cluster (clonal group) for clade weighting; None = a cluster of its own
    pub cluster: Option<String>,
}

/// Reads the sample sheet. Columns: sample, reads and/or assembly (exactly one filled per row),
/// optional train (auto/yes/no), optional cluster (clade weighting: a group's training samples
/// that share a cluster split one sample's weight), and one column per level.
pub fn load_samples(path: &Path, levels: &[String]) -> Result<Vec<Sample>> {
    let text =
        std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let base = path.parent().unwrap_or(Path::new("."));
    let mut lines = text
        .lines()
        .enumerate()
        .filter(|(_, l)| !l.trim().is_empty() && !l.starts_with('#'));
    let (_, head) = lines.next().context("empty sample sheet")?;
    let head: Vec<&str> = head.split('\t').map(str::trim).collect();
    let col = |name: &str| head.iter().position(|&h| h == name);
    let c_sample = col("sample").context("sample sheet needs a 'sample' column")?;
    let (c_reads, c_asm, c_train) = (col("reads"), col("assembly"), col("train"));
    let c_cluster = col("cluster");
    ensure!(
        c_reads.is_some() || c_asm.is_some(),
        "sample sheet needs a 'reads' or an 'assembly' column"
    );
    let c_levels: Vec<usize> = levels
        .iter()
        .map(|l| col(l).with_context(|| format!("level column '{l}' not in the sample sheet")))
        .collect::<Result<_>>()?;

    let mut out = Vec::new();
    let mut names = std::collections::HashSet::new();
    for (ln, line) in lines {
        let f: Vec<&str> = line.split('\t').map(str::trim).collect();
        let get = |c: Option<usize>| c.and_then(|c| f.get(c)).copied().filter(|s| !s.is_empty());
        let row = ln + 1;
        let name = get(Some(c_sample))
            .with_context(|| format!("line {row}: no sample name"))?
            .to_string();
        ensure!(
            names.insert(name.clone()),
            "line {row}: sample {name} listed twice"
        );
        let path = |p: &str| {
            let p = Path::new(p);
            if p.is_absolute() {
                p.to_path_buf()
            } else {
                base.join(p)
            }
        };
        let input = match (get(c_reads), get(c_asm)) {
            (Some(r), None) => {
                let files: Vec<PathBuf> = r
                    .split(',')
                    .map(str::trim)
                    .filter(|s| !s.is_empty())
                    .map(path)
                    .collect();
                ensure!(
                    (1..=2).contains(&files.len()),
                    "line {row}: reads must be one or two FASTQ files"
                );
                Input::Reads(files)
            }
            (None, Some(a)) => Input::Assembly(path(a)),
            (Some(_), Some(_)) => {
                bail!("line {row}: {name} has both reads and an assembly; give one")
            }
            (None, None) => bail!("line {row}: {name} has neither reads nor an assembly"),
        };
        let train = match get(c_train).unwrap_or("auto").to_ascii_lowercase().as_str() {
            "auto" => Train::Auto,
            "yes" => Train::Yes,
            "no" => Train::No,
            other => bail!("line {row}: train must be auto, yes or no, not '{other}'"),
        };
        let labels = c_levels
            .iter()
            .map(|&c| get(Some(c)).map(String::from))
            .collect();
        out.push(Sample {
            name,
            input,
            labels,
            train,
            cluster: get(c_cluster).map(String::from),
        });
    }
    ensure!(!out.is_empty(), "sample sheet has no samples");
    Ok(out)
}

/// One group at one level. `path` = labels from the top level down to this one, so equal labels
/// under different parents are different groups.
#[derive(Debug, Clone, PartialEq)]
pub struct Group {
    pub level: usize,
    pub path: Vec<String>,
    /// indices of all samples carrying this label path (training or not)
    pub members: Vec<usize>,
}

impl Group {
    pub fn name(&self) -> &str {
        self.path.last().unwrap()
    }
    /// unique id, e.g. "C.neoformans/VNI/VNIa"
    pub fn id(&self) -> String {
        self.path.join("/")
    }
}

/// Groups at every level, top level first, then by name. A sample counts for a level only if its
/// labels above that level are all given (an unknown parent leaves the child unplaced). A label
/// must keep the same parent everywhere, or the sheet is inconsistent.
pub fn hierarchy(samples: &[Sample], levels: &[String]) -> Result<Vec<Group>> {
    let mut parent: HashMap<(usize, &str), Option<&str>> = HashMap::new();
    for s in samples {
        for (l, lab) in s.labels.iter().enumerate() {
            let Some(lab) = lab.as_deref() else { continue };
            let par = if l == 0 {
                None
            } else {
                s.labels[l - 1].as_deref()
            };
            if l > 0 && par.is_none() {
                bail!(
                    "sample {}: {} '{}' given without a {}",
                    s.name,
                    levels[l],
                    lab,
                    levels[l - 1]
                );
            }
            match parent.insert((l, lab), par) {
                Some(prev) if prev != par => bail!(
                    "{} '{}' sits under {} '{}' and '{}': labels must be consistent",
                    levels[l],
                    lab,
                    levels[l - 1],
                    prev.unwrap_or("-"),
                    par.unwrap_or("-")
                ),
                _ => {}
            }
        }
    }
    let mut groups: BTreeMap<(usize, Vec<String>), Vec<usize>> = BTreeMap::new();
    for (i, s) in samples.iter().enumerate() {
        let mut path = Vec::new();
        for (l, lab) in s.labels.iter().enumerate() {
            let Some(lab) = lab else { break };
            path.push(lab.clone());
            groups.entry((l, path.clone())).or_default().push(i);
        }
    }
    Ok(groups
        .into_iter()
        .map(|((level, path), members)| Group {
            level,
            path,
            members,
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
        let p = dir.join(name);
        std::fs::write(&p, text).unwrap();
        p
    }

    fn tmpdir(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("cladeprint-cfg-{tag}-{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn config_defaults_paths_and_unknown_keys() {
        let d = tmpdir("c");
        let p = write(
            &d,
            "c.toml",
            "samples = \"s.tsv\"\noutput = \"out\"\nlevels = [\"species\", \"lineage\"]\n[markers]\nf_in = 0.8\n",
        );
        let c = Config::load(&p).unwrap();
        assert_eq!(c.samples, d.join("s.tsv"));
        assert_eq!(c.tmp(), d.join("out/tmp"));
        assert_eq!(
            (c.kmers.k, c.markers.f_in, c.markers.min_group_size),
            (31, 0.8, 3)
        );
        let bad = write(
            &d,
            "b.toml",
            "samples = \"s\"\noutput = \"o\"\nlevels = [\"a\"]\n[markers]\noutgroup = \"siblings\"\n",
        );
        assert!(
            Config::load(&bad).is_err(),
            "unimplemented options must be rejected"
        );
        let k = write(
            &d,
            "k.toml",
            "samples = \"s\"\noutput = \"o\"\nlevels = [\"a\"]\n[kmers]\nk = 25\n",
        );
        assert!(Config::load(&k).is_err());
        std::fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn sheet_and_hierarchy() {
        let d = tmpdir("s");
        let lv = vec![
            "species".to_string(),
            "lineage".to_string(),
            "sub".to_string(),
        ];
        let sheet = "sample\treads\tassembly\tspecies\tlineage\tsub\ttrain\tcluster\n\
            a\tr/a_1.fq.gz, r/a_2.fq.gz\t\tCn\tVNI\tVNIa\t\tc1\n\
            b\tr/b.fq.gz\t\tCn\tVNI\tVNIb\tno\tc1\n\
            c\t\tasm/c.fa\tCn\tVNII\t\tyes\t\n\
            d\t/abs/d.fq\t\tCg\t\t\tauto\n";
        let s = load_samples(&write(&d, "s.tsv", sheet), &lv).unwrap();
        assert_eq!(
            s[0].input,
            Input::Reads(vec![d.join("r/a_1.fq.gz"), d.join("r/a_2.fq.gz")])
        );
        assert_eq!(s[2].input, Input::Assembly(d.join("asm/c.fa")));
        assert_eq!(
            s.iter().map(|x| x.cluster.as_deref()).collect::<Vec<_>>(),
            vec![Some("c1"), Some("c1"), None, None]
        );
        assert_eq!(s[3].input, Input::Reads(vec![PathBuf::from("/abs/d.fq")]));
        assert_eq!(
            (s[1].train, s[2].train, s[3].train),
            (Train::No, Train::Yes, Train::Auto)
        );
        let g = hierarchy(&s, &lv).unwrap();
        let ids: Vec<(usize, String, Vec<usize>)> = g
            .iter()
            .map(|g| (g.level, g.id(), g.members.clone()))
            .collect();
        assert_eq!(
            ids,
            vec![
                (0, "Cg".into(), vec![3]),
                (0, "Cn".into(), vec![0, 1, 2]),
                (1, "Cn/VNI".into(), vec![0, 1]),
                (1, "Cn/VNII".into(), vec![2]),
                (2, "Cn/VNI/VNIa".into(), vec![0]),
                (2, "Cn/VNI/VNIb".into(), vec![1]),
            ]
        );
        std::fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn sheet_errors() {
        let d = tmpdir("e");
        let lv = vec!["species".to_string(), "lineage".to_string()];
        let cases = [
            (
                "both inputs",
                "sample\treads\tassembly\tspecies\tlineage\na\tx.fq\ty.fa\tCn\tVNI\n",
            ),
            (
                "no input",
                "sample\treads\tassembly\tspecies\tlineage\na\t\t\tCn\tVNI\n",
            ),
            (
                "duplicate",
                "sample\treads\tspecies\tlineage\na\tx.fq\tCn\tVNI\na\ty.fq\tCn\tVNI\n",
            ),
            (
                "bad train",
                "sample\treads\tspecies\tlineage\ttrain\na\tx.fq\tCn\tVNI\tmaybe\n",
            ),
            (
                "missing level column",
                "sample\treads\tspecies\na\tx.fq\tCn\n",
            ),
        ];
        for (what, text) in cases {
            assert!(
                load_samples(&write(&d, "s.tsv", text), &lv).is_err(),
                "{what}"
            );
        }
        let inconsistent = "sample\treads\tspecies\tlineage\na\tx.fq\tCn\tVNI\nb\ty.fq\tCg\tVNI\n";
        let s = load_samples(&write(&d, "s.tsv", inconsistent), &lv).unwrap();
        assert!(hierarchy(&s, &lv).is_err(), "VNI under two species");
        let orphan = "sample\treads\tspecies\tlineage\na\tx.fq\t\tVNI\n";
        let s = load_samples(&write(&d, "s.tsv", orphan), &lv).unwrap();
        assert!(hierarchy(&s, &lv).is_err(), "lineage without species");
        std::fs::remove_dir_all(&d).unwrap();
    }
}
