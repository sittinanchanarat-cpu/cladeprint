"""Exact group-specific k-mers (design step 3, no near-miss exclusion yet).

For each target group (lineage) T:
  in_count(k)  = number of T isolates containing k-mer k
  out_count(k) = number of non-T isolates containing k
  candidate    = in_count >= ceil(F_IN * |T|)
  reported: candidates with out_count = 0 (strict) and <= 1 (tolerates one contaminated or
  mislabelled out-group isolate)

All set operations run inside kmc_tools (sorted on-disk databases), so memory stays flat.
Groups come from the leave-one-out lineage calls (pilot/pilot_lineages.tsv), which matched
published labels 67/67. Excluded: mixed isolates, conflicting labels, non-neoformans, and
low-coverage isolates (they would fake absences in the target group).

Usage: micromamba run -n cladeprint python group_kmers.py [F_IN]
"""
import math
import os
import subprocess
import sys
import pandas as pd

F_IN = float(sys.argv[1]) if len(sys.argv) > 1 else 0.9
os.chdir(os.path.dirname(os.path.abspath(__file__)))
os.makedirs("groups", exist_ok=True)

lin = pd.read_csv("../pilot/pilot_lineages.tsv", sep="\t").set_index("run")
stats = pd.read_csv("stats.tsv", sep="\t").set_index("run")
s = lin.join(stats, how="inner")
reason = pd.Series("", index=s.index)
reason[s.mixed_allele_frac > 0.05] = "mixed"
reason[s.role == "conflict"] = "conflict"
reason[s.peak_cov < 4] = "low_cov"  # peak 4x isolates still hold ~17.5M of ~18M solid k-mers: kept
reason[~s.good] = reason[~s.good].where(reason[~s.good] != "", "sparse")
s["use"] = reason.where(reason != "", "train")
s["major"] = s.called.map(lambda l: "VNI" if l.startswith("VNI") and not l.startswith("VNII") else l)
s[["called", "major", "use", "peak_cov", "solid_kmers", "country", "bioproject"]].to_csv("samples.tsv", sep="\t")

train = s[s.use == "train"]
print(f"{len(train)} training isolates; excluded: {s[s.use != 'train'].use.value_counts().to_dict()}")


def kmc_tools(*args):
    subprocess.run(["kmc_tools", "-t2", *args], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def count_db(name, dbs, batch=12):
    """Union of presence databases with summed counters = number of isolates holding each k-mer.

    Merges at most `batch` databases per kmc_tools call (partial sums first, then their sum):
    a single 108-way union got the session OOM-killed on this 7 GB machine.
    """
    out = f"groups/{name}"
    if os.path.exists(out + ".kmc_suf"):
        return out
    if len(dbs) > batch:
        dbs = [count_db(f"{name}.p{i // batch}", dbs[i:i + batch], batch) for i in range(0, len(dbs), batch)]
    spec = f"groups/{name}.op"
    with open(spec, "w") as f:
        f.write("INPUT:\n")
        for i, d in enumerate(dbs):
            f.write(f"r{i} = {d}\n")
        f.write("OUTPUT:\n")
        f.write(f"{out}.tmp = " + " + sum ".join(f"r{i}" for i in range(len(dbs))) + "\n")
        f.write("OUTPUT_PARAMS:\n-cs65535\n")
    kmc_tools("complex", spec)
    for ext in (".kmc_pre", ".kmc_suf"):  # rename only on success, so a killed run never looks finished
        os.replace(out + ".tmp" + ext, out + ext)
    return out


def n_info(db):
    out = subprocess.run(["kmc_tools", "info", db], capture_output=True, text=True).stdout
    for line in out.splitlines():
        if "total k-mers" in line:
            return int(line.split(":")[-1])
    return 0


targets = {"VNII": ("major", "VNII"), "VNBI": ("major", "VNBI"), "VNBII": ("major", "VNBII"),
           "VNI": ("major", "VNI"), "VNIa": ("called", "VNIa"), "VNIb": ("called", "VNIb"),
           "VNIc": ("called", "VNIc")}
dall = count_db("ALL", [f"db/{r}" for r in train.index])
rows = []
for name, (col, val) in targets.items():
    tin = list(train.index[train[col] == val])
    tout = list(train.index[train[col] != val])
    if len(tin) < 3:
        print(f"{name}: only {len(tin)} training isolates, skipped")
        continue
    need = math.ceil(F_IN * len(tin))
    din = count_db(f"{name}.in", [f"db/{r}" for r in tin])
    dout = f"groups/{name}.out"  # out-group counts = all-isolate counts minus target counts
    kmc_tools("simple", dall, din, "counters_subtract", dout)
    cand = f"groups/{name}.cand"
    kmc_tools("transform", din, "-ci" + str(need), "reduce", cand)
    strict = f"groups/{name}.strict"
    kmc_tools("simple", cand, dout, "kmers_subtract", strict)
    one = f"groups/{name}.out1"
    kmc_tools("simple", cand, dout, "-cx1", "intersect", one, "-ocleft")
    rows.append(dict(target=name, n_in=len(tin), n_out=len(tout), min_in=need,
                     candidates=n_info(cand), out0=n_info(strict), out_le1=n_info(strict) + n_info(one)))
    print(rows[-1], flush=True)

res = pd.DataFrame(rows)
res.to_csv(f"group_kmers_fin{F_IN}.tsv", sep="\t", index=False)
print(f"\nGroup-specific 31-mers (in >= {F_IN:.0%} of target, absent from out-group):")
print(res.to_string(index=False))
