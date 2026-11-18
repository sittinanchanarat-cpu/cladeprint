"""Correct strain names that an early classify.py took from the BioSample isolate field.

classify.py names a sample by its BioSample strain attribute, then isolate, then sample name.
classified.json was written by an earlier version that took isolate first. That only matters
where a submitter filled both fields with different values, and among census samples that is
one project: PRJNA1054030 (Chinese Academy of Sciences), whose strain field holds the strain
(YQJ105) and whose isolate field holds molecular type and a serial number ("(VGI)-11").
22 census runs carried the isolate label; this restores the strain name.

The strain column is rewritten row by row, keyed on BioSample (or run, in vgi_dist.tsv), and
every other byte of each file is left as it was. Strain names are labels only: no analysis
matches samples by them, and none of these samples was merged as a duplicate, so no result
changes. Run from the repository root; safe to run again.
"""

import json
import re

JUNK = re.compile(r"^(strain|missing|not applicable|not collected|n/?a|unknown|-|none)$", re.I)

bs = json.load(open("census/biosamples.json"))
cl = json.load(open("census/classified.json"))
fix = {}
for b, rec in cl.items():
    real = (bs.get(b, {}).get("strain") or "").strip()
    if real and not JUNK.match(real) and rec["strain"].strip() != real:
        fix[b] = (rec["strain"], real)
        rec["strain"] = real
runs = {}
for line in open("census/census_runs.tsv"):
    f = line.rstrip("\n").split("\t")
    if f[0] in fix:
        runs[f[1]] = f[0]
print(f"{len(fix)} samples, {len(runs)} census runs; projects:",
      sorted({cl[b]["proj"] for b in fix}))


def rewrite(path, key, col="strain", keymap=None):
    # newline="" keeps each file's own line endings (the manifest is CRLF)
    lines = open(path, newline="").read().split("\n")
    head = lines[0].split("\t")
    k, c = head.index(key), head.index(col)
    n = 0
    for i, line in enumerate(lines[1:], 1):
        f = line.split("\t")
        if len(f) <= max(k, c) or c == len(f) - 1:  # strain is never the last column; keeps any \r intact
            continue
        b = keymap.get(f[k]) if keymap is not None else f[k]
        if b in fix and f[c] != fix[b][1]:
            f[c] = fix[b][1]
            lines[i] = "\t".join(f)
            n += 1
    open(path, "w", newline="").write("\n".join(lines))
    print(f"{path}: {n} rows")


rewrite("census/cryptococcus_sporothrix_manifest.tsv", "biosample")
rewrite("census/census_runs.tsv", "biosample")
rewrite("vnii/donor/vgi_samples.tsv", "biosample")
rewrite("vnii/donor/vgi_dist.tsv", "run", keymap=runs)
json.dump(cl, open("census/classified.json", "w"))
