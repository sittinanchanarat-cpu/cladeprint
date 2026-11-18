# Python 2 (runs inside the Kover container). Rewrite a Kover dataset's binary phenotype for
# another target, in the dataset's own genome order. Kover's from-contigs needs a phenotype at
# creation and takes > 2 h per dataset for 204 fungal genomes, so one matrix serves every target.
# usage: python kover_swap.py <copy.kover> <target.pheno> <description>
import sys
import h5py
import numpy as np

path, pheno, desc = sys.argv[1:]
lab = dict(l.split() for l in open(pheno))
f = h5py.File(path, "r+")
ids = [g.decode() if isinstance(g, bytes) else g for g in f["genome_identifiers"][...]]
old = f["phenotype"][...]
new = np.array([int(lab[g]) for g in ids], dtype=old.dtype)
f["phenotype"][...] = new
f["phenotype"].attrs["description"] = desc
f.attrs["phenotype_description"] = desc
f.attrs["phenotype_metadata_source"] = pheno
for k in list(f.keys()):
    if k.startswith("splits"):
        del f[k]
f.close()
print("%s: %d genomes, %d positive (was %d)" % (desc, len(ids), new.sum(), old.sum()))
