# Figure S6 of the paper.
"""Fig. S6  Every published assay with a product-level check, on 274 NCBI assemblies.

The full version of Fig. 6A: all row groups of assay_panel.py (species, serotyping band
sizes, serotype alleles, serotype x mating type, molecular typing), not only the two in the main text.

Inputs as assay_panel.py. Usage (from the repo root):
  python figures/msystems/figS6_assays_full.py
"""
import importlib.util
import os

import glob as _glob
import matplotlib.font_manager as _fm
for _f in _glob.glob(os.path.join(os.environ.get("CLADEPRINT_FONTS", ""), "*.ttf")):
    _fm.fontManager.addfont(_f)

os.chdir(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # assay_panel reads census/asm/
spec = importlib.util.spec_from_file_location("r5a", "figures/msystems/assay_panel.py")
os.makedirs("figures/msystems/out", exist_ok=True)
r5a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r5a)
r5a.main(full=True, stem="figures/msystems/out/FigureS6", label=None)
