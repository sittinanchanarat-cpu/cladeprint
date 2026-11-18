# Figures of the paper

Scripts that draw every main and supplemental figure of
*cladeprint: species and lineage markers from public population genomes of the Cryptococcus species complex*.
Run them from the repository root. Outputs (PDF, PNG, TIFF) are written to `figures/msystems/out/` or, for the
scripts in `figures/`, next to the script.

| Figure | Script | Inputs |
|---|---|---|
| 1 | `figures/fig1_overview.py` → `fig1_overview.*` | `census/out/qc.tsv`, `census/out/samples.tsv` |
| 2 | `figures/make_figures.py` → `fig2_markers.*` | `census/out/scores*.tsv`, `census/out/samples.tsv`, pilot tables |
| 3 | `figures/fig3_tools.py` → `fig3_tools.*` | `census/out/bench/` |
| 4 | `figures/msystems/fig4.py` | `census/ext/nn_check.tsv`, `census/ext/summary.tsv` |
| 5 | `figures/msystems/fig5_vnii.py` | `vnii/` tables, `vnii/wm276.gff`, `figures/msystems/data/` |
| 6 | `figures/msystems/fig6_assays.py` (+ `assay_panel.py`) | `census/asm/`, `census/out/r5_kelley_mama.tsv`, `census/out/design3p/pairs_screened.tsv`, `labels/assays/` |
| S1 | `figures/figS1_sweep.py` → `figS1_sweep.*` | `census/out/{groups,loo}.tsv`, `census/out/cut2/`, `census/out/sweep/` |
| S2 | `figures/msystems/figS2_lowcov.py` | `census/ext/summary.tsv` |
| S3 | `figures/msystems/figS3_coverage.py` | `census/calib/summary.tsv` |
| S4 | `figures/fig6_vnia.py` → `fig6_vnia.*` | `census/out/loo.tsv`, `census/out/subclade/`, `census/out/r3/`, `labels/literature_labels.tsv` |
| S5 | `figures/msystems/figS5_block_tree.py` | `figures/msystems/data/block_asm.contree`, `vnii/donor/vgi_samples.tsv` |
| S6 | `figures/msystems/figS6_assays_full.py` (+ `assay_panel.py`) | `census/asm/` |

Inputs that are not tracked in git (per-run analysis tables under `census/out/` and `census/asm/`, and the
public WM276 annotation `vnii/wm276.gff`) are archived in `cladeprint_figure_inputs.zip` in the Zenodo record;
unzip it at the repository root.

The published figures use Arial. To reproduce them exactly, set `CLADEPRINT_FONTS` to a directory holding the
Arial `.ttf` files (lower-case extensions); otherwise matplotlib falls back to Liberation Sans or DejaVu Sans.
Python packages: matplotlib, numpy, pandas, biopython.
