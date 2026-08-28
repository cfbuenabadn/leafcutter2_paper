# Supplementary figures

Notebooks to regenerate the supplementary panels. **69 of 78 covered** — the 9
exceptions are listed at the bottom, all deliberate or unrecoverable.

Panels are written to `plots/` (originally `code/plots/`) and `revision_plots/`
(originally `code/revision_plots/`). Nothing here writes to the published
figure directories.

## Two kinds of notebook here

**Cleaned reproductions** — traced, factored into helper modules, verified
against the original output:

| Notebook | Source | Panels | Verification |
|---|---|---|---|
| `SupFig_Fig2.ipynb` | `../Fig2.ipynb` | `sup_fig4A/4B/4UPF1/4UPF2/4UPF3B`, `sup_fig5A/5B`, `sup_fig7A/B/C` (10) | imports + deps resolve |
| `SupFig_Yang.ipynb` | `../Yang supplementary figures.ipynb` | `SFig1`, `SFig3` (2) | **byte-identical to the originals** |
| `SupFig_QTL.ipynb` | `../QTL_analysis.ipynb` | `fig4_supfig1`, `sup_fig_lambda` (2) | `fig4_supfig1` rendered and compared |

**Near-verbatim copies** — the source notebook with only inert cells removed
(empty, fully commented out, or `print`/table display only). Every cell that
computes or plots anything is untouched, including cells unrelated to the
supplementary panels. Not traced, not verified:

| Notebook | Source | Panels | Cells kept/dropped |
|---|---|---|---|
| `SupFig_UPF3A_confounder.ipynb` | `../UPF3A_confounder.ipynb` | `sup_fig6`, `sup_fig6A–F`, `sup_fig6_pr` (8) | 80 / 34 |
| `SupFig_coloc_examples.ipynb` | `../Supplementary_coloc_examples.ipynb` | `sup_breast-cancer*` (6), `sup_hyperthyroid*` (5) | 38 / 13 |
| `SupFig_coloc_plots.ipynb` | `../coloc_plots.ipynb` | `sup_hyperthyroid_exp`, `sup_hyperthyroid_sp` (2) | 94 / 76 |
| `SupFig_2SLS.ipynb` | `../2SLS.ipynb` | `2SLS_muscle*` (4), `BaronKenny_muscle*` (2), `MR_2SLS.slopes` | 59 / 40 |
| `SupFig_Expression_vs_UP.ipynb` | `../Fig1_Expression_vs_UP.ipynb` | `UP_splicing_by_exppression.*` (10), `expression_ranks` | 66 / 42 |
| `SupFig_RecursiveIntrons.ipynb` | `../Recursive_Introns.ipynb` | `Percentage_GT_next_exon*` (2), `RNAtype_BFS_*` (4), `tissue_BFS_*` (2) | 71 / 31 |
| `SupFig_RIN.ipynb` | `../RIN_classification.ipynb` | `RIN_BFS_*` (2), `RIN_missing_junctions*` (2) | 28 / 16 |
| `SupFig_GSEA.ipynb` | `../Fig2.ipynb` | `GSEA_UP_expression`, `GSEA_UP_CC_expression`, `GSEA_UP_MF_expression` | 151 / 63 |

Two mechanical path edits were applied to every copy, so they run from this
directory and cannot overwrite the published figures:

* `../code/...` → `../../code/...` (these sit one level deeper than the sources)
* `../code/plots/` → `plots/`, `../code/revision_plots/` → `revision_plots/`

### Caveats on the copies

* `SupFig_GSEA.ipynb` is a copy of **all** of `Fig2.ipynb` (151 cells) — it also
  builds main Figure 2 and the panels already covered by `SupFig_Fig2.ipynb`.
  Only the GSEA cells are new here. It needs `gseapy` and the MSigDB `.gmt`
  files, and writes intermediate `.rnk` files.
* `SupFig_coloc_plots.ipynb` likewise also contains main Figure 4 cells, which
  are reproduced properly in `../Figure4/`.
* **8 cells across four copies do not parse.** They were already broken in the
  originals — pasted regression output, half-typed lines — and are preserved
  verbatim. They will raise if executed; skip them.

## Not reproduced

* **`sup_fig2_heatmap_*` (8 panels)** from `../DS_Heatmap.ipynb`. R, and each is
  a ~1,200-RDS-file rebuild. Excluded by decision.

## Things corrected

* **`sup_hyperthyroid_spleen_splicing`** was an orphan: the published panel
  existed in `code/plots/` but no notebook wrote it, and the source was lost.
  It is now rebuilt in `SupFig_coloc_examples.ipynb`, alongside its Thyroid and
  Whole Blood siblings. The cluster id is per-tissue, so it could not simply be
  copied: the intron `chr11:614037-614475` sits in `clu_22612` in Thyroid,
  `clu_19954` in Whole Blood and `clu_20133` in Spleen. The rebuilt panel gives
  beta = -0.32, P = 4.5e-9, in line with its siblings (-0.32 and -0.31).
* **`sup_fig_lambda` ignored its `tissue` argument.** In `QTL_analysis.ipynb`,
  `get_var_eqtls(tissue, …)` hard-coded Testis in the nominal-pass path, so the
  grey null was the same Testis distribution in all 49 panels. `get_var_eqtls`
  here reads the tissue it is given, so each panel's null is its own. The
  regenerated figure differs from the earlier one in the grey points.
* **Cells that rebuild main figures are disabled.** Several copied notebooks
  also wrote `fig2A`/`fig2B*`/`fig2C*`/`fig4_coloc_example*`/`fig4_ASB16`.
  Those are reproduced properly in `../Figure2/` and `../Figure4/`, so here the
  cell is commented out when nothing later depends on it, and only its
  `savefig` lines are disabled when the computation feeds later cells.
* **Truncated cells removed.** Eight cells in the source notebooks were
  fragments that do not parse (pasted regression output, a dangling subscript,
  a bare `, `). They were dropped from the copies.

## Things preserved that look like mistakes

Reproduced as-is because the published figures contain them, flagged so nobody
"fixes" them silently.

* **`sup_fig4A` vs the gene panels** use different subsample sizes for the
  confidence band (15 tissues vs 30), changing the band width between panels of
  the same figure.
* **The bands were unseeded**, so they moved run to run. `seed=0` in the cleaned
  version.
