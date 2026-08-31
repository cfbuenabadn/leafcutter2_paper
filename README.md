# LeafCutter2: detecting and quantifying unproductive splicing from short-read RNA-seq

Code and figures for the LeafCutter2 paper: annotation-free classification of splice
junctions as productive or unproductive, and what that classification reveals about
tissue-specific gene regulation, cross-species conservation, and disease risk.

## What LeafCutter2 does

LeafCutter2 keeps the original LeafCutter approach to defining alternative splicing —
pool split reads across samples, connect introns that share a donor or acceptor into a
graph, and take the connected components as clusters — so no transcript reconstruction
is needed. It then adds rarely-used junctions back into clusters and recomputes relative
usage, which matters because unproductive junctions are often cryptic and are depleted
in polyA RNA-seq after NMD has degraded them.

The new capability is **classifying which junctions produce an unproductive transcript
using only start and stop codon annotations**. A breadth-first search over combinations
of junctions looks for a path from an annotated start codon to a stop codon that avoids
a premature termination codon (PTC). A dynamic-programming formulation makes this scale:
if a PTC-free path exists from a start codon to a given 5' splice site, every PTC-free
path onward from that 5'ss inherits it. Over a million junctions are annotated in under
an hour. On top of the PTC call, LeafCutter2 implements four rules reported to modulate
how efficiently a PTC triggers NMD — the 50-nt rule, the long-exon rule, and the number
of junctions upstream and downstream of the PTC.

Applied to 17,329 GTEx samples across 49 tissues, to seven vertebrate species, and to
three ROSMAP brain regions (DLPFC n=806, AC n=603, PCC n=449), this yields thousands of unproductive splicing events that tune host
gene expression, are conserved across amniotes, and colocalize with GWAS loci —
including AD risk genes such as *TSPAN14*, *CASS4* and *PICALM*.

## Where each figure is built

The paper's figures are split across three repositories.

| Figure | Built in | What |
|---|---|---|
| Fig. 2 | [`analysis/Figure2/`](analysis/Figure2/) | unproductive splicing across 49 GTEx tissues; UPF3A correlation; splicing-vs-expression; the tissue-pattern heatmap; *GABBR1* |
| Fig. 4 | [`analysis/Figure4/`](analysis/Figure4/) | p-sQTLs and u-sQTLs across GTEx; effects on host gene expression; GWAS colocalization; the *ASB16* / bipolar disorder locus |
| Fig. 5 | [`analysis/Figure5/`](analysis/Figure5/) | ROSMAP brain u-sQTLs; enrichment against eQTLs and pQTLs; PTWAS and colocalization; *TSPAN14* |
| Supplementary | [`analysis/SupFigures/`](analysis/SupFigures/) | the supplementary panels |
| Fig. 3 | [bfairkun/20260825_comparativesplicing_paper][comp] | comparative unproductive splicing across seven vertebrate species; conserved poison exons; *ARHGAP17* |
| Supplementary Note 1 | [bfairkun/20260825_leaf2simulation_paper][sim] | the simulation benchmark against long-read ground truth, across depth, 3' coverage bias, annotation and method |

Both of Ben's repositories publish a rendered figure notebook that needs no setup:
[comparative splicing](https://bfairkun.github.io/20260825_comparativesplicing_paper/ComparativeSplicingFigures.html)
and [simulation benchmark](https://bfairkun.github.io/20260825_leaf2simulation_paper/SimulationBenchmarkFigures.html).

Figure 1 is not built in this repository.

[comp]: https://github.com/bfairkun/20260825_comparativesplicing_paper
[sim]: https://github.com/bfairkun/20260825_leaf2simulation_paper

## Rendered figure notebooks

Every panel, with the code that made it, as a single self-contained page — no setup,
nothing to install:

| | |
|---|---|
| [Figure 1](https://cfbuenabadn.github.io/leafcutter2_paper/Figure1.html) | junction classification, usage quartiles, NMD perturbations and efficiency rules |
| [Figure 2](https://cfbuenabadn.github.io/leafcutter2_paper/Figure2.html) | unproductive splicing across GTEx tissues; UPF3A; splicing vs expression |
| [Figure 2c](https://cfbuenabadn.github.io/leafcutter2_paper/Figure2_heatmap.html) | the tissue-pattern heatmap |
| [Figure 2d](https://cfbuenabadn.github.io/leafcutter2_paper/Figure2_prepare_sashimi.html) | inputs for the *GABBR1* sashimi panel |
| [Figure 4](https://cfbuenabadn.github.io/leafcutter2_paper/Figure4.html) | p-sQTLs and u-sQTLs; GWAS colocalization; the *ASB16* locus |
| [Figure 5 (Python)](https://cfbuenabadn.github.io/leafcutter2_paper/Figure5_Python.html) | *TSPAN14* mediation |
| [Figure 5 (R)](https://cfbuenabadn.github.io/leafcutter2_paper/Figure5_R.html) | QTL enrichment and PTWAS panels |

The same files are committed under [`docs/`](docs/). To rebuild them after re-running a
notebook — save it first, then:

```bash
python3 analysis/render_notebooks.py            # all of them
python3 analysis/render_notebooks.py Figure4    # just one
```

Quarto renders the outputs stored in the notebook rather than re-executing it, so this
takes seconds and the heavy cells stay run once. The script refuses to render a notebook
whose figures are not saved to disk, rather than emit a page with no plots.

## Reproducing the figures

Each figure directory follows the same three-file layout: a `*_helpers.py` (or `.R`) that
reads pipeline output and returns plot-ready objects, a `*_plot_helpers.py` that draws,
and a notebook that calls them. The notebook's first heavy cell runs `run_all()`, which
pickles everything it computed into `figure_data/`; every later cell plots from those
pickles.

`figure_data/` is committed, so the panels can be redrawn without any pipeline output:

```bash
# in analysis/Figure4/, uncomment the load_plot_data cell and skip run_all()
jupyter lab analysis/Figure4/Figure4.ipynb
```

Panels are written to each directory's `plots/`, which is deliberately not tracked —
the figures in the paper are multi-panel composites assembled from them by hand.

Figure 2c and Figure 5a-b are R; everything else is Python. They live in separate
notebooks (`Figure2_heatmap.ipynb`, `Figure5_R.ipynb`) rather than one polyglot notebook,
so each runs under its own kernel.

## Reproducing the data

```bash
cd code
mamba env create -f envs/leafcutter2_paper.yaml
conda activate leafcutter2_paper

snakemake -n                                  # dry-run; the check to run after any edit
snakemake --profile snakemake_profiles/slurm  # UChicago RCC Midway
```

Two inputs are built by tracked notebooks rather than by rules, because that is how they
were actually made: `analysis/GetGTExTables.ipynb` (the PSI tables behind Fig. 2c) and
`analysis/GetLDmatrices.ipynb` (the LD matrices behind Fig. 4f).

## What this repository cannot rebuild

Stated plainly, since a dry-run will not tell you:

* **The GTEx splicing and expression phenotypes, and the genotypes**, were produced by an
  earlier project (SpliFi) and are read here by absolute path. Those upstream rules are
  reproduced verbatim but commented out in [`code/rules/SpliFi_upstream.smk`](code/rules/SpliFi_upstream.smk),
  so the chain is documented and could be re-run, but this workflow does not re-run it.
* **`code/Ru_plots/*.rds`**, the direct inputs to every Figure 5 panel, are delivered
  pre-computed. The upstream preprocessing that produced the ROSMAP phenotypes and
  covariates is recorded in [`analysis/sQTL_data_preprocessing_Ru/`](analysis/sQTL_data_preprocessing_Ru/);
  the association analysis and plotting that produced the `.rds` themselves are not.
* **The LD matrices in `code/results/coloc/LD/`** were written by four hard-coded
  per-locus calls, so Fig. 4f reproduces for those four loci rather than for any locus.

## Layout

| Path | Contents |
|---|---|
| `analysis/Figure{2,4,5}/` | one directory per figure: notebook, helpers, and committed `figure_data/` |
| `analysis/SupFigures/` | supplementary panels |
| `analysis/sQTL_data_preprocessing_Ru/` | collaborator record: ROSMAP preprocessing, documentation only |
| `analysis/*.ipynb` | standalone analyses kept for the record (2SLS, GO, recursive introns, RIN, UPF3A) |
| `code/Snakefile`, `code/rules/` | the workflow: differential splicing and expression, e/s/u-sQTL mapping, colocalization, GWAS preparation, sashimi plots |
| `code/scripts/` | analysis scripts called by the rules |
| `code/config/` | sample lists, tissue lists, GWAS trait table, contrasts |
| `code/envs/` | conda environments |
| `code/snakemake_profiles/slurm/` | cluster profile |

`analysis/` and `code/` each ignore everything by default and re-include only what the
repository needs, so rendered figures and multi-terabyte intermediates stay out.

## Data and software

* **LeafCutter2** itself: [github.com/leafcutter2/leafcutter2](https://github.com/leafcutter2/leafcutter2).
  The published figures were generated with `bfairkun/leafcutter2` at commit `9409dd7`,
  which predates the public package's console-script interface.
* **p-sQTL and u-sQTL summary statistics** for all 49 GTEx tissues, formatted for overlap
  and statistical-genetics analyses: [Zenodo 10.5281/zenodo.15098365](https://doi.org/10.5281/zenodo.15098365).
* GTEx v8 and ROSMAP are controlled-access and are not redistributed here.

## License

MIT — see [LICENSE](LICENSE). Applies to the code in this repository; the GTEx and
ROSMAP data it reads are governed by their own access agreements.

---

This is a [workflowr][] project; `analysis/_site.yml` and `docs/` are its site scaffolding.

[workflowr]: https://github.com/jdblischak/workflowr
