# Figure notebooks from the companion repositories

Figure 3 and the Supplementary Note 1 benchmark were produced by Benjamin Fair
in two separate repositories. The notebooks are kept here so that every figure
in the paper can be read in one place, but **they cannot be run from this
repository** — their pipelines, configuration and inputs live in the originals.

| Notebook | Figure | Repository |
|---|---|---|
| `ComparativeSplicingFigures.qmd` | Figure 3, and the conserved-splicing supplementary panels | [bfairkun/20260825_comparativesplicing_paper][comp] |
| `NMD_GroupingDiscrepancy.qmd` | footnote: two ways of calling a testis library juvenile or adult | [comp][] |
| `SimulationBenchmarkFigures.qmd` | Supplementary Note 1, the short-read simulation benchmark | [bfairkun/20260825_leaf2simulation_paper][sim] |

The rendered versions are in [`../../docs/`](../../docs/) and are linked from the
repository README. Each is self-contained: every figure is embedded, so they open
without any of the data.

To re-run one, clone the repository it came from. Both ship the plot-ready tables
in their own `output/`, so the figures rebuild in seconds without executing the
pipelines:

    quarto render analysis/ComparativeSplicingFigures.qmd
    quarto render analysis/SimulationBenchmarkFigures.qmd

[comp]: https://github.com/bfairkun/20260825_comparativesplicing_paper
[sim]: https://github.com/bfairkun/20260825_leaf2simulation_paper
