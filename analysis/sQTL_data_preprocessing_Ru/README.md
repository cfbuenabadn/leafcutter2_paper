# ROSMAP brain sQTL preprocessing (Ru)

Ru's notebooks documenting how the ROSMAP AD data behind **Figure 5** was
preprocessed, for three brain regions: anterior cingulate (`AC/`), dorsolateral
prefrontal cortex (`DLPFC/`) and posterior cingulate (`PCC/`).

They are copied here verbatim, with their stored outputs, as a record.

## These do not run here, and are not meant to

They are **SoS** (Script of Scripts) notebooks, not Python or R. Their cells are
mostly invocations of Columbia's [xqtl-protocol][] pipeline inside Singularity
containers on Columbia hardware, e.g.

```
sos run pipeline/GWAS_QC.ipynb genotype_phenotype_sample_overlap \
    --cwd output/data_preprocessing/genotype_data/ \
    --container /mnt/vast/hpc/csg/containers_xqtl/bioinfo.sif
```

The pipeline notebooks they call (`GWAS_QC.ipynb`, `splicing_normalization.ipynb`,
`phenotype_imputation.ipynb`) live in [xqtl-protocol][] and are **not** in this
repository. Every path they reference is on `/mnt/vast/hpc/csg` or `~/Work`, on a
machine this repository has no access to. Nothing here is wired into the
Snakemake workflow, and nothing should be: they are documentation of work done
elsewhere.

[xqtl-protocol]: https://github.com/cumc/xqtl-protocol

## What they cover, and where they stop

Per region:

| Notebook | What it does |
|---|---|
| `1_*_phenotype_preprocessing` | junc files -> leafcutter2 -> splicing phenotypes; qqnorm + imputation to the association-ready bed |
| `1.2_*_leafcutter2_results_QC` | QC of the leafcutter2 output (the bulk of the stored plots and tables) |
| `2_genotype_preprocessing` | DLPFC only, and markdown only -- see below |
| `3_*_covariate_preprocessing` | genotype PCA and kinship -> the covariate file for association |

The chain ends at association-ready phenotypes and covariates. The association
analysis and plotting that produced `code/Ru_plots/*.rds` -- the direct inputs to
every Figure 5 panel -- are **not** documented here. See `code/code_report.txt`,
section 5.

## Two things that look like gaps but are not

* **AC and PCC have no `2_genotype_preprocessing`.** Genotype QC was run once for
  all regions, in the [pQTL analysis][pqtl] of the same cohort. DLPFC's copy is a
  markdown pointer to it and contains no code. Nothing is missing.
* **The `3_*_covariate` notebooks are one cell each.** Thin, but that cell holds
  the exact `sos run` invocations, which is the part worth keeping.

[pqtl]: https://github.com/cumc/fungen-xqtl-analysis/blob/main/analysis/Wang_Columbia/ROSMAP/pqtl/genotype_preprocessing.ipynb
