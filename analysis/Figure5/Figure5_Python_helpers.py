"""
Data loading and processing for the Python panel of Figure 5.

Created with Figure code cleaner.
Source notebook: ../Figure5.ipynb

Figures served by this module:
  * Fig5D-right       TSPAN14 eQTL / u-sQTL boxplots by genotype, with per-genotype
                      n on the x axis and the nominal effect size + P annotated

Only standard packages are imported at module level; rpy2 (needed to read the
RDS) is imported inside the function that uses it, so `from
Figure5_Python_helpers import load_plot_data` works in an environment without it.
"""

import os
import pickle

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #

BASE = '/project/yangili1/cfbuenabadn/leafcutter2_paper'
PLOTS_DIR = f'{BASE}/analysis/Figure5/plots'

TSPAN14_RDS = f'{BASE}/code/Ru_plots/leafcutter2_fig4e_right.tspan14.25Mar.rds'

# Genotype order on the x axis (Figure5.ipynb)
GENO_ORDER = ['G/G', 'G/C', 'C/C']

# Panels, left to right, with their accent colour
PANELS = [('eQTL', 'tab:blue'), ('u-sQTL', 'tab:red')]

# Genotype coding for the effect-size regression. beta and P are computed on the
# donors actually plotted (see make_Fig5D_right_data), not taken from the
# discovery-stage QTL scan.
DOSAGE = {'G/G': 0, 'G/C': 1, 'C/C': 2}   # dosage of the C allele


# R installations rpy2 can embed, in preference order. The conda environment's
# own R is tried first so the notebook works without `module load R`.
R_HOME_CANDIDATES = [
    os.path.join(os.sys.prefix, 'lib', 'R'),
    '/software/R-4.1.0-el8-x86_64/lib64/R',
]


def ensure_r_home():
    """Point rpy2 at an R installation, if the environment has not already.

    rpy2 raises 'openrlib.R_HOME cannot be None' when R_HOME is unset, which is
    the default in this conda environment. R_HOME must be set before rpy2 is
    imported, so this runs first inside `load_tspan14`.
    """
    if os.environ.get('R_HOME'):
        return os.environ['R_HOME']
    for candidate in R_HOME_CANDIDATES:
        if os.path.exists(os.path.join(candidate, 'lib', 'libR.so')):
            os.environ['R_HOME'] = candidate
            return candidate
    raise RuntimeError(
        'No R installation found for rpy2. Set R_HOME, or skip this cell and '
        'load the pickled data with load_plot_data() instead.')


def load_tspan14(rds_path=TSPAN14_RDS):
    """The TSPAN14 genotype / normalized-phenotype table, read from Ru's RDS."""
    ensure_r_home()
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    readRDS = ro.r['readRDS']
    rds = readRDS(rds_path)
    with (ro.default_converter + pandas2ri.converter).context():
        df = ro.conversion.get_conversion().rpy2py(rds)

    # rpy2 hands back 'variable' as an R factor; this collapses it to plain strings
    df['variable'] = list(df.variable)
    return df


def qtl_stats(sub, dosage=DOSAGE):
    """Nominal effect size and P for one panel, from the donors plotted.

    Linear regression of the normalized phenotype on C-allele dosage at
    chr10:80503927 (GRCh38), two-sided t test on the slope, d.f. = n - 2. These
    are nominal (unadjusted) values describing exactly the data in the panel;
    they are NOT the discovery-stage QTL statistics, which were estimated in the
    full DLPFC cohort and permutation-adjusted (see Methods / legend).
    """
    from scipy import stats

    x = sub.geno.astype(str).map(dosage).astype(float).values
    y = sub.value.astype(float).values
    res = stats.linregress(x, y)
    return dict(beta=res.slope, se=res.stderr, t=res.slope / res.stderr,
                dof=len(y) - 2, pval=res.pvalue, n=len(y))


def make_Fig5D_right_data(df, geno_order=GENO_ORDER):
    """Plot-ready pieces for Fig5D-right.

    Returns the per-panel sub-tables, the per-genotype n used for the x tick
    labels, the genotype order, and the nominal beta / P computed from the
    plotted donors. The counts come from the eQTL panel; both panels carry the
    same donors.
    """
    counts = df.loc[df.variable == 'eQTL'].geno.astype(str).value_counts()
    panels = {variable: df.loc[df.variable == variable].copy()
              for variable, _ in PANELS}
    return {
        'panels': panels,
        'counts': {g: int(counts[g]) for g in geno_order},
        'geno_order': list(geno_order),
        'stats_by_panel': {v: qtl_stats(sub) for v, sub in panels.items()},
    }


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #

PLOT_READY_VARS = ['Fig5D_right_data']


def run_all(data_dir='figure_data'):
    """Run the pipeline, pickle every plot-ready variable, and return them."""
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(PLOTS_DIR, exist_ok=True)

    Fig5D_right_data = make_Fig5D_right_data(load_tspan14())

    data = {'Fig5D_right_data': Fig5D_right_data}
    for name, value in data.items():
        with open(os.path.join(data_dir, f'{name}.pickle'), 'wb') as fh:
            pickle.dump(value, fh)
    return data


def load_plot_data(data_dir='figure_data'):
    """Load every plot-ready variable back from `data_dir` (no rpy2 needed)."""
    data = {}
    for name in PLOT_READY_VARS:
        with open(os.path.join(data_dir, f'{name}.pickle'), 'rb') as fh:
            data[name] = pickle.load(fh)
    return data
