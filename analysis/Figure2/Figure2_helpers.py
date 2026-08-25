"""
Data loading and processing for Figure 2.

Created with Figure code cleaner.
Source notebook: ../Fig2.ipynb

Panels served by this module (letters follow the published Figure 2 caption):
  * Fig. 2a  unproductive read percentage per sample, by tissue
             (was `fig2A` in Fig2.ipynb)
  * Fig. 2b  Spearman correlation of splicing vs expression change, per tissue pair
             (was `fig2C` in Fig2.ipynb)
  * Fig. 2e  GABBR1 expression across tissues
             (was `fig2_boxplots` in Fig2.ipynb)

Fig. 2c (heatmap) lives in Figure2_heatmap.ipynb, Fig. 2d (GABBR1 sashimi) in
Figure2_prepare_sashimi.ipynb.

Pickle keys in figure_data/ follow the same letters: fig2a_panels,
fig2a_tissue_names, fig2b_series, fig2b_source_data, fig2e_boxplot_df,
fig2e_boxplot_stats. Pickles written before the relabelling used the old names
(fig2c_series, boxplot_df, ...) and were moved to ../old_pickles/.

Only standard packages are imported at module level; heavy/optional packages
(tqdm, scipy) are imported inside the functions that use them.
"""

import os
import re
import gzip
import pickle

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Paths (verbatim from Fig2.ipynb)
# --------------------------------------------------------------------------- #

GTEX_NOISY_PHENO_DIR = '/project/yangili1/cfbuenabadn/SpliFi/code/results/pheno/noisy/GTEx/'
GTEX_TPM_TABLE = ('/project2/mstephens/cfbuenabadn/gtex-stm/code/gtex_tables/'
                  'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz')
GTEX_SAMPLE_TSV = '/project2/mstephens/cfbuenabadn/gtex-stm/data/sample.tsv'
DS_DGE_TABLE = ('/project/yangili1/cfbuenabadn/leafcutter2_paper/code/results/'
                'ds_v_dge_confounder/pairwise_comparisons.tab.gz')

# Genes pulled out of the GTEx TPM table.  The original notebook also read
# ABCA5/AKAP8L/CNNM3/DLG4/EIF4A2/MRTO4/MYOM2/NOC2L/SFT2D1/SRSF3/UPF1/UPF2/UPF3B
# for supplementary panels; none of them enter Fig. 2a, 2b or 2e.
GENES = ['UPF3A', 'GABBR1']

# Ten tissues shown in Fig. 2e (and sup_fig7*)
TEN_TISSUES = ['Brain-Anteriorcingulatecortex_BA24',
               'Brain-Cortex',
               'Brain-FrontalCortex_BA9',
               'Brain-Putamen_basalganglia',
               'Heart-AtrialAppendage',
               'Liver',
               'Lung',
               'Muscle-Skeletal',
               'Skin-NotSunExposed_Suprapubic',
               'WholeBlood']

TEN_TISSUES_CLEAN = ['Brain - Anterior\ncingulate cortex (BA24)',
                     'Brain - Cortex',
                     'Brain - Frontal\ncortex (BA9)',
                     'Brain - Putamen\n(basal ganglia)',
                     'Heart - Atrial appendage',
                     'Liver',
                     'Lung',
                     'Muscle - Skeletal',
                     'Skin - Not sun exposed\n(suprapubic)',
                     'Whole blood']

BOXPLOT_PALETTE = ['#EEEE00', '#EEEE00', '#EEEE00', '#EEEE00', '#B452CD',
                   '#CDB79E', '#9ACD32', '#7A67EE', '#3A5FCD', '#FF00FF']

# Two groups compared by the paired test annotated on Fig. 2e.
BRAIN_TISSUES = [t for t in TEN_TISSUES if t.startswith('Brain-')]
NON_BRAIN_TISSUES = [t for t in TEN_TISSUES if not t.startswith('Brain-')]

# --------------------------------------------------------------------------- #
# GTEx tissue colours and display names (verbatim from Fig2.ipynb)
# --------------------------------------------------------------------------- #

gtex_colors = {
  "Adipose-Subcutaneous": {
    "tissue_abbrv": "ADPSBQ", 
    "tissue_color_hex": "FFA54F", 
    "tissue_color_rgb": "255,165,79"
  }, 
  "Adipose-Visceral_Omentum": {
    "tissue_abbrv": "ADPVSC", 
    "tissue_color_hex": "EE9A00", 
    "tissue_color_rgb": "238,154,0"
  }, 
  "AdrenalGland": {
    "tissue_abbrv": "ADRNLG", 
    "tissue_color_hex": "8FBC8F", 
    "tissue_color_rgb": "143,188,143"
  }, 
  "Artery-Aorta": {
    "tissue_abbrv": "ARTAORT", 
    "tissue_color_hex": "8B1C62", 
    "tissue_color_rgb": "139,28,98"
  }, 
  "Artery-Coronary": {
    "tissue_abbrv": "ARTCRN", 
    "tissue_color_hex": "EE6A50", 
    "tissue_color_rgb": "238,106,80"
  }, 
  "Artery-Femoral": {
    "tissue_abbrv": "ARTFMR", 
    "tissue_color_hex": "FF4500", 
    "tissue_color_rgb": "255,69,0"
  }, 
  "Artery-Tibial": {
    "tissue_abbrv": "ARTTBL", 
    "tissue_color_hex": "FF0000", 
    "tissue_color_rgb": "255,0,0"
  }, 
  "Bladder": {
    "tissue_abbrv": "BLDDER", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Brain-Amygdala": {
    "tissue_abbrv": "BRNAMY", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain-Anteriorcingulatecortex_BA24": {
    "tissue_abbrv": "BRNACC", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain-Caudate_basalganglia": {
    "tissue_abbrv": "BRNCDT", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain-CerebellarHemisphere": {
    "tissue_abbrv": "BRNCHB", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain-Cerebellum": {
    "tissue_abbrv": "BRNCHA", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain-Cortex": {
    "tissue_abbrv": "BRNCTXA", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain-FrontalCortex_BA9": {
    "tissue_abbrv": "BRNCTXB", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain-Hippocampus": {
    "tissue_abbrv": "BRNHPP", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain-Hypothalamus": {
    "tissue_abbrv": "BRNHPT", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain-Nucleusaccumbens_basalganglia": {
    "tissue_abbrv": "BRNNCC", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain-Putamen_basalganglia": {
    "tissue_abbrv": "BRNPTM", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain-Spinalcord_cervicalc-1": {
    "tissue_abbrv": "BRNSPC", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Brain-Substantianigra": {
    "tissue_abbrv": "BRNSNG", 
    "tissue_color_hex": "EEEE00", 
    "tissue_color_rgb": "238,238,0"
  }, 
  "Breast-MammaryTissue": {
    "tissue_abbrv": "BREAST", 
    "tissue_color_hex": "00CDCD", 
    "tissue_color_rgb": "0,205,205"
  }, 
  "Cells-EBV-transformedlymphocytes": {
    "tissue_abbrv": "LCL", 
    "tissue_color_hex": "EE82EE", 
    "tissue_color_rgb": "238,130,238"
  }, 
  "Cells-Culturedfibroblasts": {
    "tissue_abbrv": "FIBRBLS", 
    "tissue_color_hex": "9AC0CD", 
    "tissue_color_rgb": "154,192,205"
  }, 
  "Cervix-Ectocervix": {
    "tissue_abbrv": "CVXECT", 
    "tissue_color_hex": "EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Cervix-Endocervix": {
    "tissue_abbrv": "CVSEND", 
    "tissue_color_hex": "EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Colon-Sigmoid": {
    "tissue_abbrv": "CLNSGM", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Colon-Transverse": {
    "tissue_abbrv": "CLNTRN", 
    "tissue_color_hex": "EEC591", 
    "tissue_color_rgb": "238,197,145"
  }, 
  "Esophagus-GastroesophagealJunction": {
    "tissue_abbrv": "ESPGEJ", 
    "tissue_color_hex": "8B7355", 
    "tissue_color_rgb": "139,115,85"
  }, 
  "Esophagus-Mucosa": {
    "tissue_abbrv": "ESPMCS", 
    "tissue_color_hex": "8B7355", 
    "tissue_color_rgb": "139,115,85"
  }, 
  "Esophagus-Muscularis": {
    "tissue_abbrv": "ESPMSL", 
    "tissue_color_hex": "CDAA7D", 
    "tissue_color_rgb": "205,170,125"
  }, 
  "FallopianTube": {
    "tissue_abbrv": "FLLPNT", 
    "tissue_color_hex": "EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Heart-AtrialAppendage": {
    "tissue_abbrv": "HRTAA", 
    "tissue_color_hex": "B452CD", 
    "tissue_color_rgb": "180,82,205"
  }, 
  "Heart-LeftVentricle": {
    "tissue_abbrv": "HRTLV", 
    "tissue_color_hex": "7A378B", 
    "tissue_color_rgb": "122,55,139"
  }, 
  "Kidney-Cortex": {
    "tissue_abbrv": "KDNCTX", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Kidney-Medulla": {
    "tissue_abbrv": "KDNMDL", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Liver": {
    "tissue_abbrv": "LIVER", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Lung": {
    "tissue_abbrv": "LUNG", 
    "tissue_color_hex": "9ACD32", 
    "tissue_color_rgb": "154,205,50"
  }, 
  "MinorSalivaryGland": {
    "tissue_abbrv": "SLVRYG", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Muscle-Skeletal": {
    "tissue_abbrv": "MSCLSK", 
    "tissue_color_hex": "7A67EE", 
    "tissue_color_rgb": "122,103,238"
  }, 
  "Nerve-Tibial": {
    "tissue_abbrv": "NERVET", 
    "tissue_color_hex": "FFD700", 
    "tissue_color_rgb": "255,215,0"
  }, 
  "Ovary": {
    "tissue_abbrv": "OVARY", 
    "tissue_color_hex": "FFB6C1", 
    "tissue_color_rgb": "255,182,193"
  }, 
  "Pancreas": {
    "tissue_abbrv": "PNCREAS", 
    "tissue_color_hex": "CD9B1D", 
    "tissue_color_rgb": "205,155,29"
  }, 
  "Pituitary": {
    "tissue_abbrv": "PTTARY", 
    "tissue_color_hex": "B4EEB4", 
    "tissue_color_rgb": "180,238,180"
  }, 
  "Prostate": {
    "tissue_abbrv": "PRSTTE", 
    "tissue_color_hex": "D9D9D9", 
    "tissue_color_rgb": "217,217,217"
  }, 
  "Skin-NotSunExposed_Suprapubic": {
    "tissue_abbrv": "SKINNS", 
    "tissue_color_hex": "3A5FCD", 
    "tissue_color_rgb": "58,95,205"
  }, 
  "Skin-SunExposed_Lowerleg": {
    "tissue_abbrv": "SKINS", 
    "tissue_color_hex": "1E90FF", 
    "tissue_color_rgb": "30,144,255"
  }, 
  "SmallIntestine-TerminalIleum": {
    "tissue_abbrv": "SNTTRM", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Spleen": {
    "tissue_abbrv": "SPLEEN", 
    "tissue_color_hex": "CDB79E", 
    "tissue_color_rgb": "205,183,158"
  }, 
  "Stomach": {
    "tissue_abbrv": "STMACH", 
    "tissue_color_hex": "FFD39B", 
    "tissue_color_rgb": "255,211,155"
  }, 
  "Testis": {
    "tissue_abbrv": "TESTIS", 
    "tissue_color_hex": "A6A6A6", 
    "tissue_color_rgb": "166,166,166"
  }, 
  "Thyroid": {
    "tissue_abbrv": "THYROID", 
    "tissue_color_hex": "008B45", 
    "tissue_color_rgb": "0,139,69"
  }, 
  "Uterus": {
    "tissue_abbrv": "UTERUS", 
    "tissue_color_hex": "EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "Vagina": {
    "tissue_abbrv": "VAGINA", 
    "tissue_color_hex": "EED5D2", 
    "tissue_color_rgb": "238,213,210"
  }, 
  "WholeBlood": {
    "tissue_abbrv": "WHLBLD", 
    "tissue_color_hex": "FF00FF", 
    "tissue_color_rgb": "255,0,255"
  }
}

sorted_tissues_names = ['Testis', 'Brain - Cerebellum', 'Brain - Cerebellar Hemisphere', 'Thyroid',
       'Prostate', 'Spleen', 'Nerve - Tibial', 'Small Intestine - Terminal Ileum',
       'Pituitary', 'Ovary', 'Lung', 'Uterus', 'Kidney - Cortex',
       'Breast - Mammary Tissue', 'Vagina', 'Skin - Not Sun Exposed (Suprapubic)',
       'Adipose - Subcutaneous', 'Skin - Sun Exposed (Lower leg)', 'Colon - Sigmoid',
       'Brain - Spinal cord (cervical c-1)', 'Whole Blood', 'Brain - Cortex',
       'Esophagus - Gastroesophageal Junction', 'Brain - Hippocampus',
       'Colon - Transverse', 'Adipose - Visceral (Omentum)', 'Minor Salivary Gland',
       'Brain - Hypothalamus', 'Esophagus - Muscularis', 'Brain - Substantianigra',
       'Liver', 'Brain - Caudate (basal ganglia)',
       'Brain - Nucleus accumbens (basal ganglia)', 'Stomach',
       'Brain - Putamen (basal ganglia)', 'Artery - Tibial', 'Brain - Amygdala',
       'Brain - Frontal Cortex (BA9)', 'Brain - Anterior cingulate cortex (BA24)',
       'Artery - Coronary', 'Artery - Aorta', 'Heart - Atrial Appendage',
       'Adrenal Gland', 'Esophagus - Mucosa', 'Cells - EBV-transformed lymphocytes',
       'Pancreas', 'Heart - Left Ventricle', 'Muscle - Skeletal',
       'Cells - Cultured fibroblasts']

# --------------------------------------------------------------------------- #
# Loading / processing
# --------------------------------------------------------------------------- #

def transform_string(s):
    """GTEx tissue-name normaliser used throughout Fig2.ipynb."""
    s = s.replace(" ", "")        # Remove spaces
    s = re.sub(r"\(", "_", s)     # Replace '(' with '_'
    s = re.sub(r"\)", "", s)      # Remove ')'
    return s


def load_unproductive_pct(pheno_dir=GTEX_NOISY_PHENO_DIR):
    """Per-sample percentage of unproductive (UP) junction reads, per tissue.

    Reads the leafcutter2 `noise_by_intron` count tables for every GTEx tissue
    and returns `pct_df` indexed by `<donor>.<tissue>`.
    """
    from tqdm import tqdm

    tissues = os.listdir(pheno_dir)

    tissues_list = []
    pct_list = []
    total_counts_list = []
    up_counts_list = []
    cols_list = []

    for tissue in tissues:
        all_counts_array = 0
        up_counts_array = 0

        counts_file = f'{pheno_dir}{tissue}/leafcutter_perind_numers.counts.noise_by_intron.gz'
        with gzip.open(counts_file, 'rb') as fh:
            cols = [x.split('.')[0] + '.' + tissue
                    for x in fh.readline().decode().rstrip().split(' ')[1:]]
            cols_list += cols
            for x in tqdm(fh):
                sample_row = x.decode().rstrip().split(' ')
                itype = sample_row[0].split(':')[-1]
                sample_array = np.array([int(x) for x in sample_row[1:]])
                if itype == 'UP':
                    up_counts_array += sample_array
                all_counts_array += sample_array

        total_counts_list += list(all_counts_array)
        up_counts_list += list(up_counts_array)
        pct = 100 * (up_counts_array / all_counts_array)
        tissues_list += [tissue] * len(pct)
        pct_list += list(pct)

    pct_df = pd.DataFrame()
    pct_df['tissue'] = tissues_list
    pct_df['pct'] = pct_list
    pct_df['total_counts'] = total_counts_list
    pct_df['up_counts'] = up_counts_list
    pct_df.index = cols_list

    pct_df['log_counts'] = np.log10(pct_df.total_counts)
    pct_df['alt_id'] = pct_df.index

    return pct_df


def load_gene_tpm(genes=GENES, gtex_table=GTEX_TPM_TABLE, sample_tsv=GTEX_SAMPLE_TSV):
    """Per-sample TPM for `genes`, annotated with the GTEx tissue of each sample."""
    from tqdm import tqdm

    samples = pd.read_csv(sample_tsv, sep='\t')
    samples['tissue'] = samples.tissue_site_detail.apply(lambda x: transform_string(x))

    genes = list(genes)
    tpm = {}
    with gzip.open(gtex_table, 'rb') as fh:
        fh.readline()
        fh.readline()
        sample_list = fh.readline().decode().rstrip().split('\t')[2:]
        for line in tqdm(fh):
            row = line.decode().rstrip().split('\t')
            if row[1] in genes:
                tpm[row[1]] = [float(x) for x in row[2:]]

    gene_df = pd.DataFrame()
    for gene in genes:
        gene_df[gene] = tpm[gene]
    gene_df.index = sample_list

    samples['alt_id'] = (samples['entity:sample_id'].apply(lambda x: '-'.join(x.split('-')[:2]))
                         + '.' + samples.tissue)
    gene_df = gene_df.merge(samples[['entity:sample_id', 'tissue', 'alt_id']],
                            left_index=True, right_on='entity:sample_id')

    return gene_df


def build_upf_df(gene_df, pct_df):
    """Join per-sample gene TPM with the per-sample unproductive-read percentage."""
    return gene_df.merge(pct_df, left_on=['alt_id', 'tissue'], right_on=['alt_id', 'tissue'])


def make_fig2a_panels(upf_df, pct_df):
    """Plot-ready per-tissue series for fig2A.

    Tissues are ordered by decreasing median unproductive-read percentage and
    then reversed, exactly as in `sorted_tissues[::-1]` in Fig2.ipynb, so that
    the returned list is already in left-to-right plotting order.
    """
    sorted_tissues = pct_df.groupby('tissue').pct.median().sort_values().index[::-1]

    panels = []
    for tissue in sorted_tissues[::-1]:
        sub = upf_df.loc[upf_df.tissue == tissue, ]
        panels.append({
            'tissue': tissue,
            'color': '#' + gtex_colors[tissue]['tissue_color_hex'],
            'pct': np.array(sub.pct.sort_values()),
            'upf3a': np.array(sub.sort_values('pct').UPF3A),
        })

    tissue_names = list(sorted_tissues_names[::-1])

    return panels, tissue_names


def make_boxplot_df(upf_df, gene='GABBR1', tissues=TEN_TISSUES):
    """Plot-ready long dataframe for Fig. 2e (tissue x gene TPM).

    One row = one GTEx RNA-seq sample = one tissue of one post-mortem donor,
    i.e. the unit of study. `donor` is carried along so that sample sizes can be
    reported as both samples and independent donors.
    """
    sub = upf_df.loc[upf_df.tissue.isin(tissues)].sort_values('tissue').copy()
    sub['donor'] = sub['entity:sample_id'].apply(lambda x: '-'.join(x.split('-')[:2]))
    return sub[['tissue', 'donor', gene]]


def make_boxplot_stats(fig2e_boxplot_df, gene='GABBR1', tissues=TEN_TISSUES,
                       group_a=BRAIN_TISSUES, group_b=NON_BRAIN_TISSUES,
                       primary='wilcoxon'):
    """Sample sizes and the brain-vs-non-brain group test for Fig. 2e.

    Each box is one tissue, and within a tissue every GTEx donor contributes
    exactly one RNA-seq sample, so the per-box n is already a count of
    independent individuals (`n_per_tissue` should equal `donors_per_tissue`;
    both are returned so this can be checked). No technical replicates are
    involved: distinct brain regions from one donor are distinct biological
    samples, not repeated measurements of one sample.

    The group comparison pools four brain tissues against six non-brain tissues,
    and a donor usually contributes to both sides. Each donor is therefore first
    collapsed to one value per group -- the median across that donor's tissues
    within the group -- and the two groups are compared with:

      * `wilcoxon`     : Wilcoxon signed-rank, two-sided, over the donors present
                         in BOTH groups. Each donor serves as their own control,
                         so the groups share no individuals and n is an exact
                         count of paired individuals. This is the default.
      * `mannwhitney`  : Mann-Whitney U, two-sided, over the donor-level values
                         of each group. Retained for comparison; note the two
                         groups are not disjoint, since `donors_in_both` donors
                         appear on both sides.

    `primary` names which of the two the figure should annotate. Both are always
    computed. Each sub-dict carries `label`, `test`, `stat_name`, `statistic` and
    `pvalue` so the plotting code stays agnostic about which one is shown.
    """
    from scipy.stats import mannwhitneyu, wilcoxon

    n_per_tissue = [int((fig2e_boxplot_df.tissue == t).sum()) for t in tissues]
    donors_per_tissue = [int(fig2e_boxplot_df.loc[fig2e_boxplot_df.tissue == t, 'donor'].nunique())
                         for t in tissues]

    a = fig2e_boxplot_df.loc[fig2e_boxplot_df.tissue.isin(group_a)]
    b = fig2e_boxplot_df.loc[fig2e_boxplot_df.tissue.isin(group_b)]

    # One value per donor per group: the median across that donor's tissues.
    a_donor = a.groupby('donor')[gene].median()
    b_donor = b.groupby('donor')[gene].median()
    shared = a_donor.index.intersection(b_donor.index)

    U, mw_pvalue = mannwhitneyu(a_donor, b_donor, alternative='two-sided')
    mannwhitney = {
        'label': ('Brain (n = {} donors) v non-brain (n = {} donors)'
                  .format(len(a_donor), len(b_donor))),
        'test': 'Mann-Whitney U, two-sided',
        'stat_name': 'U',
        'statistic': float(U),
        'pvalue': float(mw_pvalue),
        'n_a': int(len(a_donor)),
        'n_b': int(len(b_donor)),
    }

    wilcoxon_result = None
    if len(shared) > 0:
        paired_a = a_donor.loc[shared]
        paired_b = b_donor.loc[shared]
        W, w_pvalue = wilcoxon(paired_a, paired_b, alternative='two-sided')
        wilcoxon_result = {
            'label': 'Brain v non-brain, n = {} paired donors'.format(len(shared)),
            'test': 'Wilcoxon signed-rank, two-sided',
            'stat_name': 'W',
            'statistic': float(W),
            'pvalue': float(w_pvalue),
            'n_pairs': int(len(shared)),
            # Descriptive effect size: median of the within-donor differences.
            'median_difference': float(np.median(paired_a - paired_b)),
        }

    if primary == 'wilcoxon' and wilcoxon_result is None:
        primary = 'mannwhitney'

    return {
        'gene': gene,
        'tissues': list(tissues),
        'n_per_tissue': n_per_tissue,
        'donors_per_tissue': donors_per_tissue,
        'group_a': list(group_a),
        'group_b': list(group_b),
        'n_samples_a': int(a.shape[0]),
        'n_samples_b': int(b.shape[0]),
        'n_donors_a': int(a_donor.shape[0]),
        'n_donors_b': int(b_donor.shape[0]),
        'donors_in_both': int(len(shared)),
        'wilcoxon': wilcoxon_result,
        'mannwhitney': mannwhitney,
        'primary': primary,
    }


def load_ds_dge(ds_dge_table=DS_DGE_TABLE):
    """Pairwise-tissue splicing-vs-expression correlations, Bladder dropped, FDR added."""
    from scipy.stats import false_discovery_control

    ds_dge = pd.read_csv(ds_dge_table, sep='\t', index_col=0)
    ds_dge = ds_dge.loc[[x for x in ds_dge.index if "Bladder" not in x]]
    ds_dge['FDR'] = false_discovery_control(ds_dge.spearman_pval)
    ds_dge['FDR_p'] = false_discovery_control(ds_dge.spearman_pval_p)
    return ds_dge


def make_fig2b_series(ds_dge, min_n=50, fdr=0.1):
    """Plot-ready (x, y) scatter series for the two Fig. 2b volcano panels."""
    def _series(rho_col, pval_col, fdr_col, significant):
        keep = ds_dge.n >= min_n
        keep &= (ds_dge[fdr_col] <= fdr) if significant else (ds_dge[fdr_col] > fdr)
        sub = ds_dge.loc[keep].sort_values(pval_col)
        return {'x': np.array(sub[rho_col][::-1]),
                'y': np.array(-np.log10(sub[pval_col])[::-1])}

    return {
        'unproductive_ns': _series('spearman', 'spearman_pval', 'FDR', False),
        'unproductive_sig': _series('spearman', 'spearman_pval', 'FDR', True),
        'productive_ns': _series('spearman_p', 'spearman_pval_p', 'FDR_p', False),
        'productive_sig': _series('spearman_p', 'spearman_pval_p', 'FDR_p', True),
    }


def spearman_ci(rho, n, alpha=0.05):
    """Two-sided confidence interval for Spearman's rho.

    Fisher z transform with the Bonett-Wright standard error,
    SE_z = sqrt((1 + rho^2 / 2) / (n - 3)), which is the standard
    correction for rank correlations. Returns (low, high) arrays;
    entries with n <= 3 are NaN.
    """
    from scipy.stats import norm

    rho = np.asarray(rho, dtype=float)
    n = np.asarray(n, dtype=float)

    with np.errstate(invalid='ignore', divide='ignore'):
        se = np.sqrt((1 + rho ** 2 / 2) / (n - 3))
        z = np.arctanh(rho)
        crit = norm.ppf(1 - alpha / 2)
        low = np.tanh(z - crit * se)
        high = np.tanh(z + crit * se)

    valid = n > 3
    return np.where(valid, low, np.nan), np.where(valid, high, np.nan)


def make_fig2b_source_data(ds_dge, min_n=50, fdr=0.1, alpha=0.05):
    """Exact per-point statistics behind Fig. 2b, for the Source Data table.

    One row per plotted point: the tissue pair, which splicing class it belongs
    to, the exact n, degrees of freedom, Spearman's rho, its confidence
    interval, the exact two-sided p-value and the Benjamini-Hochberg adjusted
    p-value that determines the point's colour.

    Note on the adjustment: `FDR` / `FDR_p` are carried over from `load_ds_dge`,
    where Benjamini-Hochberg is applied across *every* tissue pair remaining
    after Bladder removal -- including pairs with n < min_n that this table (and
    the figure) does not show. That is the adjustment actually used for the
    published panel, so it is reported as-is.
    """
    frames = []
    for label, rho_col, pval_col, fdr_col, n_col in [
            ('unproductive', 'spearman', 'spearman_pval', 'FDR', 'n'),
            ('productive', 'spearman_p', 'spearman_pval_p', 'FDR_p', 'n_p')]:

        sub = ds_dge.loc[ds_dge.n >= min_n].sort_values(pval_col)
        ci_low, ci_high = spearman_ci(sub[rho_col], sub[n_col], alpha=alpha)

        frames.append(pd.DataFrame({
            'tissue_pair': sub.index,
            'splicing_class': label,
            'n_genes': sub[n_col].to_numpy(),
            'df': sub[n_col].to_numpy() - 2,
            'spearman_rho': sub[rho_col].to_numpy(),
            f'ci_low_{int((1 - alpha) * 100)}': ci_low,
            f'ci_high_{int((1 - alpha) * 100)}': ci_high,
            'p_value_two_sided': sub[pval_col].to_numpy(),
            'bh_adjusted_p': sub[fdr_col].to_numpy(),
            'significant_at_fdr_{}'.format(fdr): sub[fdr_col].to_numpy() <= fdr,
        }))

    return pd.concat(frames, ignore_index=True)


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #

PLOT_READY_VARS = ['fig2a_panels', 'fig2a_tissue_names', 'fig2e_boxplot_df', 'fig2e_boxplot_stats',
                   'fig2b_series', 'fig2b_source_data']


def run_all(data_dir='figure_data'):
    """Run the full pipeline, pickle every plot-ready variable, and return them."""
    os.makedirs(data_dir, exist_ok=True)

    pct_df = load_unproductive_pct()
    gene_df = load_gene_tpm()
    upf_df = build_upf_df(gene_df, pct_df)

    fig2a_panels, fig2a_tissue_names = make_fig2a_panels(upf_df, pct_df)
    fig2e_boxplot_df = make_boxplot_df(upf_df)
    fig2e_boxplot_stats = make_boxplot_stats(fig2e_boxplot_df)

    ds_dge = load_ds_dge()
    fig2b_series = make_fig2b_series(ds_dge)
    fig2b_source_data = make_fig2b_source_data(ds_dge)

    data = {
        'fig2a_panels': fig2a_panels,
        'fig2a_tissue_names': fig2a_tissue_names,
        'fig2e_boxplot_df': fig2e_boxplot_df,
        'fig2e_boxplot_stats': fig2e_boxplot_stats,
        'fig2b_series': fig2b_series,
        'fig2b_source_data': fig2b_source_data,
    }

    for name, value in data.items():
        with open(os.path.join(data_dir, f'{name}.pickle'), 'wb') as fh:
            pickle.dump(value, fh)

    # Exact per-point statistics for Fig. 2b -- this is Supplementary Table 3.
    fig2b_source_data.to_csv(os.path.join(data_dir, 'fig2b_source_data.tsv'),
                             sep='\t', index=False)

    return data


def load_plot_data(data_dir='figure_data'):
    """Load every plot-ready variable back from `data_dir` (no recomputation)."""
    data = {}
    for name in PLOT_READY_VARS:
        with open(os.path.join(data_dir, f'{name}.pickle'), 'rb') as fh:
            data[name] = pickle.load(fh)
    return data
