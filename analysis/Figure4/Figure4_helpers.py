"""
Data loading and processing for Figure 4.

Created with Figure code cleaner.
Source notebooks: ../QTL_analysis.ipynb, ../hyprcoloc_results.ipynb,
                  ../coloc_plots.ipynb, ../Fig4_example.ipynb

Panels served by this module, named to match the manuscript caption:
  * fig4A          (a)  QTL_analysis.ipynb
  * fig4B_usQTL    (b)  QTL_analysis.ipynb   -- was fig4B1
  * fig4B_psQTL    (b)  QTL_analysis.ipynb   -- was fig4B2
  * fig4C          (c)  QTL_analysis.ipynb   -- was fig4B
  * fig4D          (d)  hyprcoloc_results.ipynb -- was fig4D2_coloc
  * fig4E_sQTL     (e)  coloc_plots.ipynb    -- was fig4_coloc_example_sQTL_colored
  * fig4E_eQTL     (e)  coloc_plots.ipynb    -- was fig4_coloc_example_eQTL_colored
  * fig4E_ASB16    (e)  coloc_plots.ipynb    -- was fig4_ASB16
  * fig4F          (f)  Fig4_example.ipynb   -- was fig4_coloc_example
  * fig4_colocs_sc2  not a Figure 4 panel; hyprcoloc_results.ipynb

Only standard packages are imported at module level; heavy/optional packages
(tabix, vcf, scipy, tqdm, rpy2) are imported inside the functions that use them,
so `load_plot_data` works in an environment without them.
"""

import os
import gzip
import pickle

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #

BASE = '/project/yangili1/cfbuenabadn/leafcutter2_paper'

# Panels are written here by the notebook.
PLOTS_DIR = f'{BASE}/analysis/Figure4/plots'

TOTAL_SQTLS_TABLE = f'{BASE}/code/analysis_files/Total_sQTLs.tsv.gz'
# QTL_analysis.ipynb read this as 'sQTL_stats.tsv.gz' from analysis/; the file
# now lives in code/analysis_files/.
SQTL_STATS_TABLE = f'{BASE}/code/analysis_files/sQTL_stats.tsv.gz'

SQTL_DIR = f'{BASE}/code/results/sqtl/GTEx'
EQTL_DIR = f'{BASE}/code/results/eqtl/GTEx'
NOISY_PHENO_DIR = '/project/yangili1/cfbuenabadn/SpliFi/code/results/pheno/noisy/GTEx'

HYPRCOLOC_TABLE = f'{BASE}/code/results/coloc/hyprcoloc_results/tables/hyprcoloc_results_filtered.tsv.gz'
UP_RESULTS_TABLE = f'{BASE}/code/results/coloc/hyprcoloc_results/tables/up_results.tsv.gz'
LEAD_SNPS_BED = f'{BASE}/code/resources/gwas/LeadSnpWindows.bed'

COLOC_DATA_DIR = f'{BASE}/code/results/coloc/data'
GENCODE_EXONS = ('/project2/mstephens/cfbuenabadn/gtex-stm/code/Annotations/'
                 'gencode.v44.primary_assembly.exons.bed.gz')
GTEX_VCF = ('/project/yangili1/cdai/genome_index/hs38/GTEx_v7/'
            'GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_'
            'WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz')

# Tissue used for the fibroblast scatter panels (fig4B_usQTL / fig4B_psQTL)
FIBROBLASTS = 'Cells-Culturedfibroblasts'

# --------------------------------------------------------------------------- #
# Constants recovered verbatim (GTEx colours, tissue and trait names, the ASB16
# example locus, and the GWAS accession table for fig4D).
# --------------------------------------------------------------------------- #

gtex_colors = {'Adipose-Subcutaneous': {'tissue_abbrv': 'ADPSBQ',
                          'tissue_color_hex': 'FFA54F',
                          'tissue_color_rgb': '255,165,79'},
 'Adipose-Visceral_Omentum': {'tissue_abbrv': 'ADPVSC',
                              'tissue_color_hex': 'EE9A00',
                              'tissue_color_rgb': '238,154,0'},
 'AdrenalGland': {'tissue_abbrv': 'ADRNLG',
                  'tissue_color_hex': '8FBC8F',
                  'tissue_color_rgb': '143,188,143'},
 'Artery-Aorta': {'tissue_abbrv': 'ARTAORT',
                  'tissue_color_hex': '8B1C62',
                  'tissue_color_rgb': '139,28,98'},
 'Artery-Coronary': {'tissue_abbrv': 'ARTCRN',
                     'tissue_color_hex': 'EE6A50',
                     'tissue_color_rgb': '238,106,80'},
 'Artery-Femoral': {'tissue_abbrv': 'ARTFMR',
                    'tissue_color_hex': 'FF4500',
                    'tissue_color_rgb': '255,69,0'},
 'Artery-Tibial': {'tissue_abbrv': 'ARTTBL',
                   'tissue_color_hex': 'FF0000',
                   'tissue_color_rgb': '255,0,0'},
 'Bladder': {'tissue_abbrv': 'BLDDER',
             'tissue_color_hex': 'CDB79E',
             'tissue_color_rgb': '205,183,158'},
 'Brain-Amygdala': {'tissue_abbrv': 'BRNAMY',
                    'tissue_color_hex': 'EEEE00',
                    'tissue_color_rgb': '238,238,0'},
 'Brain-Anteriorcingulatecortex_BA24': {'tissue_abbrv': 'BRNACC',
                                        'tissue_color_hex': 'EEEE00',
                                        'tissue_color_rgb': '238,238,0'},
 'Brain-Caudate_basalganglia': {'tissue_abbrv': 'BRNCDT',
                                'tissue_color_hex': 'EEEE00',
                                'tissue_color_rgb': '238,238,0'},
 'Brain-CerebellarHemisphere': {'tissue_abbrv': 'BRNCHB',
                                'tissue_color_hex': 'EEEE00',
                                'tissue_color_rgb': '238,238,0'},
 'Brain-Cerebellum': {'tissue_abbrv': 'BRNCHA',
                      'tissue_color_hex': 'EEEE00',
                      'tissue_color_rgb': '238,238,0'},
 'Brain-Cortex': {'tissue_abbrv': 'BRNCTXA',
                  'tissue_color_hex': 'EEEE00',
                  'tissue_color_rgb': '238,238,0'},
 'Brain-FrontalCortex_BA9': {'tissue_abbrv': 'BRNCTXB',
                             'tissue_color_hex': 'EEEE00',
                             'tissue_color_rgb': '238,238,0'},
 'Brain-Hippocampus': {'tissue_abbrv': 'BRNHPP',
                       'tissue_color_hex': 'EEEE00',
                       'tissue_color_rgb': '238,238,0'},
 'Brain-Hypothalamus': {'tissue_abbrv': 'BRNHPT',
                        'tissue_color_hex': 'EEEE00',
                        'tissue_color_rgb': '238,238,0'},
 'Brain-Nucleusaccumbens_basalganglia': {'tissue_abbrv': 'BRNNCC',
                                         'tissue_color_hex': 'EEEE00',
                                         'tissue_color_rgb': '238,238,0'},
 'Brain-Putamen_basalganglia': {'tissue_abbrv': 'BRNPTM',
                                'tissue_color_hex': 'EEEE00',
                                'tissue_color_rgb': '238,238,0'},
 'Brain-Spinalcord_cervicalc-1': {'tissue_abbrv': 'BRNSPC',
                                  'tissue_color_hex': 'EEEE00',
                                  'tissue_color_rgb': '238,238,0'},
 'Brain-Substantianigra': {'tissue_abbrv': 'BRNSNG',
                           'tissue_color_hex': 'EEEE00',
                           'tissue_color_rgb': '238,238,0'},
 'Breast-MammaryTissue': {'tissue_abbrv': 'BREAST',
                          'tissue_color_hex': '00CDCD',
                          'tissue_color_rgb': '0,205,205'},
 'Cells-EBV-transformedlymphocytes': {'tissue_abbrv': 'LCL',
                                      'tissue_color_hex': 'EE82EE',
                                      'tissue_color_rgb': '238,130,238'},
 'Cells-Culturedfibroblasts': {'tissue_abbrv': 'FIBRBLS',
                               'tissue_color_hex': '9AC0CD',
                               'tissue_color_rgb': '154,192,205'},
 'Cervix-Ectocervix': {'tissue_abbrv': 'CVXECT',
                       'tissue_color_hex': 'EED5D2',
                       'tissue_color_rgb': '238,213,210'},
 'Cervix-Endocervix': {'tissue_abbrv': 'CVSEND',
                       'tissue_color_hex': 'EED5D2',
                       'tissue_color_rgb': '238,213,210'},
 'Colon-Sigmoid': {'tissue_abbrv': 'CLNSGM',
                   'tissue_color_hex': 'CDB79E',
                   'tissue_color_rgb': '205,183,158'},
 'Colon-Transverse': {'tissue_abbrv': 'CLNTRN',
                      'tissue_color_hex': 'EEC591',
                      'tissue_color_rgb': '238,197,145'},
 'Esophagus-GastroesophagealJunction': {'tissue_abbrv': 'ESPGEJ',
                                        'tissue_color_hex': '8B7355',
                                        'tissue_color_rgb': '139,115,85'},
 'Esophagus-Mucosa': {'tissue_abbrv': 'ESPMCS',
                      'tissue_color_hex': '8B7355',
                      'tissue_color_rgb': '139,115,85'},
 'Esophagus-Muscularis': {'tissue_abbrv': 'ESPMSL',
                          'tissue_color_hex': 'CDAA7D',
                          'tissue_color_rgb': '205,170,125'},
 'FallopianTube': {'tissue_abbrv': 'FLLPNT',
                   'tissue_color_hex': 'EED5D2',
                   'tissue_color_rgb': '238,213,210'},
 'Heart-AtrialAppendage': {'tissue_abbrv': 'HRTAA',
                           'tissue_color_hex': 'B452CD',
                           'tissue_color_rgb': '180,82,205'},
 'Heart-LeftVentricle': {'tissue_abbrv': 'HRTLV',
                         'tissue_color_hex': '7A378B',
                         'tissue_color_rgb': '122,55,139'},
 'Kidney-Cortex': {'tissue_abbrv': 'KDNCTX',
                   'tissue_color_hex': 'CDB79E',
                   'tissue_color_rgb': '205,183,158'},
 'Kidney-Medulla': {'tissue_abbrv': 'KDNMDL',
                    'tissue_color_hex': 'CDB79E',
                    'tissue_color_rgb': '205,183,158'},
 'Liver': {'tissue_abbrv': 'LIVER',
           'tissue_color_hex': 'CDB79E',
           'tissue_color_rgb': '205,183,158'},
 'Lung': {'tissue_abbrv': 'LUNG',
          'tissue_color_hex': '9ACD32',
          'tissue_color_rgb': '154,205,50'},
 'MinorSalivaryGland': {'tissue_abbrv': 'SLVRYG',
                        'tissue_color_hex': 'CDB79E',
                        'tissue_color_rgb': '205,183,158'},
 'Muscle-Skeletal': {'tissue_abbrv': 'MSCLSK',
                     'tissue_color_hex': '7A67EE',
                     'tissue_color_rgb': '122,103,238'},
 'Nerve-Tibial': {'tissue_abbrv': 'NERVET',
                  'tissue_color_hex': 'FFD700',
                  'tissue_color_rgb': '255,215,0'},
 'Ovary': {'tissue_abbrv': 'OVARY',
           'tissue_color_hex': 'FFB6C1',
           'tissue_color_rgb': '255,182,193'},
 'Pancreas': {'tissue_abbrv': 'PNCREAS',
              'tissue_color_hex': 'CD9B1D',
              'tissue_color_rgb': '205,155,29'},
 'Pituitary': {'tissue_abbrv': 'PTTARY',
               'tissue_color_hex': 'B4EEB4',
               'tissue_color_rgb': '180,238,180'},
 'Prostate': {'tissue_abbrv': 'PRSTTE',
              'tissue_color_hex': 'D9D9D9',
              'tissue_color_rgb': '217,217,217'},
 'Skin-NotSunExposed_Suprapubic': {'tissue_abbrv': 'SKINNS',
                                   'tissue_color_hex': '3A5FCD',
                                   'tissue_color_rgb': '58,95,205'},
 'Skin-SunExposed_Lowerleg': {'tissue_abbrv': 'SKINS',
                              'tissue_color_hex': '1E90FF',
                              'tissue_color_rgb': '30,144,255'},
 'SmallIntestine-TerminalIleum': {'tissue_abbrv': 'SNTTRM',
                                  'tissue_color_hex': 'CDB79E',
                                  'tissue_color_rgb': '205,183,158'},
 'Spleen': {'tissue_abbrv': 'SPLEEN',
            'tissue_color_hex': 'CDB79E',
            'tissue_color_rgb': '205,183,158'},
 'Stomach': {'tissue_abbrv': 'STMACH',
             'tissue_color_hex': 'FFD39B',
             'tissue_color_rgb': '255,211,155'},
 'Testis': {'tissue_abbrv': 'TESTIS',
            'tissue_color_hex': 'A6A6A6',
            'tissue_color_rgb': '166,166,166'},
 'Thyroid': {'tissue_abbrv': 'THYROID',
             'tissue_color_hex': '008B45',
             'tissue_color_rgb': '0,139,69'},
 'Uterus': {'tissue_abbrv': 'UTERUS',
            'tissue_color_hex': 'EED5D2',
            'tissue_color_rgb': '238,213,210'},
 'Vagina': {'tissue_abbrv': 'VAGINA',
            'tissue_color_hex': 'EED5D2',
            'tissue_color_rgb': '238,213,210'},
 'WholeBlood': {'tissue_abbrv': 'WHLBLD',
                'tissue_color_hex': 'FF00FF',
                'tissue_color_rgb': '255,0,255'}}

FIG4A_TISSUE_NAMES = ['Testis',
 'Thyroid',
 'Nerve - Tibial',
 'Adipose - Subcutaneous',
 'Skin - Sun exposed (lower leg)',
 'Artery - Tibial',
 'Cells - Cultured fibroblasts',
 'Lung',
 'Skin - Not sun exposed (suprapubic)',
 'Muscle - Skeletal',
 'Adipose - Visceral (omentum)',
 'Esophagus - Muscularis',
 'Esophagus - Mucosa',
 'Breast - Mammary tissue',
 'Artery - Aorta',
 'Colon - Transverse',
 'Esophagus - Gastroesophageal junction',
 'Whole blood',
 'Colon - Sigmoid',
 'Heart - Atrial appendage',
 'Pituitary',
 'Stomach',
 'Spleen',
 'Brain - Cerebellum',
 'Heart - Left ventricle',
 'Adrenal gland',
 'Prostate',
 'Pancreas',
 'Brain - Cerebellar hemisphere',
 'Cells - EBV-transformed lymphocytes',
 'Artery - Coronary',
 'Brain - Cortex',
 'Ovary',
 'Small intestine - Terminal ileum',
 'Brain - Nucleus accumbens (basal ganglia)',
 'Brain - Caudate (basal ganglia)',
 'Minor salivary gland',
 'Brain - Frontal cortex (BA9)',
 'Uterus',
 'Vagina',
 'Liver',
 'Brain - Hypothalamus',
 'Brain - Putamen (basal ganglia)',
 'Brain - Anterior cingulate cortex (BA24)',
 'Brain - Hippocampus',
 'Brain - Spinalcord (cervical c-1)',
 'Brain - Amygdala',
 'Brain - Substantia nigra',
 'Kidney - Cortex']

FIG4C_TISSUE_NAMES = ['Brain - Cerebellum',
 'Spleen',
 'Brain - Cerebellar hemisphere',
 'Brain - Substantia nigra',
 'Kidney - Cortex',
 'Brain - Amygdala',
 'Small intestine - Terminal ileum',
 'Brain - Nucleus accumbens (basal ganglia)',
 'Brain - Anterior cingulate cortex (BA24)',
 'Brain - Frontal cortex (BA9)',
 'Pituitary',
 'Brain - Hippocampus',
 'Prostate',
 'Whole blood',
 'Brain - Spinal cord (cervical c-1)',
 'Brain - Hypothalamus',
 'Liver',
 'Uterus',
 'Brain - Cortex',
 'Brain - Putamen (basal ganglia)',
 'Vagina',
 'Brain - Caudate (basal ganglia)',
 'Colon - Sigmoid',
 'Thyroid',
 'Esophagus - Gastroesophageal junction',
 'Heart - Left ventricle',
 'Nerve - Tibial',
 'Artery - Coronary',
 'Ovary',
 'Heart - Atrial appendage',
 'Lung',
 'Stomach',
 'Pancreas',
 'Colon - Transverse',
 'Esophagus - Muscularis',
 'Adipose - Subcutaneous',
 'Skin - Sun exposed (lower leg)',
 'Minor salivary gland',
 'Breast - Mammary tissue',
 'Adrenal gland',
 'Adipose - Visceral (omentum)',
 'Artery - Aorta',
 'Skin - Not sun exposed (suprapubic)',
 'Testis',
 'Cells - EBV-transformed lymphocytes',
 'Artery - Tibial',
 'Esophagus - Mucosa',
 'Muscle - Skeletal',
 'Cells - Cultured fibroblasts']

FIG4D_TRAIT_NAMES = ['Bipolar disorder',
 'Multiple sclerosis',
 'Heart failure',
 'Coronary artery disease',
 'Myocardial infarction',
 'Visceral adipose tissue',
 'Age when finished education',
 'Inflammatory bowel disease',
 'Schizophrenia',
 'Asthma - childhood onset',
 'Ulcerative colitis',
 'Hypothyroidism',
 'Chronic obstructive pulmonary disease',
 'Basal cell carcinoma',
 'Breast cancer',
 "Dupuytren's disease",
 "Crohn's disease",
 'Atrial fibrillation',
 'Atopic eczema']

FIG4D_TRAIT_ACCESSIONS = {'bipolar_disorder': {'trait': 'Bipolar disorder',
                      'trait_file': 'Bipolar_disorder',
                      'accession': 'bip2021',
                      'source': 'Torino',
                      'assembly': 'GRCh38'},
 'IMSGC2019': {'trait': 'Multiple sclerosis',
               'trait_file': 'Multiple_sclerosis',
               'accession': 'IMSGC2019',
               'source': 'Ben',
               'assembly': 'GRCh38'},
 'heart_failure': {'trait': 'Heart failure',
                   'trait_file': 'Heart_failure',
                   'accession': 'GCST90162626',
                   'source': 'Torino',
                   'assembly': 'GRCh38'},
 'coronary_artery_disease': {'trait': 'Coronary artery disease',
                             'trait_file': 'Coronary_artery_disease',
                             'accession': 'GCST90132314',
                             'source': 'Torino',
                             'assembly': 'GRCh38'},
 'myocardial_infarction': {'trait': 'Myocardial infarction',
                           'trait_file': 'Myocardial_infarction',
                           'accession': 'GCST011365',
                           'source': 'Torino',
                           'assembly': 'GRCh38'},
 'Visceral_adipose_tissue_measurement': {'trait': 'Visceral adipose tissue',
                                         'trait_file': 'Visceral_adipose_tissue_measurement',
                                         'accession': 'GCST008744',
                                         'source': 'Leafcutter2',
                                         'assembly': 'GRCh37',
                                         'author': 'Karlsson T',
                                         'pubmed': '31501611'},
 'age_when_finished_full-time_education': {'trait': 'Age when finished education',
                                           'trait_file': 'Age_when_finished_full-time_education',
                                           'accession': 'GCST90267280',
                                           'source': 'Torino',
                                           'assembly': 'GRCh38'},
 'GCST004131': {'trait': 'Inflammatory bowel disease',
                'trait_file': 'Inflammatory_bowel_disease',
                'accession': 'GCST004131',
                'source': 'Ben',
                'assembly': 'GRCh38'},
 'schizophrenia': {'trait': 'Schizophrenia',
                   'trait_file': 'Schizophrenia',
                   'accession': 'scz2022',
                   'source': 'Torino',
                   'assembly': 'GRCh38'},
 'GCST007800': {'trait': 'Asthma - childhood onset',
                'trait_file': 'Asthma_childhood_onset',
                'accession': 'GCST007800',
                'source': 'Ben',
                'assembly': 'GRCh38'},
 'GCST004133': {'trait': 'Ulcerative colitis',
                'trait_file': 'Ulcerative_colitis',
                'accession': 'GCST004133',
                'source': 'Ben',
                'assembly': 'GRCh38'},
 'Hypothyroidism': {'trait': 'Hypothyroidism',
                    'trait_file': 'Hypothyroidism',
                    'accession': 'GCST90319320',
                    'source': 'Leafcutter2',
                    'assembly': 'GRCh38',
                    'author': 'Figueredo J',
                    'pubmed': '39067062'},
 'Chronic_obstructive_pulmonary_disease': {'trait': 'Chronic obstructive pulmonary '
                                                    'disease',
                                           'trait_file': 'Chronic_obstructive_pulmonary_disease',
                                           'accession': 'GCST90244098',
                                           'source': 'Leafcutter2',
                                           'assembly': 'GRCh38',
                                           'author': 'Cosentino J',
                                           'pubmed': '37069358'},
 'basal_cell_carcinoma': {'trait': 'Basal cell carcinoma',
                          'trait_file': 'Basal_cell_carcinoma',
                          'accession': 'GCST90137411',
                          'source': 'Torino',
                          'assembly': 'GRCh38'},
 'GCST004988': {'trait': 'Breast cancer',
                'trait_file': 'Breast_cancer',
                'accession': 'GCST004988',
                'source': 'Ben',
                'assembly': 'GRCh38'},
 'Dupuytrens_disease': {'trait': "Dupuytren's disease",
                        'trait_file': 'Dupuytrens_disease',
                        'accession': 'GCST90301252',
                        'source': 'Leafcutter2',
                        'assembly': 'GRCh38',
                        'author': 'Riesmeijer SA',
                        'pubmed': '38172110'},
 'GCST004132': {'trait': "Crohn's disease",
                'trait_file': 'Crohns_disease',
                'accession': 'GCST004132',
                'source': 'Ben',
                'assembly': 'GRCh38'},
 'atrial_fibrillation': {'trait': 'Atrial fibrillation',
                         'trait_file': 'Atrial_fibrillation',
                         'accession': 'GCST006061',
                         'source': 'Torino',
                         'assembly': 'GRCh38'},
 'Atopic_eczema': {'trait': 'Atopic eczema',
                   'trait_file': 'Atopic_eczema',
                   'accession': 'GCST90244787',
                   'source': 'Leafcutter2',
                   'assembly': 'GRCh37',
                   'author': 'Budu-Aggrey A',
                   'pubmed': '37794016'},
 'Rheumatoid_arthritis': {'trait': 'Rheumatoid arthritis',
                          'trait_file': 'Rheumatoid_arthritis',
                          'accession': 'GCST002318',
                          'source': 'Leafcutter2',
                          'assembly': 'GRCh38',
                          'author': 'Okada Y',
                          'pubmed': '24390342'}}

FIG4D_UNRESOLVED_ACCESSIONS = ['bip2021', 'scz2022', 'IMSGC2019']

FIG4B_usQTL_SELECTION = ('unproductive (UP) introns in clusters containing both productive and unproductive '
 "introns (ctype 'PR,UP'), with a significant sQTL (Storey q <= 0.1 from the QTLtools "
 "permutation pass) and a cluster accounting for >= 10% of the gene's junction reads")

FIG4B_psQTL_SELECTION = ("productive (PR) introns in wholly productive clusters (ctype 'PR'), with a "
 'significant sQTL (Storey q <= 0.1 from the QTLtools permutation pass) and a cluster '
 "accounting for >= 10% of the gene's junction reads")

FIG4B_UNIT = ("one intron-variant-gene triplet: an intron passing the panel's selection, its lead "
 "sQTL variant, and that variant's effect on the host gene's expression. These are NOT "
 'replicates')

FIG4B_DONOR_UNIT = ('GTEx donors with cultured fibroblasts; one donor = one primary culture = one library '
 '= one genotype, i.e. biological replicates, no technical replicates')

ASB16_TISSUE = 'Brain-CerebellarHemisphere'

ASB16_CHROM = 'chr17'

ASB16_START = 44176940

ASB16_END = 44177608

ASB16_SNP = 'chr17:44176913'

ASB16_GENE = 'ENSG00000161664'

ASB16_SPLICE_PHE = 'chr17:44176940:44177608:clu_28859_+:UP:chr17_44114525_N_N_bipolar_disorder'

ASB16_EXPR_PHE = 'ENSG00000161664.7:chr17_44114525_N_N_bipolar_disorder'

ASB16_TRANSCRIPTS = ['ENST00000293414', 'ENST00000589618']

ASB16_NOM_DIR = '/project/yangili1/cfbuenabadn/leafcutter2_paper/code/results/coloc/qtls'

ASB16_GWAS_TRAIT = 'chr17_44114525_N_N_bipolar_disorder'

ASB16_LEAD_LABEL = 'rs7212573'

BIPOLAR_STATS = '/project/yangili1/cfbuenabadn/leafcutter2_paper/code/resources/gwas/StatsForColoc/Bipolar_disorder.standardized.txt.tabix.gz'

BIPOLAR_PVALS = '/project2/mstephens/cfbuenabadn/gtex-stm/code/gwas/hg38_summary_stats/bipolar_disorder.bed.gz'

ASB16_LD = '/project/yangili1/cfbuenabadn/leafcutter2_paper/code/results/coloc/LD/chr17_44114525_N_N_bipolar_disorder.tsv.gz'


def tissue_color(tissue):
    """'#' + the GTEx hex colour for `tissue`."""
    return '#' + gtex_colors[tissue]['tissue_color_hex']


# =========================================================================== #
# fig4A -- total sQTLs per tissue                        (QTL_analysis.ipynb)
# =========================================================================== #

def load_total_sqtls(table=TOTAL_SQTLS_TABLE):
    """Per-tissue p-/u-/other sQTL counts, as tabulated by QTL_analysis.ipynb."""
    return pd.read_csv(table, sep='\t', index_col=0)


def get_sample_sizes(tissues):
    """Number of donors behind each tissue's permutation pass (dof1 + 2)."""
    from tqdm import tqdm

    n_samples = []
    for tissue in tqdm(tissues, position=0, leave=True):
        perm = pd.read_csv(
            f'{SQTL_DIR}/{tissue}/cis_100000/perm/PermutationPass.Qval.txt.gz', sep='\t')
        n_samples.append((perm.dof1 + 2).iloc[0])
    return n_samples


def get_donor_count(tissue):
    """Number of GTEx donors behind one tissue's QTL mapping (dof1 + 2).

    This is the number of **biological** replicates for that tissue: one
    post-mortem donor contributes one sample, one library and one set of
    genotypes. There are no technical replicates in this design.
    """
    perm = pd.read_csv(
        f'{SQTL_DIR}/{tissue}/cis_100000/perm/PermutationPass.Qval.txt.gz', sep='\t')
    return int((perm.dof1 + 2).iloc[0])


def get_all_donor_counts():
    """{tissue: number of donors} for all 49 GTEx tissues, from the permutation pass."""
    tissues = sorted(os.listdir(SQTL_DIR))
    return {t: n for t, n in zip(tissues, get_sample_sizes(tissues))}


def make_fig4A_data(donor_counts=None):
    """Plot-ready arrays for fig4A, ordered by total sQTLs (descending)."""
    qtl_dataframe = load_total_sqtls()
    sorted_tissues = list(qtl_dataframe.sort_values(by='total', ascending=False).index)
    if donor_counts is None:
        donor_counts = get_all_donor_counts()

    return {
        'sorted_tissues': sorted_tissues,
        'tissue_names': FIG4A_TISSUE_NAMES,
        'colors': [tissue_color(t) for t in sorted_tissues],
        'pr_counts': np.array(qtl_dataframe.loc[sorted_tissues, 'PR']),
        'up_counts': np.array(qtl_dataframe.loc[sorted_tissues, 'UP']),
        'n_samples': [donor_counts[t] for t in sorted_tissues],
        'donor_counts': donor_counts,
    }


# =========================================================================== #
# fig4C -- sQTL vs eQTL effect-size correlation per tissue (QTL_analysis.ipynb)
# =========================================================================== #

def load_sqtl_stats(table=SQTL_STATS_TABLE):
    """Per-tissue sQTL-vs-eQTL correlation statistics tabulated by QTL_analysis.ipynb."""
    return pd.read_csv(table, sep='\t')


def make_fig4C_data():
    """Plot-ready arrays for fig4C, ordered by u-sQTL Spearman p-value (weakest first)."""
    sqtl_stats = load_sqtl_stats()

    u = sqtl_stats.loc[sqtl_stats.sqtl_type == 'u_sqtl']
    sorted_tissues = list(u.sort_values('spearman_pval').tissue)[::-1]

    # Significance bins, from the most to the least significant end of the ranking
    n1 = len(u.loc[u.spearman_pval <= 1e-10].tissue)
    n2 = len(u.loc[u.spearman_pval <= 1e-4].tissue) - n1
    n3 = len(u.loc[u.spearman_pval <= 5e-2].tissue) - n1 - n2
    n4 = len(u.loc[u.spearman_pval > 5e-2].tissue)

    rho_u = np.array(
        sqtl_stats.loc[sqtl_stats.sqtl_type == 'u_sqtl'].set_index('tissue')
        .loc[sorted_tissues].spearman_rho)
    rho_pp = np.array(
        sqtl_stats.loc[sqtl_stats.sqtl_type == 'pp_sqtl'].set_index('tissue')
        .loc[sorted_tissues].spearman_rho)

    return {
        'sorted_tissues': sorted_tissues,
        'tissue_names': FIG4C_TISSUE_NAMES,
        'colors': [tissue_color(t) for t in sorted_tissues],
        'rho_u': rho_u,
        'rho_pp': rho_pp,
        'bins': (n1, n2, n3, n4),
    }


def spearman_ci95(rho, n):
    """95% CI for Spearman's rho: Fisher z with the Bonett-Wright standard error."""
    z = np.arctanh(rho)
    se = np.sqrt((1 + rho ** 2 / 2) / (n - 3))
    zc = 1.959963984540054
    return float(np.tanh(z - zc * se)), float(np.tanh(z + zc * se))


def make_fig4C_source_data(donor_counts=None):
    """Per-tissue source-data table for fig4C -- one row per tissue per sQTL class.

    fig4C plots 49 tissues x 2 classes, so the per-point n, exact P, test
    statistic and confidence interval cannot fit in a legend. This is the table
    the legend points to (Supplementary Table 9). Columns: tissue, sqtl_class,
    n_donors, n_pairs, spearman_rho, t, df, exact P, and the 95% CI for rho.
    """
    sqtl_stats = load_sqtl_stats()
    if donor_counts is None:
        donor_counts = get_all_donor_counts()

    keep = sqtl_stats.loc[sqtl_stats.sqtl_type.isin(['u_sqtl', 'pp_sqtl'])].copy()
    keep['sqtl_class'] = keep.sqtl_type.map({'u_sqtl': 'u-sQTL (unproductive intron)',
                                             'pp_sqtl': 'p-sQTL (productive intron)'})
    keep['n_donors'] = keep.tissue.map(donor_counts)
    keep['n_pairs'] = keep.n.astype(int)
    keep['df'] = keep.n_pairs - 2
    keep['t'] = keep.spearman_rho * np.sqrt(keep.df / (1 - keep.spearman_rho ** 2))
    ci = [spearman_ci95(r, n) for r, n in zip(keep.spearman_rho, keep.n_pairs)]
    keep['ci95_low'] = [c[0] for c in ci]
    keep['ci95_high'] = [c[1] for c in ci]

    out = keep[['tissue', 'sqtl_class', 'n_donors', 'n_pairs', 'spearman_rho', 't',
                'df', 'spearman_pval', 'ci95_low', 'ci95_high']]
    out = out.rename(columns={'spearman_pval': 'exact_two_sided_P'})
    return out.sort_values(['tissue', 'sqtl_class']).reset_index(drop=True)


# =========================================================================== #
# fig4B_usQTL / fig4B_psQTL -- fibroblast sQTL vs eQTL betas
#                                                        (QTL_analysis.ipynb)
# =========================================================================== #

def run_tabix_on_nom(nom_file, coords, transform_to_df=True):
    """QTLtools nominal-pass records overlapping `coords` ('chrom:start:end')."""
    import tabix

    chrom, start, end = coords.split(':')
    start = int(start)
    end = int(end)

    tb = tabix.open(nom_file)
    nom = tb.query(chrom, start, end)

    if transform_to_df:
        columns = ['#phe_id', 'phe_chr', 'phe_from', 'phe_to', 'phe_strd', 'n_var_in_cis',
                   'dist_phe_var', 'var_id', 'var_chr', 'var_from', 'var_to', 'nom_pval',
                   'r_squared', 'slope', 'slope_se', 'best_hit']

        nom_bed = []
        nom_bed_cols = []

        for idx, record in enumerate(nom):
            nom_bed.append(pd.Series(record))
            nom_bed_cols.append(f'record{str(idx)}')
        nom_bed = pd.concat(nom_bed, axis=1)
        nom_bed.columns = nom_bed_cols
        nom_bed = nom_bed.T
        nom_bed.columns = columns
        nom_bed['#phe_id'] = nom_bed['#phe_id'].apply(lambda x: x.split('.')[0])
        return nom_bed

    return nom


def get_expression_effect(tissue, perm_qtl):
    """eQTL nominal record for the gene / variant of one sQTL permutation hit."""
    var_id = perm_qtl.var_id
    chrom = perm_qtl.var_chr
    start = str(int(perm_qtl.var_from) - 1)
    end = str(int(perm_qtl.var_to) + 1)

    coords = f'{chrom}:{start}:{end}'
    gene = perm_qtl.gene_id.split('.')[0]

    nom_file = f'{EQTL_DIR}/{tissue}/cis_100000/nom/{chrom}.txt.gz'
    nom_df = run_tabix_on_nom(nom_file, coords, transform_to_df=True)
    nom_df = nom_df.loc[(nom_df['#phe_id'] == gene) & (nom_df['var_id'] == var_id)]

    return nom_df


def get_intron_list(tissue):
    """Median read count per intron, with its cluster, from the leafcutter2 noise table."""
    intron_list = []
    cluster_list = []
    median_count_list = []
    with gzip.open(
        f'{NOISY_PHENO_DIR}/{tissue}/leafcutter_perind_numers.counts.noise_by_intron.gz'
    ) as fh:
        fh.readline()
        for line in fh:
            line = line.decode().rstrip().split(' ')
            intron_list.append(line[0])
            cluster_list.append(line[0].split(':')[0] + ':' + line[0].split(':')[-2])
            median_counts = int(np.median([int(y) for y in line[1:]]))
            median_count_list.append(median_counts)

    intron_median_counts = pd.DataFrame()
    intron_median_counts['phenotype_id'] = intron_list
    intron_median_counts['cluster'] = cluster_list
    intron_median_counts['median_counts'] = median_count_list
    intron_median_counts['intron'] = intron_median_counts.phenotype_id.apply(
        lambda x: ':'.join(x.split(':')[:-1]))

    return intron_median_counts


def get_perm_counts(tissue):
    """sQTL permutation pass joined to intron counts, with intron/cluster PSI columns."""
    perm = pd.read_csv(
        f'{SQTL_DIR}/{tissue}/cis_100000/perm/PermutationPass.Qval.txt.gz', sep='\t')
    intron_median_counts = get_intron_list(tissue)
    perm_counts = pd.merge(perm, intron_median_counts,
                           right_on=['phenotype_id', 'intron', 'cluster'],
                           left_on=['phe_id', 'intron', 'cluster'])
    counts_per_cluster = perm_counts.groupby('cluster').median_counts.sum().reset_index()
    counts_per_cluster.columns = ['cluster', 'cluster_counts']

    counts_per_gene = perm_counts.groupby('gene_id').median_counts.sum().reset_index()
    counts_per_gene.columns = ['gene_id', 'gene_counts']

    perm_counts = pd.merge(perm_counts, counts_per_cluster, left_on='cluster', right_on='cluster')
    perm_counts = pd.merge(perm_counts, counts_per_gene, left_on='gene_id', right_on='gene_id')
    perm_counts['intron_psi'] = perm_counts.median_counts / (perm_counts.cluster_counts + 1e-10)
    perm_counts['cluster_psi'] = perm_counts.cluster_counts / (perm_counts.gene_counts + 1e-10)
    UP_per_cluster = perm_counts.groupby('cluster').itype.apply(
        lambda x: (x == 'UP').sum()).reset_index()
    UP_per_cluster.columns = ['cluster', 'UP_juncs']
    perm_counts = pd.merge(perm_counts, UP_per_cluster, left_on='cluster', right_on='cluster')

    return perm_counts


def get_sqtl_expression_effects(perm_counts, tissue, ctype='PR,UP', itype='UP',
                                qmax=1e-2, clu_psi_min=0.1):
    """Matched (sQTL beta, eQTL beta, eQTL nominal p) for every significant sQTL."""
    sqtl_slopes = []
    eqtl_slopes = []
    eqtl_pvals = []

    perm_select = perm_counts.loc[(perm_counts.ctype == ctype) & (perm_counts.itype == itype)
                                  & (perm_counts.q <= qmax)
                                  & (perm_counts.cluster_psi >= clu_psi_min)]

    for idx, row in perm_select.iterrows():
        try:
            eqtl_df = get_expression_effect(tissue, row)
            if eqtl_df.shape[0] > 0:
                eqtl_effect = float(eqtl_df.iloc[0].slope)
                sqtl_slopes.append(float(row.slope))
                eqtl_slopes.append(float(eqtl_effect))
                eqtl_pvals.append(float(eqtl_df.iloc[0].nom_pval))
        except:
            continue

    return sqtl_slopes, eqtl_slopes, eqtl_pvals


def make_rho_stats(x, y, n_donors=None, n_boot=10000, seed=0):
    """Every statistic the figure legend has to report for one beta-vs-beta panel.

    The test is a **two-sided Spearman rank correlation** on the full sample; no
    multiple-comparison adjustment is applied, because each panel is a single
    test (the multiple-testing correction in this analysis happens upstream, at
    sQTL selection, as a Storey q-value from the permutation pass).

    Returned keys
    -------------
    n            : number of data points entering the test (exact)
    n_donors     : number of GTEx donors behind the QTL mapping (biological units)
    rho          : Spearman's rho -- the effect size
    pvalue       : exact two-sided P value
    df           : degrees of freedom of the associated t statistic (n - 2)
    t            : test statistic, rho * sqrt(df / (1 - rho**2))
    ci95         : 95% CI for rho, Fisher z with the Bonett-Wright SE for Spearman
    ci95_boot    : 95% CI for rho, percentile bootstrap over `n_boot` resamples
    pearson_r    : Pearson r on the same points, for completeness
    pearson_p    : exact two-sided P value for Pearson r
    ols          : least-squares fit of y on x (slope, intercept, its own 95% CI,
                   two-sided P and standard error) -- reported so the panel can be
                   redrawn with a conventional regression line if an editor asks
    """
    from scipy.stats import spearmanr, pearsonr, linregress, t as t_dist

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)

    rho, pval = spearmanr(x, y)            # scipy default: alternative='two-sided'
    r_p, p_p = pearsonr(x, y)              # two-sided

    df = n - 2
    t_stat = rho * np.sqrt(df / (1 - rho ** 2))

    # Fisher z CI for Spearman's rho, Bonett & Wright (2000) standard error
    z = np.arctanh(rho)
    se_z = np.sqrt((1 + rho ** 2 / 2) / (n - 3))
    z_crit = 1.959963984540054
    ci95 = (float(np.tanh(z - z_crit * se_z)), float(np.tanh(z + z_crit * se_z)))

    # Percentile bootstrap CI over resampled data points
    rng = np.random.default_rng(seed)
    boot = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        boot[i] = spearmanr(x[idx], y[idx])[0]
    ci95_boot = (float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975)))

    lr = linregress(x, y)
    t_crit = t_dist.ppf(0.975, df)
    ols = {
        'slope': float(lr.slope),
        'intercept': float(lr.intercept),
        'stderr': float(lr.stderr),
        'pvalue': float(lr.pvalue),
        'ci95': (float(lr.slope - t_crit * lr.stderr),
                 float(lr.slope + t_crit * lr.stderr)),
    }

    return {
        'n': int(n),
        'n_donors': n_donors,
        'test': 'two-sided Spearman rank correlation, no multiple-comparison adjustment',
        'rho': float(rho),
        'pvalue': float(pval),
        'df': int(df),
        't': float(t_stat),
        'ci95': ci95,
        'ci95_boot': ci95_boot,
        'pearson_r': float(r_p),
        'pearson_p': float(p_p),
        'ols': ols,
    }


def format_stats_report(name, stats, unit, donor_unit):
    """Multi-line, legend-ready summary of `make_rho_stats` output."""
    lo, hi = stats['ci95']
    blo, bhi = stats['ci95_boot']
    ols = stats['ols']
    olo, ohi = ols['ci95']
    return '\n'.join([
        f'{name}',
        f'  unit of study : {unit}',
        f'  n (points)    : {stats["n"]}',
        f'  n (donors)    : {stats["n_donors"]}  ({donor_unit})',
        f'  test          : {stats["test"]}',
        f"  Spearman rho  : {stats['rho']:.4f}   (effect size)",
        f"  t({stats['df']})       : {stats['t']:.4f}",
        f"  exact P       : {stats['pvalue']:.4g}",
        f'  95% CI (rho)  : [{lo:.4f}, {hi:.4f}]   Fisher z, Bonett-Wright SE',
        f'  95% CI (boot) : [{blo:.4f}, {bhi:.4f}]   percentile bootstrap',
        f"  Pearson r     : {stats['pearson_r']:.4f}  (exact two-sided P = {stats['pearson_p']:.4g})",
        f"  OLS slope     : {ols['slope']:.4f} +/- {ols['stderr']:.4f} (SE), "
        f"95% CI [{olo:.4f}, {ohi:.4f}], exact P = {ols['pvalue']:.4g}",
    ])


def make_ols_band(x, y, x_grid):
    """Least-squares fit of y on x, with its 95% CI and 95% prediction interval.

    Evaluated on `x_grid` so the notebook's plot cells never compute anything.
    The CI band is the textbook one for the *fitted mean*,

        yhat(x) +/- t(0.975, n-2) * s * sqrt(1/n + (x - xbar)^2 / Sxx)

    and the prediction interval replaces `1/n` with `1 + 1/n`. At the sample
    sizes here the CI band is a hairline; the prediction interval is what shows
    how much scatter there actually is around the fit. They answer different
    questions and the legend must say which one is drawn.
    """
    from scipy.stats import t as t_dist

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x_grid = np.asarray(x_grid, dtype=float)
    n = len(x)

    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (intercept + slope * x)
    s_resid = np.sqrt(np.sum(resid ** 2) / (n - 2))

    xbar = x.mean()
    Sxx = np.sum((x - xbar) ** 2)
    t_crit = t_dist.ppf(0.975, n - 2)

    fitted = intercept + slope * x_grid
    se_mean = s_resid * np.sqrt(1 / n + (x_grid - xbar) ** 2 / Sxx)
    se_pred = s_resid * np.sqrt(1 + 1 / n + (x_grid - xbar) ** 2 / Sxx)

    return {
        'slope': float(slope),
        'intercept': float(intercept),
        's_resid': float(s_resid),
        'fitted': fitted,
        'ci_lo': fitted - t_crit * se_mean,
        'ci_hi': fitted + t_crit * se_mean,
        'pi_lo': fitted - t_crit * se_pred,
        'pi_hi': fitted + t_crit * se_pred,
    }


def make_rho_fit(x, y, n_draws=100, draw_size=100, pool_size=None, seed=None,
                 n_boot=2000, x_grid=None):
    """The dashed line and the shaded band drawn on the beta-vs-beta scatters.

    Note that neither is a conventional regression fit, and the legend has to say
    so. The dashed line has **slope equal to Spearman's rho** of the full sample
    and **intercept equal to mean(y)**; the default shaded band spans the
    **10th-90th percentile of Spearman's rho** over `n_draws` random subsamples of
    `draw_size` points drawn without replacement. That is the construction used in
    QTL_analysis.ipynb and hyprcoloc_results.ipynb, reproduced here verbatim.

    fig4B_usQTL and fig4B_psQTL are now drawn with `line='ols'` and no band, so
    for those panels only `ols` is used; the rho-slope construction is kept
    because the hyprcoloc panel (fig4_colocs_sc2) still reproduces it.

    `seed` makes the band reproducible -- the original notebooks left the RNG
    unseeded, so their band differed slightly from run to run. Set it.

    A defensible alternative band is also returned under the 'boot_2p5'/'boot_97p5'
    keys: the 2.5th-97.5th percentile of rho over `n_boot` bootstrap resamples of
    the full sample, i.e. a genuine 95% CI. `Figure4_plot_helpers.plot_beta_scatter`
    draws it when called with `band='bootstrap95'`. `n_boot` is lower here than in
    `make_rho_stats` because this one is only drawn, never quoted -- the CI that
    goes in the legend comes from `make_rho_stats`, at 10,000 resamples.
    """
    from scipy.stats import spearmanr

    x = np.array(x)
    y = np.array(y)
    if pool_size is None:
        pool_size = len(x)
    if x_grid is None:
        # matches the linspace the original notebooks drew the dashed line over
        x_grid = np.linspace(-3.2, 3.2, 100)
    x_grid = np.asarray(x_grid, dtype=float)

    rng = np.random.default_rng(seed) if seed is not None else np.random
    rho_list = []
    for i in range(n_draws):
        idx = rng.choice(np.arange(pool_size), draw_size, replace=False)
        rho_list.append(spearmanr(x[idx], y[idx])[0])

    # Genuine 95% CI on rho, for the 'bootstrap95' band
    boot_rng = np.random.default_rng(0 if seed is None else seed)
    n = len(x)
    boot = np.empty(n_boot)
    for i in range(n_boot):
        idx = boot_rng.integers(0, n, n)
        boot[i] = spearmanr(x[idx], y[idx])[0]

    return {
        'x': x,
        'y': y,
        'slope': spearmanr(x, y)[0],
        'intercept': np.mean(y),
        'slope_90': np.quantile(rho_list, 0.9),
        'slope_10': np.quantile(rho_list, 0.1),
        'boot_97p5': np.quantile(boot, 0.975),
        'boot_2p5': np.quantile(boot, 0.025),
        'x_grid': x_grid,
        'ols': make_ols_band(x, y, x_grid),
        # How the drawn band is defined -- copy this into the figure legend.
        'band': {
            'default': (f'10th-90th percentile of Spearman rho over {n_draws} random '
                        f'subsamples of {draw_size} points drawn without replacement'),
            'bootstrap95': (f'2.5th-97.5th percentile of Spearman rho over {n_boot} '
                            f'bootstrap resamples of all {n} points (95% CI)'),
            'line': ('slope = Spearman rho of the full sample; '
                     'intercept = mean of the y values'),
            'ols_ci95': ('least-squares fit of y on x, with the 95% confidence '
                         'interval of the fitted mean'),
            'ols_pi95': ('least-squares fit of y on x, with the 95% prediction '
                         'interval for a single new observation'),
            'n_draws': n_draws,
            'draw_size': draw_size,
            'n_boot': n_boot,
            'seed': seed,
        },
    }


def make_fibroblast_fits(seed=0):
    """fig4B_usQTL and fig4B_psQTL fits and legend statistics.

    Returns ``(usQTL_fit, psQTL_fit, usQTL_stats, psQTL_stats)``. The two panels
    share a sample of donors and differ only in which introns enter:
    fig4B_psQTL is the negative control for fig4B_usQTL -- productive introns,
    where no NMD-mediated coupling between splicing and expression is expected.
    """
    perm_counts = get_perm_counts(FIBROBLASTS)
    u_sqtl_effects = get_sqtl_expression_effects(
        perm_counts, FIBROBLASTS, ctype='PR,UP', itype='UP', qmax=1e-1, clu_psi_min=0.1)
    pp_sqtl_effects = get_sqtl_expression_effects(
        perm_counts, FIBROBLASTS, ctype='PR', itype='PR', qmax=1e-1, clu_psi_min=0.1)

    # QTL_analysis.ipynb hard-coded the u-sQTL band's subsample pool at 1391 points.
    usQTL_fit = make_rho_fit(u_sqtl_effects[0], u_sqtl_effects[1], pool_size=1391,
                             seed=seed, x_grid=np.linspace(-3.2, 3.2, 100))
    psQTL_fit = make_rho_fit(pp_sqtl_effects[0], pp_sqtl_effects[1], seed=seed,
                             x_grid=np.linspace(-3.2, 3.2, 100))

    n_donors = get_donor_count(FIBROBLASTS)
    usQTL_stats = make_rho_stats(u_sqtl_effects[0], u_sqtl_effects[1],
                                 n_donors=n_donors, seed=seed)
    psQTL_stats = make_rho_stats(pp_sqtl_effects[0], pp_sqtl_effects[1],
                                 n_donors=n_donors, seed=seed)

    usQTL_stats['selection'] = FIG4B_usQTL_SELECTION
    psQTL_stats['selection'] = FIG4B_psQTL_SELECTION
    for st in (usQTL_stats, psQTL_stats):
        st['unit'] = FIG4B_UNIT
        st['donor_unit'] = FIG4B_DONOR_UNIT
        st['tissue'] = FIBROBLASTS

    return usQTL_fit, psQTL_fit, usQTL_stats, psQTL_stats


# =========================================================================== #
# hyprcoloc tables -- fig4D and fig4_colocs_sc2   (hyprcoloc_results.ipynb)
# =========================================================================== #

def load_hyprcoloc(table=HYPRCOLOC_TABLE):
    """hyprcoloc results with the per-row GWAS / expression / splicing coloc flags."""
    hyprcoloc = pd.read_csv(table, sep='\t')

    gwas_colocs = []
    up_colocs = []
    expr_colocs = []
    lfc_colocs = []
    lfc_other_colocs = []
    tissue = []
    for idx, row in hyprcoloc.iterrows():
        if row.traits is np.nan:
            gwas_colocs.append(False)
            up_colocs.append(False)
            expr_colocs.append(False)
            lfc_colocs.append(False)
            lfc_other_colocs.append(False)
            tissue.append(np.nan)
        else:
            trait_list = row.traits.split(', ')
            gwas_coloc = row.gwas_trait in trait_list
            if gwas_coloc:
                gwas_colocs.append(True)
            else:
                gwas_colocs.append(False)
            if ':UP' in row.traits:
                up_colocs.append(True)
            else:
                up_colocs.append(False)

            if (':NE' in row.traits) or (':IN' in row.traits):
                lfc_other_colocs.append(True)
            else:
                lfc_other_colocs.append(False)

            if 'expression' in row.traits:
                expr_colocs.append(True)
            else:
                expr_colocs.append(False)
            if 'leafcutter' in row.traits:
                lfc_colocs.append(True)
                traits = row.traits.split(', ')
                for x in traits:
                    if 'leafcutter' in x:
                        tissue.append(x.split('.')[0])
                        break
            else:
                lfc_colocs.append(False)
                tissue.append(np.nan)

    hyprcoloc['gwas_colocs'] = gwas_colocs
    hyprcoloc['UP_colocs'] = up_colocs
    hyprcoloc['expr_colocs'] = expr_colocs
    hyprcoloc['lfc_colocs'] = lfc_colocs
    hyprcoloc['lfc_other_colocs'] = lfc_other_colocs

    hyprcoloc['gwas_loci'] = hyprcoloc.gwas_trait
    hyprcoloc['tissue'] = tissue
    hyprcoloc['gwas_trait'] = hyprcoloc.gwas_loci.apply(lambda x: x.split('_N_N_')[1])

    return hyprcoloc


def load_lead_snps(bed=LEAD_SNPS_BED):
    """GWAS lead-SNP windows -- the denominator of the fig4D percentages."""
    leadSNPs = pd.read_csv(bed, sep='\t', names=['chrom', 'start', 'end', 'gwas_loci'])
    leadSNPs['gwas_trait'] = leadSNPs.gwas_loci.apply(lambda x: x.split('_N_N_')[1])
    return leadSNPs


def make_fig4D_data(hyprcoloc, leadSNPs):
    """Per-trait percentage of GWAS loci colocalizing with a u-sQTL, by tissue."""
    hyprcoloc_up = hyprcoloc.loc[hyprcoloc.UP_colocs & hyprcoloc.gwas_colocs]

    gwas_pct_list = []
    gwas_pct_all_list = []
    trait_list = []
    for trait, df in hyprcoloc_up.groupby('gwas_trait'):
        gwas_pct = pd.DataFrame(
            100 * df.groupby(['tissue']).tissue.value_counts()
            / (leadSNPs.loc[leadSNPs.gwas_trait == trait].shape[0])
        )
        gwas_pct.columns = [trait]
        gwas_pct_list.append(gwas_pct)

        gwas_pct_all = (len(df.gwas_loci.unique())
                        / (leadSNPs.loc[leadSNPs.gwas_trait == trait].shape[0]))
        gwas_pct_all_list.append(gwas_pct_all)

        trait_list.append(trait)

    gwas_pct_df = pd.concat(gwas_pct_list, axis=1).fillna(0)
    gwas_pct_all_df = pd.DataFrame()
    gwas_pct_all_df['trait'] = trait_list
    gwas_pct_all_df['gwas_pct'] = gwas_pct_all_list

    ordered_traits = list(gwas_pct_all_df.sort_values('gwas_pct', ascending=False).trait)
    gwas_pct_all_df.index = gwas_pct_all_df.trait

    return {
        'boxes': gwas_pct_df[ordered_traits],
        'overall_pct': list(100 * gwas_pct_all_df.loc[ordered_traits].gwas_pct),
        'trait_names': FIG4D_TRAIT_NAMES,
    }


def make_fig4D_stats(hyprcoloc, leadSNPs, fig4D_data):
    """Box-by-box sample sizes for fig4D (Supplementary Table 10).

    Each box is one GWAS trait; each observation inside a box is one GTEx tissue.
    Because `gwas_pct_df` is built with `fillna(0)`, every box contains the same
    number of tissues -- tissues in which that trait had no u-sQTL colocalization
    contribute a structural zero rather than being dropped. The legend has to say
    so, and has to give the exact n.
    """
    boxes = fig4D_data['boxes']
    hyprcoloc_up = hyprcoloc.loc[hyprcoloc.UP_colocs & hyprcoloc.gwas_colocs]

    per_trait = []
    for trait in boxes.columns:
        per_trait.append({
            'trait': trait,
            'n_tissues_in_box': int(boxes.shape[0]),
            'n_tissues_nonzero': int((boxes[trait] > 0).sum()),
            'n_gwas_loci_tested': int(leadSNPs.loc[leadSNPs.gwas_trait == trait].shape[0]),
            'n_gwas_loci_colocalized': int(
                hyprcoloc_up.loc[hyprcoloc_up.gwas_trait == trait].gwas_loci.nunique()),
            'median_pct': float(boxes[trait].median()),
            'q1_pct': float(boxes[trait].quantile(0.25)),
            'q3_pct': float(boxes[trait].quantile(0.75)),
            'min_pct': float(boxes[trait].min()),
            'max_pct': float(boxes[trait].max()),
        })

    # Traits that were tested but produced no u-sQTL colocalization in any tissue
    # get no box in the panel. They are still part of the analysis and belong in
    # Source Data as an explicit zero, flagged with in_panel = False.
    for trait in sorted(set(leadSNPs.gwas_trait.unique()) - set(boxes.columns)):
        per_trait.append({
            'trait': trait,
            'n_tissues_in_box': int(boxes.shape[0]),
            'n_tissues_nonzero': 0,
            'n_gwas_loci_tested': int(leadSNPs.loc[leadSNPs.gwas_trait == trait].shape[0]),
            'n_gwas_loci_colocalized': int(
                hyprcoloc_up.loc[hyprcoloc_up.gwas_trait == trait].gwas_loci.nunique()),
            'median_pct': 0.0, 'q1_pct': 0.0, 'q3_pct': 0.0,
            'min_pct': 0.0, 'max_pct': 0.0,
        })

    per_trait = pd.DataFrame(per_trait)
    per_trait['in_panel'] = per_trait.trait.isin(boxes.columns)

    # 'trait' stays the hyprcoloc key; the human-readable name becomes trait_name
    meta = pd.DataFrame([{'trait': key,
                          'trait_name': v['trait'],
                          **{k: v[k] for k in v if k != 'trait'}}
                         for key, v in FIG4D_TRAIT_ACCESSIONS.items()])
    per_trait = per_trait.merge(meta, on='trait', how='left')
    missing = per_trait.loc[per_trait.accession.isna(), 'trait'].tolist()
    if missing:
        raise ValueError(f'no accession recorded for: {missing}')

    return {
        'n_boxes': int(boxes.shape[1]),
        'n_traits_tested': int(leadSNPs.gwas_trait.nunique()),
        'n_per_box': int(boxes.shape[0]),
        'n_gwas_loci_total': int(leadSNPs.gwas_loci.nunique()),
        'per_trait': per_trait,
        'unresolved_accessions': FIG4D_UNRESOLVED_ACCESSIONS,
        'box_definition': (
            'seaborn/matplotlib boxplot defaults: centre line, median; box bounds, '
            '25th and 75th percentiles (Q1, Q3); whiskers, most extreme value within '
            '1.5 x IQR of the nearer hinge; points beyond the whiskers are plotted '
            'individually as outliers. Minima and maxima of each box are given in the '
            'per_trait table.'),
    }


def load_up_results(table=UP_RESULTS_TABLE):
    """Per-junction GWAS-colocalizing sQTL / eQTL betas (hyprcoloc_results.ipynb)."""
    return pd.read_csv(table, sep='\t')


def make_coloc_pr_fit(nom_df):
    """PR-intron scatter fit behind hyprcoloc_results.ipynb's fig4_colocs_sc2."""
    x = np.array(nom_df.loc[(nom_df.distance_abs <= 1000) & (nom_df.itype == 'PR')
                            & (nom_df.PSI >= 0)].sQTL_beta)
    y = np.array(nom_df.loc[(nom_df.distance_abs <= 1000) & (nom_df.itype == 'PR')
                            & (nom_df.PSI >= 0)].eQTL_beta)
    return make_rho_fit(x, y, x_grid=np.linspace(-2, 2, 100))


# =========================================================================== #
# fig4E / fig4F -- the ASB16 example locus     (coloc_plots.ipynb,
#                                               Fig4_example.ipynb)
# =========================================================================== #

def run_tabix(tabix_file, chrom, start, end, file_type='other'):
    """Rows of a bgzipped/tabixed table overlapping a region, with its own header."""
    import tabix

    tabix_obj = tabix.open(tabix_file)
    result = tabix_obj.query(chrom, start, end)

    rows = [line for line in result]

    with gzip.open(tabix_file) as fh:
        cols = fh.readline().decode().rstrip().split('\t')

    df = pd.DataFrame(rows, columns=cols)
    if file_type == 'leafcutter':
        df['gid'] = df.pid.apply(lambda x: x.split(':')[3])
    elif file_type == 'expression':
        df['gid'] = df.pid.apply(lambda x: x.split('.')[0])

    return df


def get_best_snp_record(vcf_reader, snp):
    """VCF record(s) at 'chrom:pos'; the exact-position record when there is one."""
    chrom, location = snp.split(':')
    location = int(location)

    record_list = []
    for record in vcf_reader.fetch(chrom, location - 1, location + 1):
        if record.POS == location:
            return [record]
        record_list.append(record)

    return record_list


def make_genotype_boxplot_df(qqnorm_data, phenotype, record):
    """Genotype / normalized-phenotype table for one variant-phenotype pair.

    Samples are ordered hom-ref, het, hom-alt, exactly as in coloc_plots.ipynb's
    `make_boxplot_df_colored`.
    """
    qq_samples = qqnorm_data.columns[6:]
    hom_refs = [str(x.sample) for x in record.get_hom_refs() if str(x.sample) in qq_samples]
    hets = [str(x.sample) for x in record.get_hets() if str(x.sample) in qq_samples]
    hom_alts = [str(x.sample) for x in record.get_hom_alts() if str(x.sample) in qq_samples]

    samples = hom_refs + hets + hom_alts

    ref = record.REF
    alt = record.ALT[0]

    genotype = [str(ref) + ',' + str(ref)] * len(hom_refs)
    genotype += [str(ref) + ',' + str(alt)] * len(hets)
    genotype += [str(alt) + ',' + str(alt)] * len(hom_alts)

    df = pd.DataFrame()
    df['genotype'] = genotype
    df['qqnorm'] = list(qqnorm_data.loc[phenotype, samples].astype(float))
    df.index = hom_refs + hets + hom_alts

    return df


def get_asb16_qtl_stats(n_donors):
    """Exact QTLtools nominal-pass statistics for the two fig4E panels.

    Returns one dict per panel with the effect size (beta), its standard error,
    the exact two-sided nominal P value, r^2, degrees of freedom and the 95% CI on
    beta. QTLtools tests the association between the normalized phenotype and the
    genotype dosage under an additive linear model, two-sided, after regressing out
    the covariates; the P value is nominal, i.e. not adjusted for the number of
    variants in cis (that adjustment is what the permutation pass does, and it is
    the q <= 0.1 threshold used to select these loci in the first place).
    """
    from scipy.stats import t as t_dist

    def _one(qtl_type, phenotype):
        nom = run_tabix(
            f'{ASB16_NOM_DIR}/{ASB16_TISSUE}/{qtl_type}.NominalPass/{ASB16_CHROM}.txt.tabix.gz',
            ASB16_CHROM, ASB16_START - 30000, ASB16_END + 30000)
        pos = ASB16_SNP.split(':')[1]
        row = nom.loc[(nom['#phe_id'] == phenotype) & (nom.var_from == pos)].iloc[0]

        beta = float(row.slope)
        se = float(row.slope_se)
        df = n_donors - 2
        t_crit = t_dist.ppf(0.975, df)

        # Phenotype-level permutation-adjusted P (beta approximation) and Storey q,
        # from the same coloc QTL run. Reported in the legend alongside the nominal
        # P; the panel annotates the nominal beta and P, which describe the plotted
        # variant, while adj. P is a phenotype-level statistic that QTLtools earns
        # at that phenotype's own top variant.
        perm = pd.read_csv(
            f'{ASB16_NOM_DIR}/{ASB16_TISSUE}/{qtl_type}.PermutationPass.FDR_Added.txt.gz',
            sep=' ')
        prow = perm.loc[perm.phe_id == phenotype].iloc[0]

        return {
            'phenotype': phenotype,
            'variant': row.var_id,
            'beta': beta,
            'beta_se': se,
            'pvalue': float(row.nom_pval),
            'r_squared': float(row.r_squared),
            'df': int(df),
            't': beta / se,
            'ci95': (beta - t_crit * se, beta + t_crit * se),
            'adj_pvalue': float(prow.adj_beta_pval),
            'qvalue': float(prow.q),
            'perm_top_variant': str(prow.var_id),
            'test': ('two-sided t test on the slope of an additive linear model '
                     '(QTLtools nominal pass); nominal P, not adjusted for the '
                     'number of variants in cis. adj_pvalue is the '
                     'permutation-adjusted P (beta approximation) for the phenotype.'),
        }

    return (_one('leafcutter', ASB16_SPLICE_PHE), _one('expression', ASB16_EXPR_PHE))


def make_asb16_boxplot_stats(df):
    """Per-genotype n and box statistics for one fig4E genotype boxplot."""
    groups = []
    for genotype, sub in df.groupby('genotype', sort=False):
        q1, med, q3 = sub.qqnorm.quantile([0.25, 0.5, 0.75])
        iqr = q3 - q1
        inside = sub.qqnorm[(sub.qqnorm >= q1 - 1.5 * iqr) & (sub.qqnorm <= q3 + 1.5 * iqr)]
        groups.append({
            'genotype': genotype,
            'n': int(len(sub)),
            'median': float(med),
            'q1': float(q1),
            'q3': float(q3),
            'whisker_low': float(inside.min()),
            'whisker_high': float(inside.max()),
            'min': float(sub.qqnorm.min()),
            'max': float(sub.qqnorm.max()),
        })
    return {
        'n_total': int(len(df)),
        'groups': pd.DataFrame(groups),
        'box_definition': (
            'centre line, median; box bounds, 25th and 75th percentiles (Q1, Q3); '
            'whiskers, most extreme value within 1.5 x IQR of the nearer hinge; '
            'outliers are not drawn as separate markers because every donor is '
            'overlaid as an individual point (stripplot, showfliers=False).'),
    }


def make_fig4F_data():
    """LocusCompare points for fig4F: GWAS, eQTL and u-sQTL P values at ASB16.

    Follows Fig4_example.ipynb, which intersects only the four tables actually
    plotted here. coloc_plots.ipynb builds the same panel but intersects several
    further tissues as well, which shrinks the SNP set -- do not mix the two.

    Returns -log10 P values already taken, so the plot cell only draws.
    """
    win_lo = ASB16_START - 1000000
    win_hi = ASB16_END + 1000000

    gwas = run_tabix(BIPOLAR_STATS, ASB16_CHROM, win_lo, win_hi)
    gwas['snp_pos'] = gwas.chrom + ':' + gwas.start

    gwas_pval = run_tabix(BIPOLAR_PVALS, ASB16_CHROM, win_lo, win_hi)
    gwas_pval['snp_pos'] = gwas_pval['chrom'] + ':' + gwas_pval.start

    nom_bch = run_tabix(
        f'{ASB16_NOM_DIR}/{ASB16_TISSUE}/leafcutter.NominalPass/{ASB16_CHROM}.txt.tabix.gz',
        ASB16_CHROM, win_lo, win_hi)
    nom_bch['snp_pos'] = nom_bch.var_chr + ':' + nom_bch.var_from

    nom_bch_e = run_tabix(
        f'{ASB16_NOM_DIR}/{ASB16_TISSUE}/expression.NominalPass/{ASB16_CHROM}.txt.tabix.gz',
        ASB16_CHROM, win_lo, win_hi)
    nom_bch_e['snp_pos'] = nom_bch_e.var_chr + ':' + nom_bch_e.var_from

    for d in (nom_bch, nom_bch_e):
        d['gwas_trait'] = d['#phe_id'].apply(lambda x: x.split(':')[-1])

    nom_bch = nom_bch.loc[(nom_bch.gwas_trait == ASB16_GWAS_TRAIT)
                          & (nom_bch['#phe_id'] == ASB16_SPLICE_PHE)]
    nom_bch_e = nom_bch_e.loc[
        (nom_bch_e.gwas_trait == ASB16_GWAS_TRAIT)
        & (nom_bch_e['#phe_id'].apply(lambda x: (ASB16_GENE in x)
                                      and (ASB16_GWAS_TRAIT in x)))]

    snps = (pd.Index(nom_bch.snp_pos)
            .intersection(pd.Index(gwas_pval.snp_pos))
            .intersection(pd.Index(gwas.snp_pos))
            .intersection(pd.Index(nom_bch_e.snp_pos)))

    gwas_pval = gwas_pval.groupby('snp_pos').first()
    nom_bch = nom_bch.groupby('snp_pos').first()
    nom_bch_e = nom_bch_e.groupby('snp_pos').first()

    LD = pd.read_csv(ASB16_LD, sep='\t')
    LD.columns = [':'.join(x.split('_')[:2]) for x in LD.columns]
    LD.index = LD.columns
    LD = LD.astype(float) ** 2

    snps_shared = snps.intersection(LD.columns)

    points = pd.DataFrame(index=snps_shared)
    points['GWAS'] = -np.log10(np.array(gwas_pval.loc[snps_shared].P.astype(float)))
    points['sQTL'] = -np.log10(np.array(nom_bch.loc[snps_shared].nom_pval.astype(float)))
    points['eQTL'] = -np.log10(np.array(nom_bch_e.loc[snps_shared].nom_pval.astype(float)))
    points['r2'] = np.array(LD.loc[ASB16_SNP, snps_shared])

    lead = {
        'label': ASB16_LEAD_LABEL,
        'GWAS': float(-np.log10(float(gwas_pval.loc[ASB16_SNP].P))),
        'sQTL': float(-np.log10(float(nom_bch.loc[ASB16_SNP].nom_pval))),
        'eQTL': float(-np.log10(float(nom_bch_e.loc[ASB16_SNP].nom_pval))),
    }

    return {
        'points': points,                            # panels 2 and 3 draw in this order
        'points_sorted': points.sort_values('r2'),   # panel 1 draws sorted by r2
        'lead': lead,
        'n_snps': int(len(snps_shared)),
    }


def make_asb16_data():
    """Boxplot tables and exon annotation behind the three fig4E panels."""
    import vcf

    window_start = ASB16_START - 30000
    window_end = ASB16_END + 30000

    qqnorm_bch = run_tabix(
        f'{COLOC_DATA_DIR}/{ASB16_TISSUE}/{ASB16_CHROM}.leafcutter.ForGWASColoc.sorted.qqnorm.bed.gz',
        ASB16_CHROM, window_start, window_end, file_type='leafcutter')
    qqnorm_bch_e = run_tabix(
        f'{COLOC_DATA_DIR}/{ASB16_TISSUE}/{ASB16_CHROM}.expression.ForGWASColoc.sorted.qqnorm.bed.gz',
        ASB16_CHROM, window_start, window_end, file_type='expression')

    qqnorm_bch.index = qqnorm_bch.pid
    qqnorm_bch_e.index = qqnorm_bch_e.pid

    vcf_gtex = vcf.Reader(open(GTEX_VCF, 'br'))
    records = get_best_snp_record(vcf_gtex, ASB16_SNP)

    sqtl_boxplot_df = make_genotype_boxplot_df(qqnorm_bch, ASB16_SPLICE_PHE, records[0])
    eqtl_boxplot_df = make_genotype_boxplot_df(qqnorm_bch_e, ASB16_EXPR_PHE, records[0])

    annot = run_tabix(GENCODE_EXONS, ASB16_CHROM, window_start, window_end)
    asb16_annot = annot.loc[annot.transcript_id.isin(ASB16_TRANSCRIPTS)]

    return sqtl_boxplot_df, eqtl_boxplot_df, asb16_annot


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #

PLOT_READY_VARS = [
    'fig4A_data',
    'fig4B_usQTL_fit', 'fig4B_psQTL_fit',
    'fig4B_usQTL_stats', 'fig4B_psQTL_stats',
    'fig4C_data', 'fig4C_source_data',
    'fig4D_data', 'fig4D_stats',
    'fig4E_sQTL_boxplot_df', 'fig4E_eQTL_boxplot_df', 'fig4E_ASB16_annot',
    'fig4E_sQTL_stats', 'fig4E_eQTL_stats',
    'fig4E_sQTL_box_stats', 'fig4E_eQTL_box_stats',
    'fig4F_data',
    'fig4_colocs_sc2_fit',
]


def run_all(data_dir='figure_data'):
    """Run the full pipeline, pickle every plot-ready variable, and return them."""
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(PLOTS_DIR, exist_ok=True)

    donor_counts = get_all_donor_counts()

    fig4A_data = make_fig4A_data(donor_counts)
    fig4C_data = make_fig4C_data()
    fig4C_source_data = make_fig4C_source_data(donor_counts)
    (fig4B_usQTL_fit, fig4B_psQTL_fit,
     fig4B_usQTL_stats, fig4B_psQTL_stats) = make_fibroblast_fits()

    hyprcoloc = load_hyprcoloc()
    leadSNPs = load_lead_snps()
    fig4D_data = make_fig4D_data(hyprcoloc, leadSNPs)
    fig4D_stats = make_fig4D_stats(hyprcoloc, leadSNPs, fig4D_data)
    fig4_colocs_sc2_fit = make_coloc_pr_fit(load_up_results())

    (fig4E_sQTL_boxplot_df, fig4E_eQTL_boxplot_df,
     fig4E_ASB16_annot) = make_asb16_data()
    n_donors = len(fig4E_sQTL_boxplot_df)
    fig4E_sQTL_stats, fig4E_eQTL_stats = get_asb16_qtl_stats(n_donors)
    fig4E_sQTL_box_stats = make_asb16_boxplot_stats(fig4E_sQTL_boxplot_df)
    fig4E_eQTL_box_stats = make_asb16_boxplot_stats(fig4E_eQTL_boxplot_df)

    fig4F_data = make_fig4F_data()

    # Source-data tables the legends point to
    fig4C_source_data.to_csv(os.path.join(data_dir, 'fig4C_source_data.tsv'),
                             sep='\t', index=False)
    fig4D_stats['per_trait'].to_csv(os.path.join(data_dir, 'fig4D_source_data.tsv'),
                                    sep='\t', index=False)

    data = {
        'fig4A_data': fig4A_data,
        'fig4B_usQTL_fit': fig4B_usQTL_fit,
        'fig4B_psQTL_fit': fig4B_psQTL_fit,
        'fig4B_usQTL_stats': fig4B_usQTL_stats,
        'fig4B_psQTL_stats': fig4B_psQTL_stats,
        'fig4C_data': fig4C_data,
        'fig4C_source_data': fig4C_source_data,
        'fig4D_data': fig4D_data,
        'fig4D_stats': fig4D_stats,
        'fig4E_sQTL_boxplot_df': fig4E_sQTL_boxplot_df,
        'fig4E_eQTL_boxplot_df': fig4E_eQTL_boxplot_df,
        'fig4E_ASB16_annot': fig4E_ASB16_annot,
        'fig4E_sQTL_stats': fig4E_sQTL_stats,
        'fig4E_eQTL_stats': fig4E_eQTL_stats,
        'fig4E_sQTL_box_stats': fig4E_sQTL_box_stats,
        'fig4E_eQTL_box_stats': fig4E_eQTL_box_stats,
        'fig4F_data': fig4F_data,
        'fig4_colocs_sc2_fit': fig4_colocs_sc2_fit,
    }

    for name, value in data.items():
        with open(os.path.join(data_dir, f'{name}.pickle'), 'wb') as fh:
            pickle.dump(value, fh)

    return data


def load_plot_data(data_dir='figure_data'):
    """Load every plot-ready variable back from `data_dir` (no recomputation)."""
    data = {}
    for name in PLOT_READY_VARS:
        with open(os.path.join(data_dir, f'{name}.pickle'), 'rb') as fh:
            data[name] = pickle.load(fh)
    return data
