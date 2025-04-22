suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(BEDMatrix))
suppressMessages(library(glue))
suppressMessages(library(bedr))
suppressMessages(library(readr))

suppressMessages(library(vcfR))
suppressMessages(library(genetics))
suppressMessages(library(hyprcoloc))

library(tidyverse)
library(dplyr)
library(bedr)
library(readr)
library(glue)
library(qvalue)

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(data.table))
suppressMessages(library(bedtoolsr))
suppressMessages(library(GenomicRanges))

for(file in list.files("colocboost/v10/", #"LeafCutter2_ColocBoost_GTEx/final_checking", 
                       pattern = ".R", full.names = T)){ source(file) }
source("LeafCutter2_ColocBoost_GTEx/preprocess_functions.R")



getPermutationList <- function(tissue_list, pheno){
  Perm_list <- lapply(setNames(tissues, tissues), function(tissue) {
  file_path <- str_glue("results/coloc/qtls/{tissue}/{pheno}.PermutationPass.FDR_Added.txt.gz")
  
  # Read and filter the data
  if (pheno == 'leafcutter'){
      perm <- read_delim(file_path, delim = " ") %>%
        filter(q <= 0.1) %>%
        mutate(
          gwas_loci = str_split(phe_id, ":", simplify = TRUE)[, 6], 
          gwas_trait = str_split(phe_id, "_N_N_", simplify = TRUE)[, 2],
          intron = apply(str_split(phe_id, ':', simplify=TRUE)[,1:5], 1, paste, collapse = ":")
        )
    } else if (pheno == 'expression') {
    perm <- read_delim(file_path, delim = " ") %>%
    filter(q <= 0.1) %>%
    mutate(
      gwas_loci = str_split(phe_id, ":", simplify = TRUE)[, 2], 
      gwas_trait = str_split(phe_id, "_N_N_", simplify = TRUE)[, 2],
      gene = str_split(phe_id, '\\.', simplify=TRUE)[,1]
    )
  }
  
  return(perm)
})
}


search_tissues <- function(gwas_loci_id, data_list, p=0.1) {
  matching_tissues <- names(data_list)[sapply(data_list, function(df) {
    any(df$gwas_loci == gwas_loci_id & df$q <= p)
  })]
  
  # Create a list where the key is the tissue name and value is the corresponding dof1
  tissue_dof1_list <- lapply(matching_tissues, function(tissue) {
    # Find the row(s) where gwas_loci == gwas_loci_id and q <= p
    matching_rows <- data_list[[tissue]][data_list[[tissue]]$gwas_loci == gwas_loci_id & data_list[[tissue]]$q <= p, ]
    
    # Return the dof1 values (assuming dof1 exists and should be numeric)
    return(as.numeric(matching_rows$dof1[1]))
  })
  
  names(tissue_dof1_list) <- matching_tissues
  return(tissue_dof1_list)
}

getSummaryStatsPheno <- function(gwas_loci, perm_list, perm_name, p, gtf){
    # Get a list of tissues and their corresponding dof1 values
    select_tissues <- search_tissues(gwas_loci, perm_list, p)
    
    # If no tissues match, return NULL
    if (length(select_tissues) == 0) {
        return(NULL)
    } else {
        combined_list <- do.call(c, lapply(names(select_tissues), function(tissue) {
            dof1 <- select_tissues[[tissue]]
            print(perm_name)
            getTissueRegionStats(region, gwas_loci, tissue, dof1, pheno=perm_name, gtf=gtf)
        }))
        
        return(combined_list)
    }
}


getColocSummaryStats <- function(region, gwas_loci, leafcutter_list, expression_list, gtf, p=0.1) {
    
    leafcutter_stats <- getSummaryStatsPheno(gwas_loci, leafcutter_list, 'leafcutter', p, gtf=gtf)
    expression_stats <- getSummaryStatsPheno(gwas_loci, expression_list, 'expression', p, gtf=gtf)
    
    if (is.null(leafcutter_stats) & is.null(expression_stats)) {return (NULL)}

    # Get the GWAS stats
    gwas_name <- getGWASName(gwas_loci)
    gwas_stats <- getGWASstats(region, gwas_loci, gwas_name)                                   

    # Combine the GWAS stats with the results from tissues
    combined_list <- c(gwas_stats, leafcutter_stats, expression_stats)                                   
    return(combined_list)
    
}

getGWASName <- function(gwas_loci){
    gwas_loci_name <- (gwas_loci %>% str_split(., '_N_N_'))[[1]][2]
    if (gwas_loci_name == 'GCST004131') {return('Inflammatory_bowel_disease')}
    else if (gwas_loci_name == 'GCST004132') {return('Crohns_disease')}
    else if (gwas_loci_name == 'GCST004133') {return('Ulcerative_colitis')}
    else if (gwas_loci_name == 'GCST004988') {return('Breast_cancer')}
    else if (gwas_loci_name == 'GCST007800') {return('Asthma_childhood_onset')}
    else if (gwas_loci_name == 'age_when_finished_full-time_education') {return('Age_when_finished_full-time_education')}
    else if (gwas_loci_name == 'coronary_artery_disease') {return('Coronary_artery_disease')}
    else if (gwas_loci_name == 'atrial_fibrillation') {return('Atrial_fibrillation')}
    else if (gwas_loci_name == 'myocardial_infarction') {return('Myocardial_infarction')}
    else if (gwas_loci_name == 'heart_failure') {return('Heart_failure')}
    else if (gwas_loci_name == 'basal_cell_carcinoma') {return('Basal_cell_carcinoma')}
    else if (gwas_loci_name == 'bipolar_disorder') {return('Bipolar_disorder')}
    else if (gwas_loci_name == 'schizophrenia') {return('Schizophrenia')}
    else if (gwas_loci_name == 'IMSGC2019') {return('Multiple_sclerosis')}
    else {return(gwas_loci_name)}
    
}

getGWASn <- function(gwas_name){
    if (gwas_name == 'Chronic_obstructive_pulmonary_disease') {return(325027)}
    else if (gwas_name == 'Inflammatory_bowel_disease') {return(59957)}
    else if (gwas_name == 'Crohns_disease') {return(40266)}
    else if (gwas_name == 'Ulcerative_colitis') {return(45975)}
    else if (gwas_name == 'Breast_cancer') {return(139274)}
    else if (gwas_name == 'Asthma_childhood_onset') {return(314633)}
    else if (gwas_name == 'Basal_cell_carcinoma') {return(307684)}
    else if (gwas_name == 'Multiple_sclerosis') {return(115803)}
    else if (gwas_name == 'Age_when_finished_full-time_education') {return(283749)}
    else if (gwas_name == 'Coronary_artery_disease') {return(1165690)}
    else if (gwas_name == 'Atrial_fibrillation') {return(588190)}
    else if (gwas_name == 'Myocardial_infarction') {return(639221)}
    else if (gwas_name == 'Heart_failure') {return(1665481)}
    else if (gwas_name == 'Bipolar_disorder') {return(413466)}
    else if (gwas_name == 'Schizophrenia') {return(130644)}
    else if (gwas_name == 'Rheumatoid_arthritis') {return(80799)}
    else if (gwas_name == 'Dupuytrens_disease') {return(58343)}
    else if (gwas_name == 'Hypothyroidism') {return(691986)}
    else if (gwas_name == 'Visceral_adipose_tissue_measurement') {return(325153)}
    else if (gwas_name == 'Atopic_eczema') {return(864982)}
}

getGWASstats <- function(region, gwas_loci_id, gwas_name){
    gwas.file <- glue('resources/gwas/StatsForColoc/{gwas_name}.standardized.txt.tabix.gz')
    tmpDir_gwas <- glue('tmp/hyprcoloc/{gwas_name}/{gwas_loci_id}/')
    dir.create(tmpDir_gwas, recursive=TRUE)
    ibs <- tabix(region, gwas.file, tmpDir=tmpDir_gwas)
    colnames(ibs) <- c('gwas_loci', 'chrom', 'pos', 'beta', 'sebeta', 'A1', 'A2')
    ibs <- ibs %>% filter(gwas_loci == gwas_loci_id) %>%
        mutate(variant = str_c(chrom, pos, A2, A1, sep = "_", "b38")) %>%
        dplyr::select(c(beta, sebeta, variant))
    
    ibs$beta <- as.numeric(ibs$beta)
    ibs$sebeta <- as.numeric(ibs$sebeta)
    gwas_n <- getGWASn(gwas_name)
    ibs$n <- gwas_n
    
    ibs <- ibs %>% dplyr::select(c(beta, sebeta, n, variant)) %>% filter(sebeta > 0) %>% na.omit()
    ibs <- ibs[!duplicated(ibs$variant), ]

    if (gwas_name %in% c('Schizophrenia', 'Bipolar_disorder', 'Basal_cell_carcinoma', 
                         'Myocardial_infarction', 'Visceral_adipose_tissue_measurement')) {
        ibs$variant <- ibs$variant %>%  gsub("(chr\\d+_\\d+)_([A-Z])_([A-Z])_(b\\d+)", "\\1_\\3_\\2_\\4", .)
    } else if (gwas_name == 'Coronary_artery_disease') {
        ibs$variant <- ibs$variant %>%  gsub(
          "(chr\\d+_\\d+)_([a-z])_([a-z])_(b\\d+)", 
          "\\1_\\U\\2_\\U\\3_b38", 
          ., 
          perl = TRUE
        )
    }
    
    print(dim(ibs))
    
    ibs <- list(ibs)
    names(ibs) <- c(gwas_loci_id)
    return(ibs)

}

getTissueRegionSNPs <- function(region, gwas_loci_id, tissue, pheno){
    chrom <- str_split(region, ':', simplify=TRUE)[1,][1]
    tmpDir <- glue('tmp/hyprcoloc/{tissue}/{pheno}/{chrom}/{gwas_loci}/')
    dir.create(tmpDir, recursive=TRUE)
    X <- tabix(region, glue('results/coloc/qtls/{tissue}/{pheno}.NominalPass/{chrom}.txt.tabix.gz'), tmpDir=tmpDir)
    if (pheno == 'leafcutter'){
        X <- X %>%
            mutate(gwas_loci = str_split(phe_id, ':', simplify=TRUE)[,6], 
                     gwas_trait = str_split(phe_id, "_N_N_", simplify = TRUE)[, 2],
                  intron = apply(str_split(phe_id, ':', simplify=TRUE)[,1:5], 1, paste, collapse = ":")) %>%
            filter(gwas_loci == gwas_loci_id) %>% 
            dplyr::select(c(slope, slope_se, var_id, intron))
        colnames(X) <- c('beta', 'sebeta', 'variant', 'intron')
    } else if (pheno == 'expression') {
        X <- X %>%
        mutate(gwas_loci = str_split(phe_id, ':', simplify=TRUE)[,2], 
                 gwas_trait = str_split(phe_id, "_N_N_", simplify = TRUE)[, 2],
              gene = str_split(phe_id, '\\.', simplify=TRUE)[,1]) %>%
        filter(gwas_loci == gwas_loci_id) %>% 
        dplyr::select(c(slope, slope_se, var_id, gene))
        colnames(X) <- c('beta', 'sebeta', 'variant', 'gene')
    }
    
    X$beta <- as.numeric(X$beta)
    X$sebeta <- as.numeric(X$sebeta)
    
    X <- X %>% filter(sebeta > 0) %>% na.omit()
    X <- X[!duplicated(X$variant), ]
    print(dim(X))
    
    return (X)
    
}

getTissueRegionStats <- function(region, gwas_loci_id, tissue, dof, pheno, gtf) {
    X_stats <- getTissueRegionSNPs(region, gwas_loci_id, tissue, pheno)
    n <- dof+2
    X_stats$n <- n
    
    if (pheno == 'leafcutter'){
    
    X_list <- X_stats %>%
      group_by(intron) %>%
      group_split() %>% lapply(., function(df) {
      df %>% dplyr::select(beta, sebeta, n, variant) %>% as.data.frame()
    })
    pheno_name <- paste(tissue, pheno, sep='.')
    names(X_list) <- paste(pheno_name, unique(X_stats$intron), sep='_N_N_') %>% sapply(., annotate.intron, gtf = gtf)
    } else if (pheno == 'expression') {
    
    X_list <- X_stats %>%
      group_by(gene) %>%
      group_split() %>% lapply(., function(df) {
      df %>% dplyr::select(beta, sebeta, n, variant) %>% as.data.frame()
    })
    pheno_name <- paste(tissue, pheno, sep='.')
    names(X_list) <- paste(pheno_name, unique(X_stats$gene), sep='_N_N_') 
    }
    
    
    return(X_list)
}

getLDMatrix <- function(region){
    chrom <- (region %>% str_split(., pattern=':'))[[1]][1]
    vcf_file <- "/project/yangili1/cdai/genome_index/hs38/GTEx_v7/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz"  
    tmpDir_ld <- glue('tmp/hyprcoloc/{region}/')
    dir.create(tmpDir_ld, recursive=TRUE)
    gt_region <- bedr::tabix(region = region, file = vcf_file, tmpDir=tmpDir_ld) %>% 
        filter(FORMAT == 'GT')

    rownames(gt_region) <- gt_region$ID
    
    gt_filtered <- gt_region %>%
      mutate(AF = as.numeric(str_extract(INFO, "(?<=AF=)[0-9.]+"))) %>%
      filter(AF >= 0.05)

    # Select columns starting with "GTEX-"
    gt_region <- gt_filtered[, grepl("^GTEX-", colnames(gt_filtered))]

    gt_region <- gt_region[((gt_region == '.|.') %>% rowSums()) <= 838/20,]
    
    gt_region[gt_region == "0|0"] <- 0
    gt_region[(gt_region == "1|0") | (gt_region == "0|1")] <- 1
    gt_region[gt_region == "1|1"] <- 2
    gt_region[gt_region == ".|."] <- NA
    
    
#     return(gt_region)

    matrix_data_int <- gt_region %>% mutate_all(as.integer) #suppressWarnings(as.integer(gt_region))
    

    # Reshape into the original matrix structure
#     matrix_data_int <- matrix(matrix_data_int, nrow = nrow(gt_region), dimnames = dimnames(gt_region))
    

    imputed_data <- apply(matrix_data_int, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
                          
    X <- imputed_data %>% t() %>% cor()

    row_variances <- apply(imputed_data, 1, var)

    # Filter rows with non-zero variance
    filtered_matrix <- imputed_data[row_variances != 0, , drop = FALSE]
                          

    LD <- filtered_matrix %>% t() %>% cor()
    return (LD)
                          
}


                          
getRelevantTissues <- function(gwas.trait){
if (gwas.trait == 'Inflammatory_bowel_disease') { 
    
    tissues <- c('Adipose-Visceral_Omentum', 'Colon-Sigmoid', 'Colon-Transverse', 'Esophagus-GastroesophagealJunction','Esophagus-Mucosa',
                 'Esophagus-Muscularis', 'Spleen', 'WholeBlood', 'AdrenalGland', 'Pituitary', 'SmallIntestine-TerminalIleum', 'Stomach')
    }

else if (gwas.trait == 'Crohns_disease') { 
    tissues <- c('Adipose-Visceral_Omentum', 'Colon-Sigmoid', 'Colon-Transverse', 'Esophagus-GastroesophagealJunction','Esophagus-Mucosa',
                 'Esophagus-Muscularis', 'Spleen', 'WholeBlood', 'AdrenalGland', 'Pituitary')
}

else if (gwas.trait == 'Ulcerative_colitis') { 
    tissues <- c('Adipose-Visceral_Omentum', 'Colon-Sigmoid', 'Colon-Transverse', 'Esophagus-Mucosa',
                 'Spleen', 'WholeBlood', 'Stomach')}

else if (gwas.trait == 'Breast_cancer') { tissues <- c('Breast-MammaryTissue', 'Ovary', 'Pituitary', 'AdrenalGland', 
                                                       'WholeBlood', 'Spleen')}

else if (gwas.trait == 'Asthma_childhood_onset') { tissues <- c('Lung', 'Nerve-Tibial', 'MinorSalivaryGland', 
                                                                'WholeBlood', 'Spleen', 'Cells-EBV-transformedlymphocytes',
                                                               'AdrenalGland', 'Pituitary')}

else if (gwas.trait == 'Basal_cell_carcinoma') { tissues <- c('Skin-NotSunExposed_Suprapubic', 'Skin-SunExposed_Lowerleg', 
                                                              'WholeBlood', 'Spleen')}


else if (gwas.trait == 'Multiple_sclerosis') { tissues <- c('Brain-Anteriorcingulatecortex_BA24', 'Brain-Cerebellum', 
                                                            'Brain-Cortex', 'Brain-FrontalCortex_BA9', 
                        'Brain-Spinalcord_cervicalc-1', 'WholeBlood', 'Cells-EBV-transformedlymphocytes', 'Muscle-Skeletal', 
                        'Adipose-Subcutaneous', 'Spleen', 'Thyroid')}

else if (gwas.trait == 'Age_when_finished_full-time_education') { tissues <- c('Brain-Anteriorcingulatecortex_BA24', 
                                                                              'Brain-Cerebellum', 'Brain-Cortex', 
                                          'Brain-FrontalCortex_BA9', 
                        'Brain-Spinalcord_cervicalc-1', 'WholeBlood', 'Cells-EBV-transformedlymphocytes')}

else if (gwas.trait == 'Coronary_artery_disease') { tissues <- c('Artery-Aorta','Artery-Coronary','Artery-Tibial', 
                                                                 'Adipose-Subcutaneous',
                             'Adipose-Visceral_Omentum', 'Liver', 'Heart-AtrialAppendage', 'Heart-LeftVentricle', 'Kidney-Cortex',
                            'WholeBlood')}

else if (gwas.trait == 'Atrial_fibrillation') { tissues <- c('Artery-Aorta','Artery-Coronary','Artery-Tibial', 
                                                             'Adipose-Subcutaneous',
                             'Adipose-Visceral_Omentum', 'Liver', 'Heart-AtrialAppendage', 'Heart-LeftVentricle', 'Kidney-Cortex',
                            'WholeBlood', 'Muscle-Skeletal')}

else if (gwas.trait == 'Myocardial_infarction') { tissues <- c('Artery-Aorta','Artery-Coronary','Artery-Tibial', 
                                                               'Adipose-Subcutaneous',
                             'Adipose-Visceral_Omentum', 'Liver', 'Heart-AtrialAppendage', 'Heart-LeftVentricle', 'Kidney-Cortex',
                            'WholeBlood', 'Muscle-Skeletal')}

else if (gwas.trait == 'Heart_failure') { tissues <- c('Artery-Aorta','Artery-Coronary','Artery-Tibial', 
                                                       'Adipose-Subcutaneous',
                             'Adipose-Visceral_Omentum', 'Liver', 'Heart-AtrialAppendage', 'Heart-LeftVentricle', 'Kidney-Cortex',
                            'WholeBlood', 'Muscle-Skeletal')}

else if (gwas.trait == 'Bipolar_disorder') { tissues <- c('Brain-Anteriorcingulatecortex_BA24', 'Brain-Cerebellum', 
                                                          'Brain-Cortex', 
                                          'Brain-FrontalCortex_BA9', 
                        'Brain-Spinalcord_cervicalc-1', 'WholeBlood')}

else if (gwas.trait == 'Schizophrenia') { tissues <- c('Brain-Anteriorcingulatecortex_BA24', 'Brain-Cerebellum', 
                                                       'Brain-Cortex', 
                                          'Brain-FrontalCortex_BA9', 
                        'Brain-Spinalcord_cervicalc-1', 'WholeBlood')}

else if (gwas.trait == 'Rheumatoid_arthritis') { tissues <- c('WholeBlood', 'Cells-EBV-transformedlymphocytes', 
                                                              'Adipose-Subcutaneous', 'Spleen', 'Thyroid',
                         'Liver', 'Kidney-Cortex')}

else if (gwas.trait == 'Dupuytrens_disease') { tissues <- c('Skin-NotSunExposed_Suprapubic', 'Skin-SunExposed_Lowerleg', 
                                                       'WholeBlood', 'Cells-Culturedfibroblasts')}

else if (gwas.trait == 'Chronic_obstructive_pulmonary_disease') { tissues <- c('Lung', 'WholeBlood', 'Adipose-Subcutaneous', 
                                                                          'Adipose-Visceral_Omentum',
                                          'Cells-EBV-transformedlymphocytes', 'Cells-Culturedfibroblasts')}

else if (gwas.trait == 'Hypothyroidism') { tissues <- c('Thyroid', 'Heart-AtrialAppendage', 'Heart-LeftVentricle', 
                                                        'Kidney-Cortex', 'Muscle-Skeletal',
                   'Adipose-Subcutaneous', 'Adipose-Visceral_Omentum', 'Brain-Cortex',
                            'Brain-FrontalCortex_BA9', 'Brain-Hippocampus', 'WholeBlood')}

else if (gwas.trait == 'Visceral_adipose_tissue_measurement') { tissues <- c('Adipose-Subcutaneous', 
                                                                             'Adipose-Visceral_Omentum', 'Liver', 'Pancreas',
                                        'Heart-AtrialAppendage', 'Heart-LeftVentricle', 'Muscle-Skeletal', 'Artery-Aorta',
                                         'Artery-Coronary', 'WholeBlood')}

else if (gwas.trait == 'Atopic_eczema') { tissues <- c('Skin-NotSunExposed_Suprapubic', 'Skin-SunExposed_Lowerleg', 
                                                       'WholeBlood', 
                  'Adipose-Subcutaneous', 'Adipose-Visceral_Omentum', 'Cells-Culturedfibroblasts')}
else {error = function(e) {
      message(glue("Unrecognized GWAS trait: {gwas.trait}"))
    }}
return(tissues)
}
                          
annotate.intron <- function(intron, gtf){
  tryCatch({
  strings <- c(intron)
    bed_df <- data.frame(
    chrom = str_extract(strings, "chr[0-9XY]+"),
    start = as.numeric(str_match(strings, "chr[0-9XY]+:(\\d+)")[,2]),
    end = as.numeric(str_match(strings, "chr[0-9XY]+:\\d+:(\\d+)")[,2]),
    cluster = str_match(strings, "clu_(\\d+)")[,1],
    intron = strings,
    strand = str_extract(strings, "[-+](?=:)")
  )
  intron_ranges <- makeGRangesFromDataFrame(
      bed_df, 
      keep.extra.columns = TRUE,
      ignore.strand = FALSE,
      starts.in.df.are.0based = TRUE)
  olaps <- findOverlaps(intron_ranges, gtf, type="within", select="all", ignore.strand=FALSE)
  an.intron <- intron_ranges[olaps@from] # annotated introns
  mcols(an.intron) <- cbind(mcols(gtf[olaps@to]), mcols(an.intron))
  an.intron <- as.data.table(an.intron)[
      , .(gene_name, gene_id, rk = frank(gene_name)), 
      by = .(seqnames, start, end, strand, intron, cluster)
      ][, .(gene_name, gene_id, maxrk = max(rk)), 
        by = .(seqnames, start, end, strand, intron, cluster)
      ][maxrk == 1, -c("maxrk")]
  annotated.intron <- paste0(an.intron$intron, ':', an.intron$gene_id)
  return(annotated.intron)
}, error = function(e) {
      annotated.intron <- paste0(an.intron$intron, ':NOGENE.1')
  return(annotated.intron)}
           )
           }                          

args <- commandArgs(trailingOnly=TRUE)
region <- args[1]
gwas_loci <- args[2]

# tissues <- readLines('config/tissues.txt') 
                                                    
                          
gwas.name <- getGWASName(gwas_loci)

tissues <- getRelevantTissues(gwas.name)
                          
                          

leafcutter_list <- getPermutationList(tissues, 'leafcutter')
expression_list <- getPermutationList(tissues, 'expression')

leafcutter_tissues <- search_tissues(gwas_loci, leafcutter_list, p=1e-1) %>% names()
expression_tissues <- search_tissues(gwas_loci, expression_list, p=1e-1) %>% names()
                          
tissues <- union(leafcutter_tissues, expression_tissues)                          

if (length(tissues) == 0){                          
    file.create(glue('results/coloc/hyprcoloc_results/tables/temp/{gwas_loci}.tsv'))                          
                          
    saveRDS(list(betas = NULL, ses=NULL, traits=NULL, rsid=NULL, hyprcoloc_results=NULL), 
           glue('results/coloc/hyprcoloc_results/rds/{gwas_loci}.rds'))   
    } else {
    
gtf.f <- "/project/yangili1/cfbuenabadn/leafcutter2_paper/code/annotations/gencode.v26.GRCh38.genes.csv"
gtf <- fread(gtf.f) %>% 
  .[feature == 'gene' & gene_type == 'protein_coding',
    .(seqname, start, end, gene_name, gene_id, strand)] %>%
  unique()

gtf %>% head()

gtf <- makeGRangesFromDataFrame(gtf,
    keep.extra.columns = TRUE,
    ignore.strand = FALSE
)

    
X_SummaryStats <- getColocSummaryStats(region, gwas_loci, leafcutter_list, expression_list, gtf=gtf)
                          
if (is.null(X_SummaryStats) | (length(X_SummaryStats) == 0))  {
    file.create(glue('results/coloc/hyprcoloc_results/tables/temp/{gwas_loci}.tsv'))                          
                          
    saveRDS(list(betas = NULL, ses=NULL, traits=NULL, rsid=NULL, hyprcoloc_results=NULL), 
           glue('results/coloc/hyprcoloc_results/rds/{gwas_loci}.rds'))   
}  else {                      
                          
    common_variants <- Reduce(intersect, lapply(X_SummaryStats, function(df) df$variant)) 


    X_SummaryStats <- lapply(X_SummaryStats, function(df) {
      df[df$variant %in% common_variants, ]
    })

    print('Lets see')                                            

    betas_hypr <- data.frame(row.names = X_SummaryStats[[1]]$variant)

    for (name in names(X_SummaryStats)) {
      betas_hypr[[name]] <- X_SummaryStats[[name]]$beta
    }

    ses_hypr <- data.frame(row.names = X_SummaryStats[[1]]$variant)

    for (name in names(X_SummaryStats)) {
      ses_hypr[[name]] <- X_SummaryStats[[name]]$sebeta
    }

    traits <- colnames(betas_hypr)
    rsid <- rownames(betas_hypr)
    res <- hyprcoloc(as.matrix(betas_hypr), as.matrix(ses_hypr), trait.names=traits, snp.id=rsid )

    res$results['gwas_trait'] <- gwas_loci

    write_tsv(res$results, glue('results/coloc/hyprcoloc_results/tables/temp/{gwas_loci}.tsv'))                          

    saveRDS(list(betas = betas_hypr, ses=ses_hypr, traits=traits, rsid=rsid, hyprcoloc_results=res), 
           glue('results/coloc/hyprcoloc_results/rds/{gwas_loci}.rds'))                          


}                          
}                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          