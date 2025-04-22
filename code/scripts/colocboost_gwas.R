suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(BEDMatrix))
suppressMessages(library(glue))
suppressMessages(library(bedr))

suppressMessages(library(vcfR))
suppressMessages(library(genetics))
for(file in list.files("colocboost/v11/", #"LeafCutter2_ColocBoost_GTEx/final_checking", 
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

getSummaryStatsPheno <- function(gwas_loci, perm_list, perm_name, p){
    # Get a list of tissues and their corresponding dof1 values
    select_tissues <- search_tissues(gwas_loci, perm_list, p)
    
    # If no tissues match, return NULL
    if (length(select_tissues) == 0) {
        return(NULL)
    } else {
        combined_list <- do.call(c, lapply(names(select_tissues), function(tissue) {
            dof1 <- select_tissues[[tissue]]
            print(perm_name)
            getTissueRegionStats(region, gwas_loci, tissue, dof1, pheno=perm_name)
        }))
        
        return(combined_list)
    }
}


getColocSummaryStats <- function(region, gwas_loci, leafcutter_list, expression_list, p=0.1) {
    
    leafcutter_stats <- getSummaryStatsPheno(gwas_loci, leafcutter_list, 'leafcutter', p)
    expression_stats <- getSummaryStatsPheno(gwas_loci, expression_list, 'expression', p)

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
    print(gwas.file)
    ibs <- tabix(region, gwas.file)
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
    X <- tabix(region, glue('results/coloc/qtls/{tissue}/{pheno}.NominalPass/{chrom}.txt.tabix.gz'))
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
    print(dim(X))
    
    return (X)
    
}

getTissueRegionStats <- function(region, gwas_loci_id, tissue, dof, pheno) {
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
    names(X_list) <- paste(pheno_name, unique(X_stats$intron), sep='_N_N_')
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
    gt_region <- bedr::tabix(region = region, file = vcf_file) %>% 
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
                          
getGenotype <- function(region, tissue){
    chrom <- (region %>% str_split(., pattern=':'))[[1]][1]
    
    if (tissue == 'all'){
        vcf_file <- "/project/yangili1/cdai/genome_index/hs38/GTEx_v7/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files_GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz"  
    
    } else{
        vcf_file <- glue("/project/yangili1/cfbuenabadn/SpliFi/code/results/geno/GTEx/{tissue}/{chrom}.vcf.gz") 
    }
    
    gt_region <- bedr::tabix(region = region, file = vcf_file) %>% 
        filter(FORMAT == 'GT')

    rownames(gt_region) <- gt_region$ID
    
    gt_filtered <- gt_region %>%
      mutate(AF = as.numeric(str_extract(INFO, "(?<=AF=)[0-9.]+"))) %>%
      filter(AF >= 0.05)

    # Select columns starting with "GTEX-"
    gt_region <- gt_filtered[, grepl("^GTEX-", colnames(gt_filtered))]

    #gt_region <- gt_region[((gt_region == '.|.') %>% rowSums()) <= 838/20,]
    
    gt_region[gt_region == "0|0"] <- 0
    gt_region[(gt_region == "1|0") | (gt_region == "0|1")] <- 1
    gt_region[gt_region == "1|1"] <- 2
    gt_region[gt_region == ".|."] <- NA
    
    
#     return(gt_region)

    matrix_data_int <- gt_region %>% mutate_all(as.integer) #suppressWarnings(as.integer(gt_region))
    
    return(matrix_data_int)
                          
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
args <- commandArgs(trailingOnly=TRUE)
region <- args[1]
gwas.loci <- args[2]
                          
gwas.name <- getGWASName(gwas.loci)

tissues <- getRelevantTissues(gwas.name)
                          

leafcutter_list <- getPermutationList(tissues, 'leafcutter')
expression_list <- getPermutationList(tissues, 'expression')
                          
leafcutter_tissues <- search_tissues(gwas.loci, leafcutter_list, p=1e-1) %>% names()
expression_tissues <- search_tissues(gwas.loci, expression_list, p=1e-1) %>% names()
                          
tissues <- union(leafcutter_tissues, expression_tissues)

if (length(tissues) == 0){
    saveRDS(list(data = NULL), 
           glue('results/coloc/colocboost/data/{gwas_loci}.rds'))                            
                          
    saveRDS(list(res = NULL), 
           glue('results/coloc/colocboost/results/{gwas_loci}.rds'))  
} else {                         
                          
# gwas.name = str_split(gwas.loci, "_N_N_", simplify = TRUE)[, 2]

gwas.stats <- getGWASstats(region, gwas.loci, gwas.name)
                          
pheno <- c('leafcutter', 'expression')

getTissueIndividualData <- function(region, tissue, pheno, gwas.loci){
    chrom <- (region %>% str_split(':'))[[1]][1]
    file.name = glue('results/coloc/data/{tissue}/{chrom}.{pheno}.ForGWASColoc.sorted.significant_only.qqnorm.bed.gz')
    Y <- tabix(region = region, file.name = file.name)  
    Y <- Y %>% filter(gid == gwas.loci)
    rownames(Y) <- Y$pid
    Y <- Y[, grepl("^GTEX-", colnames(Y))] %>% mutate_all(as.numeric)
    
    X <- getGenotype(region, tissue)
    X <- t(X %>% na.omit())
    
    return(list(X=X, Y=Y))
}

X_list <- list()
Y_list <- list()

# Iterate through all combinations of tissues and pheno
for (t in tissues) {
  for (p in pheno) {
    tryCatch({
      result <- getTissueIndividualData(region, t, p, gwas.loci)
      X_list[[paste(t, p, sep = "_")]] <- result$X
      Y_list[[paste(t, p, sep = "_")]] <- result$Y
    }, error = function(e) {
      message(glue("Failed for tissue: {t}, pheno: {p}. Error: {e$message}"))
    })
  }
}


LD <- getLDMatrix(region)
                          
gwas.stats[[1]] <- gwas.stats[[1]] %>% filter(variant %in% colnames(LD))
                                                    
# extract all phenotype names for individual level Y
phenotype_names <- unlist(sapply(1:length(Y_list), function(iy){
  paste0(names(Y_list)[iy], "_", rownames(Y_list[[iy]]))
}))
                          
Y <- lapply(Y_list, function(y) {
  lapply(1:nrow(y), function(j) as.matrix(y[j,]))
})
Y <- Reduce("c", Y)
names(Y) <- phenotype_names

tmp <- unlist(sapply(1:length(Y_list), function(iy){
  rep(iy, nrow(Y_list[[iy]]))
}))
dict_YX <- cbind(1:length(tmp), tmp)  
    
data_reorganized <- list("X" = X_list,
                         "Y" = Y,
                         "dict_YX" = dict_YX,
                         "gwas.stats" = gwas.stats, 
                         "LD" = LD)

for(file in list.files("colocboost/v11", full.names = T)) {source(file)}
res <- colocboost(X = data_reorganized$X, # individual level genotype matrix
                  Y = data_reorganized$Y, # individual level phenotypes
                  dict_YX = data_reorganized$dict_YX, # the dictionary to link different phenotypes to genotypes
                  sumstat = data_reorganized$gwas.stats, # sumstat for GWAS
                  LD = data_reorganized$LD, # LD for GWAS
                  phenotypes = c(names(Y), names(data_reorganized$gwas.stats)), # recommendate: provide the phenotype names
                  target_idx = length(Y)+1 # target should be the index of sumstat, here is the number of individual Y + 1
                  )                                         
                                            
saveRDS(list(data = data_reorganized), glue('results/coloc/colocboost/data/{gwas.loci}.rds'))                                            

saveRDS(list(res_colocboost=res), glue('results/coloc/colocboost/results/{gwas.loci}.rds'))                                            
}