library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

f_in <- args[1]
stats <- args[2]
f_out <- args[3]

chrom_list <- c(paste0('chr', 1:22), 'chrX')

# if (stats == 'beta_se') {
dat <- read_tsv(f_in, 
            col_names = c('chrom', 'start', 'end', 'var_id', 'var_alt_id', 'P', 'beta', 'SE')) %>%
            filter(chrom %in% chrom_list)
beta <- dat$beta
SE <- dat$SE   

dat_out <- dat %>% 
    select(chrom, start, end, P, var_id, var_alt_id) %>%
    mutate(beta = beta) %>%
    mutate(SE=SE) %>%
    write_tsv(f_out)
