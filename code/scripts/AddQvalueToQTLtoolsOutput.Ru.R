#!/usr/bin/env Rscript

######################################################################
# @author      : bjf79 (bjf79@midway2-login1.rcc.local)
# @file        : AddQvalueToQtlToolsOutput
# @created     : Monday May 17, 2021 19:36:28 CDT
#
# @description : Add qvalue to QTLtools output
######################################################################

#Use hard coded arguments in interactive R session, else use command line args
if(interactive()){
    args <- scan(text=
                 " QTLs/QTLTools/chRNA.IR/PermutationPass.txt.gz scratch/Qvals.txt.gz   ", what='character')
} else{
    args <- commandArgs(trailingOnly=TRUE)
}

#Tissue <- args[1]
Pheno <- args[1]
FileOut <- args[2]

library(tidyverse)
library(qvalue)
library(glue)

files.in <- Pheno#glue("results/coloc/qtls/{Tissue}/{Pheno}.PermutationPass/chr{1:22}.txt.gz")

dat.in <- read_delim(files.in, delim=' ', col_types=cols()) %>% na.omit()#map_dfr(c(files.in, read_delim, delim = " ", col_types = cols()) %>% na.omit()

dat.in$q <- signif(qvalue(dat.in$adj_beta_pval)$qvalues, 5)

write_delim(dat.in, FileOut, delim=' ')
