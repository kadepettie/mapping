#!/usr/bin/R

library(argparse)
library(tidyverse)
library(readxl)

# traceback to written code
options(error=function() { traceback(3); if(!interactive()) quit("no", status = 1, runLast = FALSE) })


parser <- ArgumentParser(description='Aggregate info about SNPs in peaks into single dataframe')

parser$add_argument('--sample_info',
                    type='character',
                    default=NULL,
                    help="File with individual IDs corresponding to 1000G populations.")

parser$add_argument('--pops',
                    type='character',
                    default=NULL,
                    help="Population names separated by underscores to obtain individual ids for.")

parser$add_argument('-n', '--name',
                    type='character',
                    default='fst_snps',
                    help="Full name for output files.")

parser$add_argument('-o', '--outdir',
                    type='character',
                    default='.')

opt <- parser$parse_args()

sinfo_fname <- opt$sample_info
pops <- str_split(opt$pops, '_')[[1]]
outname <- file.path(opt$outdir, opt$name)

sinfo <- read_excel(sinfo_fname, sheet=1, col_types='text')

samps <- sinfo %>%
    dplyr::filter(Population %in% pops) %>%
    pull(Sample) %>%
    write(., outname)
