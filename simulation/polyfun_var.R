options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(tidyr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)

# setwd("simulation/polyfun_output/UKBB")
setwd("simulation/polyfun_output/1kg")

N1 = 2e5
N2_seq = c(2e4,2e5)

suffix = c("N1-200000", sprintf("N2-%d", N2_seq))

registerDoParallel(22)

SNPs = foreach(i = 1:22, .combine = "c") %dopar%
{
  read_delim(gzfile(sprintf(
    "setting_%d/%s.%d.snpvar_ridge_constrained.gz", 
    1, suffix[1], i)), delim = '\t')$SNP
}

for(s in 1:3)
{
  polyfun_var = foreach(i = 1:22, .combine = "rbind") %dopar%
  {
    foreach(k = 1:5, .combine = "cbind") %do%
    {
      read_delim(gzfile(sprintf(
        "setting_%d/%s.%d.snpvar_ridge_constrained.gz", 
        s, suffix[k], i)), delim = '\t')$SNPVAR
    }
  }
  
  colnames(polyfun_var) = suffix
  rownames(polyfun_var) = SNPs
  saveRDS(polyfun_var, sprintf("setting_%d/polyfun_var.RData", s))
}
