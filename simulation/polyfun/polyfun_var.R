options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(tidyr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)

setwd("simulation/polyfun/output")

N1 = 2e5
N2_seq = c(2e4,2e5)
suffix = c("N1-200000", sprintf("N2-%d", N2_seq), sprintf("N3-%d", N2_seq))

registerDoParallel(22)

SNPs = foreach(i = 1:22, .combine = "c") %dopar%
{
  read_delim(sprintf(
    "UKBB/setting_%d/%s.%d.snpvar_ridge_constrained.gz", 
    1, suffix[1], i), delim = '\t')$SNP
}

for(t in c("UKBB","1kg"))
{
  for(s in 1:3)
  {
    polyfun_var = foreach(i = 1:22, .combine = "rbind") %dopar%
    {
      foreach(k = 1:5, .combine = "cbind") %do%
        {
          read_delim(gzfile(sprintf(
            "%s/setting_%d/%s.%d.snpvar_ridge_constrained.gz", 
            t, s, suffix[k], i)), delim = '\t')$SNPVAR
        }
    }
    
    colnames(polyfun_var) = suffix
    rownames(polyfun_var) = SNPs
    saveRDS(polyfun_var, sprintf("%s/setting_%d/polyfun_var.RData", t, s))
  }
}
