options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)

setwd("real_data/glgc/data")

dbs = c("UKBB", "META")
phenos = c("HDL", "LDL", "TG", "TC")
pids = c("EUR", "AFR", "EAS")

registerDoParallel(10)

for(pheno in phenos[3:4])
{
  min_p = foreach(db = dbs) %do%
  {
    read.delim(sprintf("%s/%s/min_p.txt", db, pheno))
  }
  names(min_p) = dbs
  
  Ls = which((rowSums(min_p[[1]][,7:9] >= 5) > 0 | rowSums(
    min_p[[2]][,7:9] >= 5) > 0) & rowSums(min_p[[1]][,4:6] > 0) > 1)
  
  for(db in dbs)
  {
    write_delim(min_p[[db]][Ls,], 
                sprintf("%s/%s/min_p_sub.txt", db, pheno), delim = '\t')
  }
}
