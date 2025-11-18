options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(tidyr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)

setwd("~/Documents/GWAS/project_4/real_data/scz2022")

registerDoParallel(22)

pops = c("EUR" = 1, "EAS" = 5)
N = c("EUR" = 130644, "EAS" = 30761)

num_blocks = read.table("~/Documents/GWAS/data/UKBB/num_blocks.txt")[,1]

for(pid in names(pops))
{
  left_right = foreach(i = 1:22) %dopar%
  {
    res = readRDS(sprintf("gLDSC_results/%s/chr/chr%s.Rdata", pid, i))
    left = foreach(block = 1:length(res), .combine = "+") %do%
    {
      res[[block]][[1]]
    }
    right = foreach(block = 1:length(res), .combine = "+") %do%
    {
      res[[block]][[2]]
    }
    list(left, right)
  }
  
  left = foreach(i = 1:22, .combine = "+") %dopar%
  {
    left_right[[i]][[1]]
  }
  
  right = foreach(i = 1:22, .combine = "+") %dopar%
  {
    left_right[[i]][[2]]
  }
  
  Taus = solve(left, right, system = "LDLt")
  gLDSC_res = list("Taus" = Taus[-1] / N[pid],
                   "intercept" = 1+Taus[1])
  saveRDS(gLDSC_res, sprintf("gLDSC_results/%s/gLDSC_res.RData", pid))
  
  gLDSC_var = foreach(i = 1:22, .combine = "c") %do%
  {
    data_all = foreach(block = 1:num_blocks[i], .combine = "rbind") %dopar%
    {
      if(file.exists(sprintf("data/%s/chr%d/block_%d.txt.gz", pid, i, block)))
      {
        read_delim(sprintf("data/%s/chr%d/block_%d.txt.gz", pid, i, block),
                   delim = '\t', show_col_types = F, progress = F,
                   col_types = list("first_allele" = "c", 
                                    "alternative_alleles" = "c"))
      }
    }
    
    Amatrix = read_delim(sprintf(
      "%s/baseline_bed_intersect/maf_g_input_53_%s/Amatrix.%s.annot", 
      "../../annotations", pops[pid], i), delim = '\t')
    Amatrix = Amatrix[match(data_all$rsid, Amatrix$SNP),]
    
    v = (as.matrix(Amatrix[,3:56]) %*% gLDSC_res$Taus)[,1]
    names(v) = data_all$rsid
    v
  }
  
  saveRDS(gLDSC_var, sprintf("gLDSC_results/%s/gLDSC_var.RData", pid))
}
