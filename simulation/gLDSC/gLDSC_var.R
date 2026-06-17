options(stringsAsFactors = F, check.names = F)

library(readr)
library(foreach)
library(doParallel)

registerDoParallel(20)

setwd("simulation")

pops = c("EUR" = 1, "AFR" = 4, "EAS" = 5)

select_block = readRDS("data/select_block.RData")

N1 = 2e5
N2_seq = c(2e4,2e5)
suffix = c("N1-200000", sprintf("N2-%d", N2_seq), sprintf("N3-%d", N2_seq))
idx = c(1,2,2,3,3)

variant = readRDS("data/variant.RData")

for(t in c("UKBB", "1kg"))
{
  for(s in 1:3)
  {
    gLDSC_var = foreach(k = 1:5, .combine = "cbind") %do%
    {
      Amatrix = read_delim(sprintf("%s/maf_g_input_53_%s/Amatrix.1.annot",
                                   "../annotations/baseline_bed_intersect", pops[idx[k]]),
                           delim = '\t')
      Amatrix = Amatrix[match(variant$rsid, Amatrix$SNP),]
      
      left_right = readRDS(sprintf(
        "gLDSC/gLDSC_results/%s/setting_%s/%s_left_right.RData",
        t, s, suffix[k]))
      
      left = foreach(block = 1:200, .combine = "+") %dopar%
      {
        left_right[[block]][[1]]
      }
      right = foreach(block = 1:200, .combine = "+") %dopar%
      {
        left_right[[block]][[2]]
      }
      tau = solve(left, right, system = "LDLt") / c(N1,N2_seq,N2_seq)[k]
      
      x = as.matrix(Amatrix[,3:56]) %*% tau[-1]
    }
    
    rownames(gLDSC_var) = variant$rsid
    colnames(gLDSC_var) = suffix
    
    saveRDS(gLDSC_var, sprintf("gLDSC/gLDSC_results/%s/setting_%d/gLDSC_var.RData", 
                               t, s))
  }
}
