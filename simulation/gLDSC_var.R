options(stringsAsFactors = F, check.names = F)

library(readr)
library(foreach)
library(doParallel)

registerDoParallel(20)

setwd("simulation")

select_block = readRDS("data/select_block.RData")

N1 = 2e5
N2_seq = c(2e4,2e5)
suffix = c("N1-200000", sprintf("N2-%d", N2_seq))

Amatrix = read_delim("../annotations/baseline_bed_intersect/g_input_52_1/Amatrix.1.annot",
                     delim = '\t')
SNPs = foreach(block = select_block, .combine = "c") %dopar%
{
  readRDS(sprintf("data/setting_%s/%d.RData", 1, block))$rsid
}
Amatrix = Amatrix[match(SNPs, Amatrix$SNP),]

for(s in 1:3)
{
  gLDSC_var = foreach(k = 1:5, .combine = "cbind") %do%
  {
    # left_right = readRDS(sprintf(
    #   "gLDSC_results/UKBB/setting_%s/%s_left_right.RData",
    #   s, suffix[k]))
    left_right = readRDS(sprintf(
      "gLDSC_results/1kg/setting_%s/%s_left_right.RData",
      s, suffix[k]))
    
    left = foreach(block = 1:200, .combine = "+") %dopar%
    {
      left_right[[block]][[1]]
    }
    right = foreach(block = 1:200, .combine = "+") %dopar%
    {
      left_right[[block]][[2]]
    }
    tau = solve(left, right, system = "LDLt") / c(N1,N2_seq)[k]
    
    as.matrix(Amatrix[,3:55]) %*% tau[-1]
  }
  
  rownames(gLDSC_var) = SNPs
  colnames(gLDSC_var) = suffix
  
  # saveRDS(gLDSC_var, sprintf("gLDSC_results/UKBB/setting_%d/gLDSC_var.RData", s))
  saveRDS(gLDSC_var, sprintf("gLDSC_results/1kg/setting_%d/gLDSC_var.RData", s))
}
