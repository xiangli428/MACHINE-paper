library(foreach)
library(doParallel)

source('gLDSC/make_left_right_seperate3.R')

setwd("simulation")

select_block = readRDS("data/select_block.RData")

pids = c("EUR","EAS")
N1 = 2e5
N2_seq = c(2e4,2e5)
suffix = c("N1-200000", sprintf("N2-%d", N2_seq))
idx = c(1,2,2)

registerDoParallel(50)

for(s in 1:3)
{
  # dir.create(sprintf("gLDSC/gLDSC_results/UKBB/setting_%s", s), recursive = T)
  dir.create(sprintf("gLDSC/gLDSC_results/1kg/setting_%s", s), recursive = T)
  
  data = foreach(block = select_block, .combine = "rbind") %dopar%
  {
    readRDS(sprintf("data/setting_%s/%d.RData", s, block))
  }
  
  for(k in 1:3)
  {
    gwas = data.frame("SNP" = data$rsid,
                      "A1" = data$first_allele,
                      "A2" = data$alternative_alleles,
                      "N" = c(N1, N2_seq)[k],
                      "Z" = data[,12+k])
    
    # panel = sprintf("~/Documents/GWAS/project_4/LDSM/%s/chr1", pids[idx[k]])
    panel = sprintf("~/Documents/GWAS/project_4/simulation/gLDSC/LDSM_1kg/%s/chr1",
                    pids[idx[k]])
    
    result = gldsc2(panel = panel, gwas = gwas, jackknife = F, numCores = 50)
    # saveRDS(result, sprintf(
    #   "gLDSC_results/UKBB/setting_%s/%s_left_right.RData", 
    #   s, suffix[k]))
    saveRDS(result, sprintf(
      "gLDSC_results/1kg/setting_%s/%s_left_right.RData", 
      s, suffix[k]))
  }
}
