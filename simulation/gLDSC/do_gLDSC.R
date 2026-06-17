library(foreach)
library(doParallel)

setwd("simulation")

source('../gLDSC/make_left_right_seperate3.R')

select_block = readRDS("data/select_block.RData")

pids = c("EUR","AFR","EAS")
N1 = 2e5
N2_seq = c(2e4,2e5)
suffixes = c("N1-200000", sprintf("N2-%d", N2_seq), sprintf("N3-%d", N2_seq))
idx = c(1,2,2,3,3)

registerDoParallel(50)

for(t in c("UKBB", "1kg"))
{
  for(s in 1:3)
  {
    dir.create(sprintf("gLDSC/gLDSC_results/%s/setting_%s", t, s), recursive = T)
    
    data = foreach(block = select_block, .combine = "rbind") %dopar%
    {
      readRDS(sprintf("data/setting_%s/%d.RData", s, block))
    }
    
    for(k in 1:5)
    {
      gwas = data.frame("SNP" = data$rsid,
                        "A1" = data$first_allele,
                        "A2" = data$alternative_alleles,
                        "N" = c(N1, N2_seq, N2_seq)[k],
                        "Z" = data[,20+k])
      
      if(t == "UKBB")
      {
        panel = sprintf("~/Documents/GWAS/project_4/annotations/LDSM/%s/chr1", 
                        pids[idx[k]])
      } else {
        panel = sprintf("gLDSC/LDSM_1kg/%s/chr1", pids[idx[k]])
      }
      
      result = gldsc2(panel = panel, gwas = gwas, jackknife = F, numCores = 50)
      saveRDS(result, sprintf(
        "gLDSC/gLDSC_results/%s/setting_%s/%s_left_right.RData", t, s, suffixes[k]))
    }
  }
}
