options(stringsAsFactors = F, check.names = F)

library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)

setwd("~/Documents/GWAS/project_4/simulation_new/data")

registerDoParallel(50)

select_block = readRDS("select_block.RData")
variant = readRDS("variant.RData")

pids = c("EUR", "AFR", "EAS")

foreach(id = 1:200, .combine = "c", .inorder = F) %dopar%
{
  R = readRDS(sprintf("LD/%d/LD_1kg.RData", select_block[id]))
  
  for(k in 1:3)
  {
    dir.create(sprintf("LD_1kg_batch/%s/chr1/%s", pids[k], id),
               recursive = T)
    
    writeMM(R[[k]], sprintf("LD_1kg_batch/%s/chr1/%s/LD.mtx", pids[k], id))
    system(sprintf("gzip LD_1kg_batch/%s/chr1/%s/LD.mtx", pids[k], id))
    
    write.table(variant$rsid[variant$block == select_block[id]], sprintf(
      "LD_1kg_batch/%s/chr1/%s/variant.txt", pids[k], id),
      sep = '\t', row.names = F, col.names = F, quote = F)
  }
}