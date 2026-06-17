library(readr)
library(foreach)
library(doParallel)
library(Matrix)

setwd("simulation/data")

select_block = readRDS("select_block.RData")

pids = c("EUR", "AFR", "EAS")

variant = readRDS("variant.RData")

registerDoParallel(40)

foreach(block = unique(variant$block), .combine = "c") %dopar%
{
  sub = variant[variant$block == block,]
  R_UKBB = readRDS(sprintf("LD/%s/LD_UKBB.RData", block))
  names(R_UKBB) = pids
  
  LD_1kg = foreach(pid = pids) %do%
  {
    R = as.matrix(read.delim(sprintf("LD/%s/LD_1kg_%s.ld.gz", block, pid), 
                             header = F))
    idx = which(variant[variant$block == block, sprintf("REFisA2_%s", pid)])
    R[idx,-idx] %<>% multiply_by(-1)
    R[-idx,idx] %<>% multiply_by(-1)
    
    diag(R) = 0
    
    if(sum(sub[,sprintf("MAF_%s", pid)] == 0.5) > 0)
    {
      idx = foreach(j = which(sub[,sprintf("MAF_%s", pid)] == 0.5), .combine = "c") %do%
      {
        if(cor(R[,j], R_UKBB[[pid]][,j]) < 0) j
      }
      
      if(length(idx > 0))
      {
        R[idx,-idx] %<>% multiply_by(-1)
        R[-idx,idx] %<>% multiply_by(-1)
      }
    }
    
    as(R, "dsCMatrix")
  }
  
  saveRDS(LD_1kg, sprintf("LD/%s/LD_1kg.RData", block))
  
  R_eig_list = foreach(k = 1:3) %do%
  {
    R_eig = eigen(LD_1kg[[k]], symmetric = T)
    R_eig$values = R_eig$values + 1
    R_eig
  }
  
  saveRDS(R_eig_list, sprintf("LD/%d/LD_1kg_eig.RData", block))
  
  NULL
}
