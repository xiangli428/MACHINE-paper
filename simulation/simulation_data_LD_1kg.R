library(readr)
library(foreach)
library(doParallel)
library(Matrix)

setwd("simulation/data/1kg")
dir.create("../LD_1kg", recursive = T)

select_block = readRDS("../select_block.RData")

pids = c("EUR", "EAS")

freq = foreach(pid = pids) %dopar%
{
  read_table(sprintf(
    "~/Documents/GWAS/data/1000GENOME/phase3/%s/chr1.frq.gz", pid))
}
names(freq) = pids

registerDoParallel(40)

foreach(block = select_block, .combine = "c") %dopar%
{
  SNP = readRDS(sprintf("../noninf/setting_%s/%d.RData", 1, block))[,1:5]
  
  SNP_1kg = foreach(pid = pids) %do%
  {
    read.delim(sprintf("%s/%s.txt.gz", pid, block))
  }
  LD_1kg = foreach(pid = pids) %do%
  {
    as.matrix(read.delim(sprintf("LD/%s/%s.ld.gz", block, pid), 
                         header = F))
  }
  names(SNP_1kg) = names(LD_1kg) = pids
  
  R_1kg = foreach(pid = pids) %do%
  {
    rownames(LD_1kg[[pid]]) = colnames(LD_1kg[[pid]]) = SNP_1kg[[pid]]$rsid
    freq_sub = freq[[pid]][is.element(freq[[pid]]$SNP, SNP_1kg[[pid]]$SNP),]
    idx = which(freq_sub$A1 != SNP_1kg[[pid]]$alternative_alleles)
    
    LD_1kg[[pid]][-idx,idx] = -LD_1kg[[pid]][-idx,idx]
    LD_1kg[[pid]][idx,-idx] = -LD_1kg[[pid]][idx,-idx]
    
    diag(LD_1kg[[pid]]) = 0
    as(LD_1kg[[pid]], "dsCMatrix")
  }
  
  saveRDS(R_1kg, sprintf("../LD_1kg/%d.RData", block))
}
