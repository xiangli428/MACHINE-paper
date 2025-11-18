options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)

setwd("real_data/scz2022")

pops = c("EUR" = 1, "EAS" = 5)

registerDoParallel(40)

min_p = read.delim("data/min_p_sub.txt")

foreach(pid = names(pops), .combine = "c") %do%
{
  UKBB_dir = sprintf(
    "%s/%d/imputed_genotype_unique_info-0.8_maf-0.001_hwe-1e-6",
    "~/Documents/GWAS/data/UKBB", pops[pid])
  
  foreach(l = 1:nrow(min_p), .combine = "c", .inorder = F) %dopar%
  {
    i = min_p$chromosome[l]
    block = min_p$block[l]
    
    data = read.delim(sprintf("data/%s/chr%d/block_%d.txt.gz", pid, i, block))

    dir.create(sprintf("LD/%s/chr%d/%d", pid, i, block), recursive = T)
    
    SNPs = data$rsid
    M = length(SNPs)
    LD = readMM(gzfile(sprintf("%s/LD/chr%d/%d/LD.mtx.gz", UKBB_dir, i, block)))
    variants = read.delim(sprintf("%s/LD/chr%d/%d/variant.txt", UKBB_dir, i, 
                                  block), header = F)[,1]
    id = which(is.element(variants,SNPs))
    variant2SNP = c()
    variant2SNP[id] = 0:(M-1)
    idx = which(is.element(LD@i, id-1) & is.element(LD@j, id-1))
    
    sub = new("dsTMatrix",
              "i" = variant2SNP[LD@i[idx]+1],
              "j" = variant2SNP[LD@j[idx]+1],
              "Dim" = c(M,M),
              "x" = LD@x[idx],
              "uplo" = "L")
    
    writeMM(sub, sprintf("LD/%s/chr%d/%d/LD.mtx", pid, i, block))
    system(sprintf("gzip LD/%s/chr%d/%d/LD.mtx", pid, i, block))
    
    write.table(SNPs, sprintf("LD/%s/chr%d/%d/variant.txt", pid, i, block),
                sep = '\t', row.names = F, col.names = F, quote = F)
    
    R = readMM(sprintf("LD/%s/chr%d/%d/LD.mtx.gz", pid, i, block)) %>% as.matrix()
    R_eig = eigen(R, symmetric = T)
    R_eig$values = R_eig$values + 1
    saveRDS(R_eig, sprintf("LD/%s/chr%d/%d/LD_eig.RData", pid, i, block))
    
    NULL
  }
}
