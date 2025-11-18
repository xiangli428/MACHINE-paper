options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)

setwd("real_data/glgc")

phenos = c("HDL", "LDL", "TG", "TC")
pops = c("EUR" = 1, "AFR" = 4, "EAS" = 5)

registerDoParallel(20)

for(pheno in phenos)
{
  min_p_sub = read.delim(sprintf("data/UKBB/%s/min_p_sub.txt", pheno))
  
  for(pid in names(pops))
  {
    UKBB_dir = sprintf(
      "%s/%d/imputed_genotype_unique_info-0.8_maf-0.001_hwe-1e-6",
      "~/Documents/GWAS/data/UKBB", pops[pid])
    
    foreach(l = which(min_p_sub[,sprintf("num_variants.%s", pid)] > 0),
            .combine = "c") %dopar%
    {
      i = min_p_sub$chromosome[l]
      block = min_p_sub$block[l]
      
      dir.create(sprintf("LD/%s/%s/chr%s/%s", pheno, pid, i, block), 
                 recursive = T)
      
      data = read_delim(gzfile(sprintf(
        "data/UKBB/%s/%s/chr%d/block_%d.txt.gz", pheno, pid, i, block)),
        delim = '\t', show_col_types = F, progress = F,
        col_types = list("first_allele" = "c", "alternative_alleles" = "c"))
      
      R = readMM(gzfile(sprintf("%s/LD/chr%d/%d/LD.mtx.gz", UKBB_dir, i, block)))
      rownames(R) = colnames(R) = read.delim(sprintf(
        "%s/LD/chr%d/%d/variant.txt", UKBB_dir, i, block), header = F)[,1]
      R = R[data$rsid, data$rsid]
      
      writeMM(R, sprintf("LD/%s/%s/chr%d/%d/LD.mtx", pheno, pid, i, block))
      system(sprintf("gzip LD/%s/%s/chr%d/%d/LD.mtx", pheno, pid, i, block))
      
      R_eig = eigen(R, symmetric = T)
      R_eig$values = R_eig$values + 1
      saveRDS(R_eig, sprintf("LD/%s/%s/chr%d/%d/LD_eig.RData", 
                             pheno, pid, i, block))
      
      NULL
    }
  }
}
