options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(ggh4x)
library(latex2exp)
library(scales)
library(Hmisc)

setwd("~/Documents/GWAS/project_4/real_data/glgc")

dbs = c("UKBB", "META")
phenos = c("HDL", "LDL", "TG", "TC")
pids = c("EUR", "AFR", "EAS")

registerDoParallel(30)


for(i in 1:22)
{
  eQTLs = read_delim(sprintf(
    "/home/STGShare/GTEx/V10/SuSiE_chr_b37/chr%d.txt.gz", i), delim = '\t')
  
  for(pheno in phenos)
  {
    min_p = read.delim(sprintf("data/%s/%s/min_p_sub.txt", "META", pheno))
    
    dir.create(sprintf("GTEx_eQTL_SuSiE_V10/%s", pheno), recursive = T)
    
    data = foreach(block = filter(min_p, chromosome == i)$block, 
                   .combine = "rbind") %dopar%
    {
      read.delim(sprintf("data/%s/%s/merge/chr%d/block_%d.txt.gz", 
                         "META", pheno, i, block))
    }
    
    data_merge_1 = merge(data[,c(2:7,11,15,19:21)], eQTLs[,-c(3,4)], by.x = c(
      "position", "first_allele", "alternative_alleles"),
      by.y = c("variant_pos_b37", "ref_b37", "alt_b37"))
    
    data_merge_2 = merge(data[,c(2:7,11,15,19:21)], eQTLs[,-c(3,4)], by.x = c(
      "position", "first_allele", "alternative_alleles"),
      by.y = c("variant_pos_b37", "alt_b37", "ref_b37"))
    
    data_merge = rbind(data_merge_1, data_merge_2) %>% arrange(position)
    
    write_delim(data_merge, gzfile(sprintf("GTEx_eQTL_SuSiE_V10/%s/chr%d.txt.gz", 
                                           pheno, i)), delim = '\t')
  }
}
