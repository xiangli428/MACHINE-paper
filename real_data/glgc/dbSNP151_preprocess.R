options(stringsAsFactors = F, check.names = F)

library(readr)
library(tidyr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)

setwd("~/Documents/GWAS/project_4/real_data/glgc")

dbs = c("UKBB", "META")
phenos = c("HDL", "LDL", "TG", "TC")
pids = c("EUR", "AFR", "EAS")

registerDoParallel(30)

for(anno in c("NSF","NSM","NSN","R3","R5","SYN","U3","U5"))
{
  for(i in 1:22)
  {
    anno_SNPs = read_delim(sprintf(
      "%s/human_9606_b151_GRCh37p13/coding_txt/%s/00-All_%s_chr%d.txt.gz",
      "~/Documents/GWAS/data/dbSNP", anno, anno, i))
    
    for(pheno in phenos)
    {
      min_p = read.delim(sprintf("data/%s/%s/min_p_sub.txt", "META", pheno))
      
      dir.create(sprintf("dbSNP151/%s/%s", pheno, anno), recursive = T)
      
      data = foreach(block = filter(min_p, chromosome == i)$block, 
                     .combine = "rbind") %dopar%
      {
        read.delim(sprintf("data/%s/%s/merge/chr%d/block_%d.txt.gz", 
                           "META", pheno, i, block))
      }
      
      data_merge_1 = merge(data[,c(2:7,11,15,19:21)], anno_SNPs, by.x = c(
        "position", "first_allele", "alternative_alleles"),
        by.y = c("POS", "REF", "ALT"))
      
      data_merge_2 = merge(data[,c(2:7,11,15,19:21)], anno_SNPs, by.x = c(
        "position", "first_allele", "alternative_alleles"),
        by.y = c("POS", "ALT", "REF"))
      
      data_merge = rbind(data_merge_1, data_merge_2) %>% arrange(position)
      
      write_delim(data_merge, gzfile(sprintf(
        "dbSNP151/%s/%s/chr%d.txt.gz", pheno, anno, i)), delim = '\t')
    }
  }
}