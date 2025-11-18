options(stringsAsFactors = F, check.names = F)

library(readr)
library(tidyr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)

setwd("real_data/scz2022")

min_p = read.delim("data/min_p_sub.txt")

registerDoParallel(25)

for(anno in c("NSF","NSM","NSN","R3","R5","SYN","U3","U5"))
{
  dir.create(sprintf("dbSNP151/%s", anno), recursive = T)
}

for(i in 1:22)
{
  data = foreach(block = filter(min_p, chromosome == i)$block, 
                 .combine = "rbind") %dopar%
  {
    read_delim(sprintf("data/merge/chr%d/block_%d.txt.gz", i, block), 
               delim = '\t', progress = F, show_col_types = F)
  }
  
  for(anno in c("NSF","NSM","NSN","R3","R5","SYN","U3","U5"))
  {
    anno_SNPs = read_delim(sprintf(
      "%s/human_9606_b151_GRCh37p13/coding_txt/%s/00-All_%s_chr%d.txt.gz",
      "~/Documents/GWAS/data/dbSNP", anno, anno, i))
    
    data_merge_1 = merge(data[,c(2:7,11,15)], anno_SNPs, by.x = c(
      "position", "first_allele", "alternative_alleles"),
      by.y = c("POS", "REF", "ALT"))
    
    data_merge_2 = merge(data[,c(2:7,11,15)], anno_SNPs, by.x = c(
      "position", "first_allele", "alternative_alleles"),
      by.y = c("POS", "ALT", "REF"))
    
    data_merge = rbind(data_merge_1, data_merge_2) %>% arrange(position)
    
    write_delim(data_merge, gzfile(sprintf(
      "dbSNP151/%s/chr%d.txt.gz", anno, i)), delim = '\t')
  }
}