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

dir.create("GTEx_eQTL_SuSiE_V10", recursive = T)

for(i in 1:22)
{
  eQTLs = read_delim(sprintf(
    "/home/STGShare/GTEx/V10/SuSiE_chr_b37/chr%d.txt.gz", i), delim = '\t')

  data = foreach(block = filter(min_p, chromosome == i)$block,
                 .combine = "rbind") %dopar%
  {
    read_delim(sprintf("data/merge/chr%d/block_%d.txt.gz", i, block),
               delim = '\t', progress = F, show_col_types = F)
  }

  data_merge_1 = merge(data[,c(2:7,11,15)], eQTLs[,-c(3,4)], by.x = c(
    "position", "first_allele", "alternative_alleles"),
    by.y = c("variant_pos_b37", "ref_b37", "alt_b37"))

  data_merge_2 = merge(data[,c(2:7,11,15)], eQTLs[,-c(3,4)], by.x = c(
    "position", "first_allele", "alternative_alleles"),
    by.y = c("variant_pos_b37", "alt_b37", "ref_b37"))

  data_merge = rbind(data_merge_1, data_merge_2) %>% arrange(position)

  write_delim(data_merge, gzfile(sprintf("GTEx_eQTL_SuSiE_V10/chr%d.txt.gz", i)),
              delim = '\t')
}
