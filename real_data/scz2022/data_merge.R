options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)

setwd("real_data/scz2022/data")

pops = c("EUR" = 1, "EAS" = 5)
num_blocks = read.table("~/Documents/GWAS/data/UKBB/num_blocks.txt")[,1]

registerDoParallel(40)

min_p = foreach(i = 1:22, .combine = "rbind") %do%
{
  dir.create(sprintf("merge/chr%s", i), recursive = T)
  
  data = foreach(pid = names(pops)) %do%
  {
    foreach(block = 1:num_blocks[i], .combine = "rbind") %dopar%
    {
      if(file.exists(sprintf("%s/chr%d/block_%s.txt.gz", pid, i, block)))
      {
        read_delim(sprintf("%s/chr%d/block_%s.txt.gz", pid, i, block),
                   delim = '\t', show_col_types = F, progress = F,
                   col_types = list("first_allele" = "c", 
                                    "alternative_alleles" = "c"))
      }
    }
  }
  
  data_merge = merge(data[[1]][,c(1:8,13:15)], data[[2]][,c(1:8,13:15)],
                     by = 1:7, all.x = T, all.y = T, suffixes = c(".EUR",".EAS"))
  
  data_merge %<>% arrange(position)
  
  foreach(block = unique(data_merge$block), .combine = "rbind") %dopar%
  {
    sub = data_merge[data_merge$block == block,]
    
    write_delim(sub, gzfile(sprintf("merge/chr%d/block_%d.txt.gz", 
                                    i, block)), delim = '\t')
    
    data.frame("chromosome" = i,
               "block" = block,
               "num_variants" = nrow(sub),
               "num_variants.EUR" = sum(!is.na(sub$PVAL.EUR)),
               "num_variants.EAS" = sum(!is.na(sub$PVAL.EAS)),
               "minP.EUR" = min(sub$PVAL.EUR, na.rm = T),
               "minP.EAS" = min(sub$PVAL.EAS, na.rm = T))
  }
}

write_delim(min_p, "min_p.txt", delim = '\t')

min_p_sub = filter(min_p, minP.EUR < 1e-5 | minP.EAS < 1e-5)
min_p_sub = min_p_sub[-c(278,363),]
write_delim(min_p_sub, "min_p_sub.txt", delim = '\t')
