options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)

setwd("~/Documents/GWAS/project_4/real_data/glgc/data/UKBB")

phenos = c("HDL", "LDL", "TG", "TC")
pids = c("EUR", "AFR", "EAS")

num_blocks = read.table("~/Documents/GWAS/data/UKBB/num_blocks.txt")[,1]

registerDoParallel(40)

for(pheno in phenos)
{
  min_p = foreach(i = 1:22, .combine = "rbind") %do%
  {
    dir.create(sprintf("%s/merge/chr%d", pheno, i), recursive = T)
    
    data = foreach(pid = pids) %do%
    {
      foreach(block = 1:num_blocks[i], .combine = "rbind") %dopar%
      {
        if(file.exists(sprintf("%s/%s/chr%d/block_%d.txt.gz", pheno, pid, i, block)))
          read_delim(gzfile(sprintf("%s/%s/chr%d/block_%d.txt.gz", 
                                    pheno, pid, i, block)),
                     delim = '\t', show_col_types = F, progress = F,
                     col_types = list("first_allele" = "c", 
                                      "alternative_alleles" = "c"))
      }
    }
    
    data_merge = merge(
      data[[1]][,c(1:8,10:12)], data[[2]][,c(1:8,10:12)],
      by.x = 1:7, by.y = 1:7, all.x = T, all.y = T, suffixes = c(".EUR",".AFR")) %>%
      merge(data[[3]][,c(1:8,10:12)], 
            by.x = 1:7, by.y = 1:7, all.x = T, all.y = T)
    names(data_merge)[16:19] %<>% sprintf("%s.EAS", .)
    
    data_merge %<>% arrange(position)
    
    foreach(block = unique(data_merge$block), .combine = "rbind") %dopar%
    {
      sub = data_merge[data_merge$block == block,]
      
      write_delim(sub, gzfile(sprintf("%s/merge/chr%d/block_%d.txt.gz", 
                                      pheno, i, block)), delim = '\t')
      
      data.frame("chromosome" = i,
                 "block" = block,
                 "num_variants" = nrow(sub),
                 "num_variants.EUR" = sum(!is.na(sub$neglog10_pval.EUR)),
                 "num_variants.AFR" = sum(!is.na(sub$neglog10_pval.AFR)),
                 "num_variants.EAS" = sum(!is.na(sub$neglog10_pval.EAS)),
                 "nlog10minP.EUR" = max(sub$neglog10_pval.EUR, na.rm = T),
                 "nlog10minP.AFR" = max(sub$neglog10_pval.AFR, na.rm = T),
                 "nlog10minP.EAS" = max(sub$neglog10_pval.EAS, na.rm = T))
    }
  }
  
  write_delim(min_p, sprintf("%s/min_p.txt", pheno), delim = '\t')
}

