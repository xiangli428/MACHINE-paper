options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)

source('gLDSC/make_left_right_seperate3.R')

setwd("real_data/glgc")

dbs = c("UKBB", "GLGC")
phenos = c("HDL", "LDL", "TG", "TC")
pids = c("EUR", "AFR", "EAS")

num_blocks = read.table("~/Documents/GWAS/data/UKBB/num_blocks.txt")[,1]

nc = 90
registerDoParallel(nc)

for(db in dbs)
{
  sample_size = read.delim(sprintf("data/%s/sample_size.txt", db))
  N = matrix(sample_size$n, 3, 4)
  rownames(N) = pids
  colnames(N) = phenos
  
  for(pheno in phenos)
  {
    for(pid in pids)
    {
      dir.create(sprintf("gLDSC_results/%s/%s/%s", pheno, db, pid), 
                 recursive = T, showWarnings = F)
      
      for(i in 1:22)
      {
        if(!file.exists(sprintf("gLDSC_results/%s/%s/%s/chr%s.RData", 
                                pheno, db, pid, i)))
        {
          print(sprintf("%s %s %s chr%s start", pheno, db, pid, i))
          
          panel = sprintf("~/Documents/GWAS/project_4/annotations/LDSM/%s/chr%s",
                          pid, i)
          
          data_all = foreach(block = 1:num_blocks[i], .combine = "rbind") %dopar%
          {
            if(file.exists(sprintf("data/%s/%s/%s/chr%d/block_%d.txt.gz", 
                                   db, pheno, pid, i, block)))
            {
              data = read_delim(sprintf("data/%s/%s/%s/chr%d/block_%d.txt.gz", 
                                        db, pheno, pid, i, block),
                                delim = '\t', show_col_types = F, progress = F,
                                col_types = list("first_allele" = "c", 
                                                 "alternative_alleles" = "c"))
              data = data[data$keep,-14]
              data$tags[is.na(data$tags)] = ""
              data
            }
          }
          
          n = N[pid, pheno]
          
          if(db == "UKBB")
          {
            z = data_all$beta / data_all$se
            z = z * sqrt(n / (n + z^2))
          } else {
            z = data_all$Z
          }
          
          gwas = data.frame("SNP" = data_all$rsid,
                            "A1" = data_all$first_allele,
                            "A2" = data_all$alternative_alleles,
                            "N" = n,
                            "Z" = z)
          
          
          result = gldsc2(panel = panel, gwas = gwas, jackknife = F, numCores = nc)
          saveRDS(result, sprintf("gLDSC_results/%s/%s/%s/chr%s.RData", 
                                  pheno, db, pid, i))
        }
      }
    }
  }
}
