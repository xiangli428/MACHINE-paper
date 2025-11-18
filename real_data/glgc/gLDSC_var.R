options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)

setwd("real_data/glgc")

dbs = c("UKBB", "GLGC")
phenos = c("HDL", "LDL", "TG", "TC")
pids = c("EUR", "AFR", "EAS")
pops = c("EUR" = 1, "AFR" = 4, "EAS" = 5)

num_blocks = read.table("~/Documents/GWAS/data/UKBB/num_blocks.txt")[,1]

registerDoParallel(60)

for(db in dbs[2])
{
  sample_size = read.delim(sprintf("data/%s/sample_size.txt", db))
  N = matrix(sample_size$n, 3, 4)
  rownames(N) = pids
  colnames(N) = phenos
  
  for(pheno in phenos)
  {
    for(pid in pids[2])
    {
      left_right = foreach(i = 1:22) %dopar%
      {
        readRDS(sprintf("gLDSC_results/%s/%s/%s/chr%s.RData", 
                        pheno, db, pid, i))
      }
      
      left = foreach(i = 1:22, .combine = "+") %dopar%
      {
        foreach(block = 1:length(left_right[[i]]), .combine = "+") %do%
        {
          left_right[[i]][[block]][[1]]
        }
      }
      
      right = foreach(i = 1:22, .combine = "+") %dopar%
      {
        foreach(block = 1:length(left_right[[i]]), .combine = "+") %do%
        {
          left_right[[i]][[block]][[2]]
        }
      }
      
      Taus = solve(left, right, system = "LDLt")
      gLDSC_res = list("Taus" = Taus[-1] / N[pid, pheno],
                       "intercept" = 1+Taus[1])
      saveRDS(gLDSC_res, sprintf("gLDSC_results/%s/%s/%s/gLDSC_res.RData", 
                                 pheno, db, pid))
      
      gLDSC_var = foreach(i = 1:22, .combine = "c") %do%
      {
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
        
        Amatrix = read_delim(sprintf(
          "%s/baseline_bed_intersect/maf_g_input_53_%s/Amatrix.%s.annot", 
          "../../annotations", pops[pid], i), delim = '\t', show_col_types = F)
        Amatrix = Amatrix[match(data_all$rsid, Amatrix$SNP),]
        
        v = (as.matrix(Amatrix[,3:56]) %*% gLDSC_res$Taus)[,1]
        names(v) = data_all$rsid
        v
      }
      
      saveRDS(gLDSC_var, sprintf("gLDSC_results/%s/%s/%s/gLDSC_var.RData", 
                                 pheno, db, pid))
    }
  }
}
