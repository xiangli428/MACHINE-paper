options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(tidyr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(scales)

setwd("real_data/glgc")

dbs = c("UKBB", "GLGC")
phenos = c("HDL", "LDL", "TG", "TC")
pids = c("EUR", "AFR", "EAS")

registerDoParallel(30)

Dir = "RSparsePro"

for(db in dbs)
{
  sample_size = read.delim(sprintf("data/%s/sample_size.txt", db))
  N = matrix(sample_size$n, 3, 4)
  rownames(N) = pids
  colnames(N) = phenos
  
  for(pheno in phenos)
  {
    min_p = read.delim(sprintf("data/%s/%s/min_p_sub.txt", db, pheno))[,-3]
    
    for(k in 1:3)
    {
      pid = pids[k]
      min_p_sub = min_p[min_p[9+k] > 0, c(1,2,9+k,5+k)]
      names(min_p_sub)[3] = "num_variants"
      
      foreach(l = 1:nrow(min_p_sub), .combine = "c", .inorder = F) %dopar%
      {
        i = min_p_sub$chromosome[l]
        block = min_p_sub$block[l]
        
        dir.create(sprintf("%s/%s/%s/%s/chr%d/%d", Dir, pheno, db, pid, i, block), 
                   recursive = T)
        
        data = read.delim(sprintf("data/%s/%s/%s/chr%d/block_%d.txt.gz", 
                                  db, pheno, pid, i, block))
        data = data[data$keep,-14]
        data$tags[is.na(data$tags)] = ""
        
        R = readMM(gzfile(sprintf("LD/%s/%s/chr%d/%d/LD.mtx.gz", pheno, pid, 
                                  i, block))) %>% as.matrix()
        diag(R) = 1
        
        M = nrow(data)
        n = N[pid, pheno]
        
        if(db == "UKBB")
        {
          z = data$beta / data$se
        } else {
          z = data$Z
        }
        
        write_delim(data.frame("RSID" = data$rsid, "Z" = z), sprintf(
          "%s/%s/%s/%s/chr%d/%d/z.txt", Dir, pheno, db, pid, i, block), 
          delim = '\t')
        write.table(R, sprintf("%s/%s/%s/%s/chr%d/%d/ld.txt", Dir, pheno, db, 
                               pid, i, block),
                    sep = '\t', row.names = F, col.names = F, quote = F)
        
        K = 5
        tic = Sys.time()
        system(paste("/home/r9user9/anaconda3/bin/python",
                     "~/Documents/GWAS/packages/RSparsePro_LD/src/rsparsepro_ld.py",
                     "--z", sprintf("%s/%s/%s/%s/chr%d/%d/z.txt", 
                                    Dir, pheno, db, pid, i, block),
                     "--ld", sprintf("%s/%s/%s/%s/chr%d/%d/ld.txt", Dir, pheno, 
                                     db, pid, i, block),
                     "--save", sprintf("%s/%s/%s/%s/chr%d/%d/output", Dir, pheno, 
                                       db, pid, i, block),
                     sprintf("--K %d", K)))
        output = read.delim(sprintf("%s/%s/%s/%s/chr%d/%d/output.rsparsepro.txt", 
                                    Dir, pheno, db, pid, i, block))
        
        while((max(output$cs) == K) & (K < 20))
        {
          K = K + 5
          system(paste("/home/r9user9/anaconda3/bin/python",
                       "~/Documents/GWAS/packages/RSparsePro_LD/src/rsparsepro_ld.py",
                       "--z", sprintf("%s/%s/%s/%s/chr%d/%d/z.txt", 
                                      Dir, pheno, db, pid, i, block),
                       "--ld", sprintf("%s/%s/%s/%s/chr%d/%d/ld.txt", Dir,  
                                       pheno, db, pid, i, block),
                       "--save", sprintf("%s/%s/%s/%s/chr%d/%d/output", Dir, 
                                         pheno, db, pid, i, block),
                       sprintf("--K %d", K)))
          output = read.delim(sprintf("%s/%s/%s/%s/chr%d/%d/output.rsparsepro.txt", 
                                      Dir, pheno, db, pid, i, block))
        }
        toc = Sys.time()
        
        system(sprintf("rm %s/%s/%s/%s/chr%d/%d/z.txt", Dir, pheno, db, pid, i, 
                       block))
        system(sprintf("rm %s/%s/%s/%s/chr%d/%d/ld.txt", Dir, pheno, db, pid, i, 
                       block))
        
        data$PIP = output$PIP
        data$CS = NA
        CS_info = data.frame("CS" = "",
                             "CS.size" = 0,
                             "CS.nlog10minp" = 0,
                             "CS.purity" = 0,
                             "CS.SNP" = "",
                             "CS.tags" = "")[-1,]
        n_CS = max(output$cs)
        if(n_CS > 0)
        {
          for(s in 1:n_CS)
          {
            CS_id = sprintf("CS:%s-%d-%d-%s", pid, i, block, s)
            set = which(output$cs == s)
            purity = min(abs(R[set,set]))
            
            if(purity >= 0.5)
            {
              data$CS[set] = CS_id
              tags = foreach(j = set, .combine = "union") %do%
              {
                strsplit(data$tags[j], ",")[[1]]
              }
              CS_info %<>% rbind(data.frame(
                "CS" = CS_id,
                "CS.size" = length(set),
                "CS.nlog10minp" = max(data$neglog10_pval[set]),
                "CS.purity" = purity,
                "CS.SNP" = paste(data$rsid[set], collapse = ','),
                "CS.tags" = paste(tags, collapse = ',')))
            }
          }
        }
        
        write_delim(data, gzfile(sprintf("%s/%s/%s/%s/chr%d/%d/variant_info.txt.gz", 
                                         Dir, pheno, db, pid, i, block)), 
                    delim = '\t')
        write_delim(CS_info, sprintf("%s/%s/%s/%s/chr%d/%d/CS_info.txt", Dir, 
                                     pheno, db, pid, i, block), delim = '\t')
        
        result = data.frame(
          "pheno" = pheno,
          "db" = db,
          "pid" = pid,
          "chromosome" = i,
          "block" = block,
          "num_variants" = M,
          "L" = K,
          "time" = difftime(toc, tic, units = "secs"),
          "num_variants_0.95" = sum(data$PIP >= 0.95),
          "num_variants_0.5" = sum(data$PIP >= 0.5),
          "num_CS95" = nrow(CS_info),
          "num_variants_CS95" = sum(!is.na(data$CS)))
        
        write_delim(result, sprintf("%s/%s/%s/%s/chr%d/%d/result.txt", Dir,  
                                    pheno, db, pid, i, block), delim = '\t')
        
        NULL
        
      }
    }
  }
}

# results
for(pheno in phenos)
{
  for(db in dbs)
  {
    min_p = read.delim(sprintf("data/%s/%s/min_p_sub.txt", db, pheno))[,-3]

    for(pid in pids)
    {
      dir.create(sprintf("results/%s/%s/%s/%s", Dir, pheno, db, pid),
                 recursive = T)

      results = foreach(l = 1:nrow(min_p), .combine = "rbind") %dopar%
      {
        i = min_p$chromosome[l]
        block = min_p$block[l]

        if(file.exists(sprintf("%s/%s/%s/%s/chr%d/%d/result.txt",
                               Dir, pheno, db, pid, i, block)))
          read.delim(sprintf("%s/%s/%s/%s/chr%d/%d/result.txt",
                             Dir, pheno, db, pid, i, block))
      }

      write_delim(results, sprintf("results/%s/%s/%s/%s/results.txt", Dir,
                                   pheno, db, pid), delim = '\t')

      CS_info_all = foreach(l = which(results$num_CS95 > 0),
                        .combine = "rbind") %dopar%
      {
        i = results$chromosome[l]
        block = results$block[l]

        read.delim(sprintf("%s/%s/%s/%s/chr%d/%d/CS_info.txt",
                           Dir, pheno, db, pid, i, block))
      }
      CS_info_all$CS.tags[is.na(CS_info_all$CS.tags)] = ""
      write_delim(CS_info_all, sprintf("results/%s/%s/%s/%s/CS_info.txt", Dir,
                                   pheno, db, pid), delim = '\t')

      variants_CS95 = foreach(l = which(results$num_CS95 > 0),
                              .combine = "rbind") %dopar%
      {
        i = results$chromosome[l]
        block = results$block[l]

        data = read.delim(sprintf("data/%s/%s/%s/chr%d/block_%d.txt.gz",
                                  db, pheno, pid, i, block))
        CS_info = read.delim(sprintf("%s/%s/%s/%s/chr%d/%d/CS_info.txt",
                                     Dir, pheno, db, pid, i, block))
        CS_info$CS.tags[is.na(CS_info$CS.tags)] = ""
        SNPs = foreach(s = 1:nrow(CS_info), .combine = "union") %do%
        {
          c(strsplit(CS_info$CS.SNP[s], ',')[[1]],
            strsplit(CS_info$CS.tags[s], ',')[[1]])
        }
        filter(data, rsid %in% SNPs)
      }

      write_delim(variants_CS95, gzfile(sprintf(
        "results/%s/%s/%s/%s/variants_CS95.txt.gz", Dir, pheno, db, pid)),
        delim = '\t')
      
      variants_0.9 = foreach(l = which(results$num_variants_0.5 > 0), 
                             .combine = "rbind") %dopar%
      {
        i = results$chromosome[l]
        block = results$block[l]
        
        data = read.delim(sprintf("%s/%s/%s/%s/chr%d/%d/variant_info.txt.gz", 
                                  Dir, pheno, db, pid, i, block))
        
        filter(data, PIP >= 0.9)
      }
      
      if(!is.null(variants_0.9))
      {
        write_delim(variants_0.9, gzfile(sprintf(
          "results/%s/%s/%s/%s/variants_0.9.txt.gz", Dir, pheno, db, pid)),
          delim = '\t')
      }
    }
    
    results = foreach(pid = pids) %dopar%
    {
      read.delim(sprintf("results/%s/%s/%s/%s/results.txt", Dir, pheno,
                         db, pid))
    }

    num_CS = merge(results[[1]][,c(4,5,11)], results[[2]][,c(4,5,11)],
                   by = c("chromosome","block"), all = T) %>%
      merge(results[[3]][,c(4,5,11)], by = c("chromosome","block"), all = T)
    names(num_CS)[3:5] = pids
    num_CS %<>% arrange(chromosome, block)
    num_CS[,3:5][is.na(num_CS[,3:5])] = 0

    CS_overlap = foreach(l = which(rowSums(num_CS[,3:5] > 0) > 1),
                         .combine = "rbind") %dopar%
    {
      i = num_CS$chromosome[l]
      block = num_CS$block[l]

      Ks = which(num_CS[l,3:5] > 0)
      data_pop = foreach(k = Ks) %do%
      {
        read.delim(sprintf("%s/%s/%s/%s/chr%d/%d/variant_info.txt.gz",
                           Dir, pheno, db, pids[k], i, block))
      }
      names(data_pop) = pids[Ks]

      data = read.delim(sprintf("data/%s/%s/merge/chr%d/block_%d.txt.gz",
                                db, pheno, i, block))
      for(pid in pids)
      {
        data[sprintf("CS.%s", pid)] = NA
      }
      for(k in Ks)
      {
        data[match(data_pop[[pids[k]]]$rsid,data$rsid),sprintf("CS.%s", pids[k])] =
          data_pop[[pids[k]]]$CS
      }

      if(length(Ks) == 2)
      {
        k1 = Ks[1]
        k2 = Ks[2]

        foreach(s1 = sort(unique(na.omit(data[,sprintf("CS.%s", pids[k1])]))),
                .combine = "rbind") %do%
        {
          set1 = which(data[,sprintf("CS.%s", pids[k1])] == s1)
          foreach(s2 = sort(unique(na.omit(data[set1,sprintf("CS.%s", pids[k2])]))),
                  .combine = "rbind") %do%
          {
            set2 = which(data[,sprintf("CS.%s", pids[k2])] == s2)
            un = sort(union(set1,set2))
            its = sort(intersect(set1,set2))
            tmp = data.frame("pheno" = pheno,
                             "db" = db,
                             "chromosome" = i,
                             "block" = block,
                             "CS.EUR" = "",
                             "CS.AFR" = "",
                             "CS.EAS" = "",
                             "union.size" = length(un),
                             "union.SNP" = paste(data$rsid[un], collapse = ','),
                             "intersect.size" = length(its),
                             "intersect.SNP" = paste(data$rsid[its], collapse = ','))
            tmp[,sprintf("CS.%s", pids[k1])] = s1
            tmp[,sprintf("CS.%s", pids[k2])] = s2
            tmp
          }
        }
      } else {
        rbind(
        foreach(s1 = sort(unique(na.omit(data[,sprintf("CS.%s", pids[1])]))),
                .combine = "rbind") %do%
        {
          set1 = which(data[,sprintf("CS.%s", pids[1])] == s1)

          rbind(
          foreach(s2 = sort(unique(na.omit(data[set1,sprintf("CS.%s", pids[2])]))),
                  .combine = "rbind") %do%
          {
            set2 = which(data[,sprintf("CS.%s", pids[2])] == s2)
            un = sort(union(set1, set2))
            its = sort(intersect(set1, set2))
            if(length(unique(na.omit(data[its, sprintf("CS.%s", pids[3])]))) > 0)
            {
              foreach(s3 = sort(unique(na.omit(data[its, sprintf(
                "CS.%s", pids[3])]))), .combine = "rbind") %do%
                {
                  set3 = which(data[,sprintf("CS.%s", pids[3])] == s3)
                  un = sort(union(union(set1, set2), set3))
                  its = sort(intersect(intersect(set1, set2), set3))
                  data.frame("pheno" = pheno,
                             "db" = db,
                             "chromosome" = i,
                             "block" = block,
                             "CS.EUR" = s1,
                             "CS.AFR" = s2,
                             "CS.EAS" = s3,
                             "union.size" = length(un),
                             "union.SNP" = paste(data$rsid[un], collapse = ','),
                             "intersect.size" = length(its),
                             "intersect.SNP" = paste(data$rsid[its], collapse = ','))
                }
            } else {
              data.frame("pheno" = pheno,
                         "db" = db,
                         "chromosome" = i,
                         "block" = block,
                         "CS.EUR" = s1,
                         "CS.AFR" = s2,
                         "CS.EAS" = "",
                         "union.size" = length(un),
                         "union.SNP" = paste(data$rsid[un], collapse = ','),
                         "intersect.size" = length(its),
                         "intersect.SNP" = paste(data$rsid[its], collapse = ','))
            }
          },
          foreach(s3 = sort(unique(na.omit(data[set1,sprintf("CS.%s", pids[3])]))),
                  .combine = "rbind") %do%
          {
            set3 = which(data[,sprintf("CS.%s", pids[3])] == s3)
            un = sort(union(set1, set3))
            its = intersect(set1, set3)
            if(length(unique(na.omit(data[its,sprintf("CS.%s", pids[2])]))) == 0)
            {
              data.frame("pheno" = pheno,
                         "db" = db,
                         "chromosome" = i,
                         "block" = block,
                         "CS.EUR" = s1,
                         "CS.AFR" = "",
                         "CS.EAS" = s3,
                         "union.size" = length(un),
                         "union.SNP" = paste(data$rsid[un], collapse = ','),
                         "intersect.size" = length(its),
                         "intersect.SNP" = paste(data$rsid[its], collapse = ','))
            }
          })
        },
        foreach(s2 = sort(unique(na.omit(data[,sprintf("CS.%s", pids[2])]))),
                .combine = "rbind") %do%
        {
          set2 = which(data[,sprintf("CS.%s", pids[2])] == s2)

          foreach(s3 = sort(unique(na.omit(data[set2, sprintf("CS.%s", pids[3])]))),
                  .combine = "rbind") %do%
          {
            set3 = which(data[,sprintf("CS.%s", pids[3])] == s3)
            un = sort(union(set2, set3))
            its = intersect(set2, set3)

            if(length(unique(na.omit(data[its,sprintf("CS.%s", pids[1])]))) == 0)
            {
              data.frame("pheno" = pheno,
                         "db" = db,
                         "chromosome" = i,
                         "block" = block,
                         "CS.EUR" = "",
                         "CS.AFR" = s2,
                         "CS.EAS" = s3,
                         "union.size" = length(un),
                         "union.SNP" = paste(data$rsid[un], collapse = ','),
                         "intersect.size" = length(its),
                         "intersect.SNP" = paste(data$rsid[its], collapse = ','))
            }
          }
        })
      }
    }

    if(!is.null(CS_overlap))
    {
      write_delim(CS_overlap, sprintf("results/%s/%s/%s/CS_overlap.txt", Dir,
                                      pheno, db), delim = '\t')
    }
  }
}

# joint analysis
for(pheno in phenos)
{
  for(pid in pids)
  {
    dir.create(sprintf("results/%s/%s/joint/%s", Dir, pheno, pid), recursive = T)

    results = foreach(db = dbs) %do%
    {
      read.delim(sprintf("results/%s/%s/%s/%s/results.txt", Dir, pheno, db,
                         pid))
    }

    CS_common = foreach(l = 1:nrow(results[[1]]), .combine = "rbind") %dopar%
    {
      i = results[[1]]$chromosome[l]
      block = results[[1]]$block[l]

      data = foreach(db = dbs) %do%
      {
        read.delim(sprintf("%s/%s/%s/%s/chr%d/%d/variant_info.txt.gz",
                           Dir, pheno, db, pid, i, block))
      }

      data.frame("pheno" = pheno,
                 "pid" = pid,
                 "chromosome" = i,
                 "block" = block,
                 "UKBB_nCS" = length(unique(na.omit(data[[1]]$CS))),
                 "GLGC_nCS" = length(unique(na.omit(data[[2]]$CS))),
                 "UKBB_in_GLGC" = length(unique(na.omit(
                   data[[1]]$CS[!is.na(data[[2]]$CS)]))),
                 "GLGC_in_UKBB" = length(unique(na.omit(
                   data[[2]]$CS[!is.na(data[[1]]$CS)]))))
    }

    write_delim(CS_common, sprintf("results/%s/%s/joint/%s/CS_common_num.txt",
                                   Dir, pheno, pid), delim = '\t')
    
    variant_0.9_UKBB = foreach(l = 1:nrow(results[[1]]), .combine = "rbind") %dopar%
    {
      i = results[[1]]$chromosome[l]
      block = results[[1]]$block[l]
      
      data = foreach(db = dbs) %do%
      {
        read.delim(sprintf("%s/%s/%s/%s/chr%d/%d/variant_info.txt.gz",
                           Dir, pheno, db, pid, i, block))
      }
      
      data.frame("pheno" = pheno,
                 "pid" = pid,
                 "chromosome" = i,
                 "block" = block,
                 "num_variants_0.9.UKBB" = sum(data[[1]]$PIP >= 0.9),
                 "num_variants_0.9.UKBB_0.1.GLGC" = sum(
                   data[[1]]$PIP >= 0.9 & data[[2]]$PIP < 0.1),
                 check.names = F)
    }
    
    print(sum(variant_0.9_UKBB$num_variants_0.9.UKBB))
    print(sum(variant_0.9_UKBB$num_variants_0.9.UKBB_0.1.GLGC))
    
    write_delim(variant_0.9_UKBB, sprintf("results/%s/%s/joint/%s/variant_0.9_UKBB.txt",
                                          Dir, pheno, pid), delim = '\t')
  }
}
