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
library(h2D2)

setwd("real_data/glgc")

dbs = c("UKBB", "GLGC")
phenos = c("HDL", "LDL", "TG", "TC")
pids = c("EUR", "AFR", "EAS")

registerDoParallel(45)

Dir = "h2D2-gLDSC"

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
      
      gLDSC_var = readRDS(sprintf("gLDSC_results/%s/%s/%s/gLDSC_var.RData", 
                                  pheno, db, pid)) %>% pmax(max(.) / 20)
      
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
                                  i, block))) %>% as("CsparseMatrix")
        R_eig = readRDS(sprintf("LD/%s/%s/chr%d/%d/LD_eig.RData", pheno, pid, 
                                i, block))
        
        M = nrow(data)
        n = N[pid, pheno]
        
        if(db == "UKBB")
        {
          z = data$beta / data$se
        } else {
          z = data$Z
        }
        
        a = gLDSC_var[data$rsid]
        a = a / mean(a) * 0.005
        
        h2D2 = Createh2D2Object(z,
                                R,
                                n,
                                data$rsid,
                                in_sample_LD = F,
                                R_eig = R_eig,
                                a = a)
        
        tic = Sys.time()
        h2D2 = h2D2_MCMC(h2D2, mcmc_n = 5500, burn_in = 500, thin = 2, 
                         get_CS = F)
        while(max(h2D2@mcmc_samples$PSRF_beta) > 1.2)
        {
          h2D2 = h2D2_MCMC(h2D2, mcmc_n = 1000, burn_in = 1000, thin = 2, 
                           get_CS = F)
        }
        h2D2@CS = h2D2_CS(h2D2)
        toc = Sys.time()
        
        saveRDS(h2D2, file = sprintf("%s/%s/%s/%s/chr%d/%d/h2D2.RData", Dir, 
                                     pheno, db, pid, i, block))
        
        data$CL = h2D2@CL
        data$beta_mean = h2D2@mcmc_mean$beta
        
        CS_info = data.frame("CS" = "",
                             "CS.size" = 0,
                             "CS.nlog10minp" = 0,
                             "CS.purity" = 0,
                             "CS.SNP" = "",
                             "CS.tags" = "")[-1,]
        
        data$CS = NA
        n_CS = length(h2D2@CS$sets)
        if(n_CS > 0)
        {
          for(s in 1:n_CS)
          {
            CS_id = sprintf("CS:%s-%d-%d-%s", pid, i, block, s)
            set = h2D2@CS$sets[[s]]
            data$CS[set] = CS_id
            tags = foreach(j = set, .combine = "union") %do%
            {
              strsplit(data$tags[j], ",")[[1]]
            }
            CS_info %<>% rbind(data.frame(
              "CS" = CS_id,
              "CS.size" = length(set),
              "CS.nlog10minp" = max(data$neglog10_pval[set]),
              "CS.purity" = h2D2@CS$purity[s,1],
              "CS.SNP" = paste(data$rsid[set], collapse = ','),
              "CS.tags" = paste(tags, collapse = ',')))
          }
        }
        
        write_delim(data, gzfile(sprintf("%s/%s/%s/%s/chr%d/%d/variant_info.txt.gz", 
                                         Dir, pheno, db, pid, i, block)), delim = '\t')
        write_delim(CS_info, sprintf("%s/%s/%s/%s/chr%d/%d/CS_info.txt", Dir, 
                                     pheno, db, pid, i, block), delim = '\t')
        
        result = data.frame(
          "pheno" = pheno,
          "db" = db,
          "pid" = pid,
          "chromosome" = i,
          "block" = block,
          "num_variants" = M,
          "max_PSRF" = max(h2D2@mcmc_samples$PSRF_beta),
          "mcmc_n" = h2D2@mcmc_samples$n_samples + h2D2@mcmc_samples$n_burnin,
          "time" = difftime(toc, tic, units = "secs"),
          "lambda" = h2D2@lambda,
          "h2_est" = sum(h2D2@a) / (sum(h2D2@a) + h2D2@b),
          "h2_mean" = h2D2@mcmc_mean$h2,
          "h2_sd" = h2D2@mcmc_sd$h2,
          "h2_beta_mean" = h2D2@mcmc_mean$h2_beta,
          "h2_beta_sd" = h2D2@mcmc_sd$h2_beta,
          "num_variants_0.95" = sum(h2D2@CL >= 0.95),
          "num_variants_0.5" = sum(h2D2@CL >= 0.5),
          "num_CS95" = length(h2D2@CS$sets),
          "num_variants_CS95" = sum(!is.na(data$CS)))
        
        write_delim(result, sprintf("%s/%s/%s/%s/chr%d/%d/result.txt", Dir, 
                                    pheno, db, pid, i, block), delim = '\t')
        
        NULL
      }
    }
  }
}

# results
for(db in dbs)
{
  for(pheno in phenos)
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

      write_delim(results, sprintf("results/%s/%s/%s/%s/results.txt", Dir, pheno,
                                   db, pid), delim = '\t')

      CS_info_all = foreach(l = which(results$num_CS95 > 0),
                            .combine = "rbind") %dopar%
      {
        i = results$chromosome[l]
        block = results$block[l]

        if(file.exists(sprintf("%s/%s/%s/%s/chr%d/%d/CS_info.txt",
                               Dir, pheno, db, pid, i, block)))
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
        data_keep = read.delim(sprintf("%s/%s/%s/%s/chr%d/%d/variant_info.txt.gz", 
                                       Dir, pheno, db, pid, i, block))
        data$CL[data$keep] = data_keep$CL
        data$beta_mean[data$keep] = data_keep$beta_mean
        data$CS[data$keep] = data_keep$CS
        
        keep_CS = data_keep$CS
        names(keep_CS) = data_keep$rsid
        data$CS[!data$keep] = foreach(j = which(!data$keep), .combine = "c") %do%
        {
          tags = strsplit(data$tags[j], ',')[[1]]
          foreach(snp = tags, .combine = "union") %do%
          {
            if(!is.na(keep_CS[snp])) keep_CS[snp]
          } %>% paste(collapse = ',')
        }
        data$CS[which(data$CS == '')] = NA
        
        filter(data, !is.na(CS))
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

        filter(data, CL >= 0.9)
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

    num_CS = merge(results[[1]][,c(4,5,18)], results[[2]][,c(4,5,18)],
                   by = c("chromosome","block"), all = T) %>%
      merge(results[[3]][,c(4,5,18)], by = c("chromosome","block"), all = T)
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
      read.delim(sprintf("results/%s/%s/%s/%s/results.txt", Dir, pheno, db, pid))
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
                 "num_variants_0.9.UKBB" = sum(data[[1]]$CL >= 0.9),
                 "num_variants_0.9.UKBB_0.1.GLGC" = sum(
                   data[[1]]$CL >= 0.9 & data[[2]]$CL < 0.1),
                 check.names = F)
    }
    
    print(sum(variant_0.9_UKBB$num_variants_0.9.UKBB))
    print(sum(variant_0.9_UKBB$num_variants_0.9.UKBB_0.1.GLGC))
    
    write_delim(variant_0.9_UKBB, sprintf("results/%s/%s/joint/%s/variant_0.9_UKBB.txt",
                                          Dir, pheno, pid), delim = '\t')
  }
}
