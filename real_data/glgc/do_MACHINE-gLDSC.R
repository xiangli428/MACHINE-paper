options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(scales)
library(MACHINE)

setwd("real_data/glgc")

dbs = c("UKBB", "GLGC")
phenos = c("HDL", "LDL", "TG", "TC")
pids = c("EUR", "AFR", "EAS")

registerDoParallel(34)

Dir = "MACHINE-gLDSC"

for(db in dbs)
{
  sample_size = read.delim(sprintf("data/%s/sample_size.txt", db))
  N = matrix(sample_size$n, 3, 4)
  rownames(N) = pids
  colnames(N) = phenos
  
  for(pheno in phenos)
  {
    min_p = read.delim(sprintf("data/%s/%s/min_p_sub.txt", db, pheno))[
      ,c(1,2,10:13,7:9)]
    
    gLDSC_var = foreach(pid = pids) %dopar%
    {
      readRDS(sprintf("gLDSC_results/%s/%s/%s/gLDSC_var.RData", 
                      pheno, db, pid)) %>% pmax(max(.) / 20)
    }
    names(gLDSC_var) = pids
    
    foreach(l = 1:nrow(min_p), .combine = "c", .inorder = F) %dopar%
    {
      i = min_p$chromosome[l]
      block = min_p$block[l]
      
      dir.create(sprintf("%s/%s/%s/chr%d/%d", Dir, pheno, db, i, block),
                 recursive = T)
      
      data = read.delim(sprintf("data/%s/%s/merge/chr%d/block_%d.txt.gz", 
                                db, pheno, i, block))
      M = nrow(data)
      
      used_pids = foreach(pid = pids, .combine = "c") %do%
      {
        if(!all(is.na(data[,sprintf("neglog10_pval.%s", pid)]))) pid
      }
      
      R = foreach(pid = used_pids) %do%
      {
        readMM(sprintf("LD/%s/%s/chr%d/%d/LD.mtx.gz", pheno, pid, i, block))
      }
      R_eig = foreach(pid = used_pids) %do%
      {
        readRDS(sprintf("LD/%s/%s/chr%d/%d/LD_eig.RData", pheno, pid, i, block))
      }
      
      if(db == "UKBB")
      {
        Z = foreach(pid = used_pids, .combine = "cbind") %do%
        {
          z = data[,sprintf("beta.%s", pid)] / data[,sprintf("se.%s", pid)]
        }
      } else {
        Z = data[,sprintf("Z.%s", used_pids)]
      }
      
      gLDSC_var_sub = foreach(pid = used_pids, .combine = "cbind") %do%
      {
        gLDSC_var[[pid]][data$rsid] / mean(gLDSC_var[[pid]][data$rsid], na.rm = T)
      }
      a = rowMeans(gLDSC_var_sub, na.rm = T)
      a = a / mean(a) * 0.005
      c = gLDSC_var_sub / rowSums(gLDSC_var_sub, na.rm = T) * 
        rowSums(!is.na(gLDSC_var_sub), na.rm = T) * 0.2
      
      MACHINE = CreateMACHINEObject(Z, 
                                    R, 
                                    N[used_pids, pheno],
                                    SNP_ID = data$rsid,
                                    POP_ID = used_pids,
                                    in_sample_LD = F,
                                    R_eig = R_eig,
                                    a = a,
                                    c = c)
      rm(R, R_eig)
      tic = Sys.time()
      MACHINE = MACHINE_MCMC(MACHINE, mcmc_n = 5500, burn_in = 500, thin = 1,
                             get_CS = F)
      while(max(MACHINE@mcmc_samples$PSRF_beta) > 1.2)
      {
        MACHINE = MACHINE_MCMC(MACHINE, mcmc_n = 1000, burn_in = 1000, thin = 1,
                               get_CS = F)
      }
      MACHINE@CS = MACHINE_CS(MACHINE)
      toc = Sys.time()
      
      saveRDS(MACHINE, file = sprintf("%s/%s/%s/chr%d/%d/MACHINE.RData", Dir, 
                                      pheno, db, i, block))
        
      for(pid in pids)
      {
        data[sprintf("sigma2_mean.%s", pid)] = NA
        data[sprintf("CL.%s", pid)] = NA
        data[sprintf("beta_mean.%s", pid)] = NA
        data[sprintf("CS.%s", pid)] = NA
      }
      
      n_CS = foreach(pid = pids, .combine = "c") %do%
      {
        length(MACHINE@CS[[pid]]$sets)
      }
      names(n_CS) = pids
      
      CS_info = data.frame("CS" = "",
                           "CS.size" = 0,
                           "CS.neglog10_pval" = 0,
                           "CS.purity" = 0,
                           "CS.SNP" = "",
                           "CS.tags" = "")[-1,]
      
      for(pid in used_pids)
      {
        idx = which(!is.na(data[[sprintf("neglog10_pval.%s", pid)]]))
        data[idx, sprintf("sigma2_mean.%s", pid)] = 
          MACHINE@mcmc_mean$sigma2[idx, pid]
        data[idx, sprintf("CL.%s", pid)] = 
          MACHINE@CL[idx, pid]
        data[idx, sprintf("beta_mean.%s", pid)] = 
          MACHINE@mcmc_mean$beta[idx, pid]
        
        if(n_CS[pid] > 0)
        {
          for(s in 1:n_CS[pid])
          {
            CS_id = sprintf("CS:%s-%d-%d-%s", pid, i, block, s)
            set = MACHINE@CS[[pid]]$sets[[s]]
            tags = foreach(j = set, .combine = "union") %do%
              {
                strsplit(data$tags[j], ",")[[1]]
              }
            data[set, sprintf("CS.%s", pid)] = CS_id
            CS_info %<>% rbind(data.frame(
              "CS" = CS_id,
              "CS.size" = length(set),
              "CS.nlog10minp" = max(data[set, sprintf("neglog10_pval.%s", pid)]),
              "CS.purity" = MACHINE@CS[[pid]]$purity[s,1],
              "CS.SNP" = paste(data$rsid[set], collapse = ','),
              "CS.tags" = paste(tags, collapse = ',')))
          }
        }
      }
      
      write_delim(data, gzfile(sprintf("%s/%s/%s/chr%d/%d/variant_info.txt.gz", 
                                       Dir, pheno, db, i, block)), 
                  delim = '\t')
      write_delim(CS_info, sprintf("%s/%s/%s/chr%d/%d/CS_info.txt", Dir, pheno, 
                                   db, i, block), delim = '\t')
      
      result = data.frame(
        "pheno" = pheno,
        "db" = db,
        "chromosome" = i,
        "block" = block,
        "num_variants" = MACHINE@M,
        "max_PSRF" = max(MACHINE@mcmc_samples$PSRF_beta),
        "max_beta" = lapply(MACHINE@mcmc_samples$beta, abs) %>% sapply(max) %>% max,
        "max_betamean" = max(abs(MACHINE@mcmc_mean$beta)),
        "mcmc_n" = MACHINE@mcmc_samples$n_samples + MACHINE@mcmc_samples$n_burnin,
        "time" = difftime(toc, tic, units = "secs"))
      
      for(pid in pids)
      {
        if(pid %in% used_pids)
        {
          result[[sprintf("h2_beta_mean.%s", pid)]] =
            MACHINE@mcmc_mean$h2_beta[pid]
          result[[sprintf("h2_beta_sd.%s", pid)]] =
            MACHINE@mcmc_sd$h2_beta[pid]
          result[[sprintf("num_variants_0.95.%s", pid)]] = 
            sum(MACHINE@CL[,pid] >= 0.95)
          result[[sprintf("num_variants_0.5.%s", pid)]] = 
            sum(MACHINE@CL[,pid] >= 0.5)
          result[[sprintf("num_CS95.%s", pid)]] = 
            length(MACHINE@CS[[pid]]$sets)
          result[[sprintf("num_variants_CS95.%s", pid)]] = 
            if(length(MACHINE@CS[[pid]]$sets) > 0) sum(
              sapply(MACHINE@CS[[pid]]$sets, length)) else 0
        } else {
          result[[sprintf("h2_beta_mean.%s", pid)]] = 0
          result[[sprintf("h2_beta_sd.%s", pid)]] = 0
          result[[sprintf("num_variants_0.95.%s", pid)]] = 0
          result[[sprintf("num_variants_0.5.%s", pid)]] = 0
          result[[sprintf("num_CS95.%s", pid)]] = 0
          result[[sprintf("num_variants_CS95.%s", pid)]] = 0
        }
      }
      
      write_delim(result, sprintf("%s/%s/%s/chr%d/%d/result.txt", Dir, pheno, 
                                  db, i, block), delim = '\t')
      
      NULL
    }
  }
}

# results
for(db in dbs)
{
  for(pheno in phenos)
  {
    min_p = read.delim(sprintf("data/%s/%s/min_p_sub.txt", db, pheno))

    dir.create(sprintf("results/%s/%s/%s", Dir, pheno, db), recursive = T)

    results = foreach(l = 1:nrow(min_p), .combine = "rbind") %dopar%
    {
      i = min_p$chromosome[l]
      block = min_p$block[l]

      if(file.exists(sprintf("%s/%s/%s/chr%d/%d/result.txt",
                             Dir, pheno, db, i, block)))
        read.delim(sprintf("%s/%s/%s/chr%d/%d/result.txt",
                           Dir, pheno, db, i, block))
    }

    write_delim(results, sprintf("results/%s/%s/%s/results.txt", Dir, pheno,
                                 db), delim = '\t')

    summary = foreach(pid = pids, .combine = "rbind") %do%
    {
      data.frame(
        "pid" = pid,
        "num_variants_0.95" = sum(results[,sprintf("num_variants_0.95.%s", pid)]),
        "num_variants_0.5" = sum(results[,sprintf("num_variants_0.5.%s", pid)]),
        "num_CS95" = sum(results[,sprintf("num_CS95.%s", pid)]),
        "num_variants_CS95" = sum(results[,sprintf("num_variants_CS95.%s", pid)]))
    }

    write_delim(summary, sprintf("results/%s/%s/%s/summary.txt", Dir, pheno, db),
                delim = '\t')

    CS_info_all = foreach(l = 1:nrow(results), .combine = "rbind") %dopar%
    {
      i = results$chromosome[l]
      block = results$block[l]

      if(file.exists(sprintf("%s/%s/%s/chr%d/%d/CS_info.txt",
                             Dir, pheno, db, i, block)))
        read.delim(sprintf("%s/%s/%s/chr%d/%d/CS_info.txt",
                           Dir, pheno, db, i, block))
    }
    CS_info_all$CS.tags[is.na(CS_info_all$CS.tags)] = ""

    write_delim(CS_info_all, sprintf("results/%s/%s/%s/CS_info.txt", Dir, pheno,
                                     db), delim = '\t')
    
    variants_CS95 = foreach(l = which(rowSums(results[,c(15,21,27)]) > 0),
                            .combine = "rbind") %dopar%
    {
      i = results$chromosome[l]
      block = results$block[l]
      
      data = read.delim(sprintf("data/%s/%s/merge/chr%d/block_%d.txt.gz",
                                db, pheno, i, block))
      data_keep = read.delim(sprintf("%s/%s/%s/chr%d/%d/variant_info.txt.gz", 
                                     Dir, pheno, db, i, block))
      for(pid in pids)
      {
        data[sprintf("sigma2_mean.%s", pid)] = NA
        data[sprintf("CL.%s", pid)] = NA
        data[sprintf("beta_mean.%s", pid)] = NA
        data[sprintf("CS.%s", pid)] = NA
      }
      data[data$keep, 22:33] = data_keep[,21:32]
      
      for(pid in pids)
      {
        if(any(!is.na(data_keep[,sprintf("CS.%s", pid)])))
        {
          keep_CS = data_keep[,sprintf("CS.%s", pid)]
          names(keep_CS) = data_keep$rsid
          
          data[!data$keep, sprintf("CS.%s", pid)] = foreach(
            j = which(!data$keep), .combine = "c") %do%
          {
            tags = strsplit(data$tags[j], ',')[[1]]
            foreach(snp = tags, .combine = "union") %do%
            {
              if(!is.na(keep_CS[snp])) keep_CS[snp]
            } %>% paste(collapse = ',')
          }
          data[which(data[, sprintf("CS.%s", pid)] == ''), 
               sprintf("CS.%s", pid)] = NA
        }
      }
      
      filter(data, !is.na(CS.EUR) | !is.na(CS.AFR) | !is.na(CS.EAS))
    }
    
    write_delim(variants_CS95, gzfile(sprintf(
      "results/%s/%s/%s/variants_CS95.txt.gz", Dir, pheno, db)), delim = '\t')
    
    CS_overlap = foreach(l = which(rowSums(results[,c(15,21,27)] > 0) > 1),
                         .combine = "rbind") %dopar%
    {
      i = results$chromosome[l]
      block = results$block[l]

      data = read.delim(sprintf("%s/%s/%s/chr%d/%d/variant_info.txt.gz",
                                Dir, pheno, db, i, block))

      Ks = which(results[l,c(15,21,27)] > 0)
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

    write_delim(CS_overlap, sprintf("results/%s/%s/%s/CS_overlap.txt", Dir,
                                    pheno, db), delim = '\t')
    
    variants_0.9 = foreach(l = which(
      results$num_variants_0.5.EUR > 0 | results$num_variants_0.5.AFR > 0 |
        results$num_variants_0.5.EAS > 0), .combine = "rbind") %dopar%
    {
      i = results$chromosome[l]
      block = results$block[l]

      data = read.delim(sprintf("%s/%s/%s/chr%d/%d/variant_info.txt.gz",
                                Dir, pheno, db, i, block))

      filter(data, CL.EUR >= 0.9 | CL.AFR >= 0.9 | CL.EAS >= 0.9)
    }

    write_delim(variants_0.9, gzfile(sprintf(
      "results/%s/%s/%s/variants_0.9.txt.gz", Dir, pheno, db)), delim = '\t')
  }
}

# joint analysis
for(pheno in phenos)
{
  dir.create(sprintf("results/%s/%s/joint", Dir, pheno), recursive = T)

  results = foreach(db = dbs) %do%
  {
    read.delim(sprintf("results/%s/%s/%s/results.txt", Dir, pheno, db))
  }

  CS_common = foreach(l = 1:nrow(results[[1]]), .combine = "rbind") %dopar%
  {
    i = results[[1]]$chromosome[l]
    block = results[[1]]$block[l]

    if(rowSums(results[[1]][l,c(15,21,27)]) > 0 & rowSums(results[[2]][l,c(15,21,27)]) > 0)
    {
      data = foreach(db = dbs) %do%
      {
        read.delim(sprintf("%s/%s/%s/chr%d/%d/variant_info.txt.gz",
                           Dir, pheno, db, i, block))
      }

      res = data.frame("pheno" = pheno, "chromosome" = i, "block" = block)
      for(pid in pids)
      {
        res[[sprintf("UKBB_nCS.%s",pid)]] = results[[1]][l,sprintf("num_CS95.%s",pid)]
        res[[sprintf("GLGC_nCS.%s",pid)]] = results[[2]][l,sprintf("num_CS95.%s",pid)]
        res[[sprintf("UKBB_in_GLGC.%s",pid)]] = length(unique(na.omit(
          data[[1]][,sprintf("CS.%s",pid)][!is.na(data[[2]][,sprintf("CS.%s",pid)])])))
        res[[sprintf("GLGC_in_UKBB.%s",pid)]] = length(unique(na.omit(
          data[[2]][,sprintf("CS.%s",pid)][!is.na(data[[1]][,sprintf("CS.%s",pid)])])))
      }
    } else {
      res = data.frame("pheno" = pheno, "chromosome" = i, "block" = block)
      for(pid in pids)
      {
        res[[sprintf("UKBB_nCS.%s",pid)]] = results[[1]][l,sprintf("num_CS95.%s",pid)]
        res[[sprintf("GLGC_nCS.%s",pid)]] = results[[2]][l,sprintf("num_CS95.%s",pid)]
        res[[sprintf("UKBB_in_GLGC.%s",pid)]] = 0
        res[[sprintf("GLGC_in_UKBB.%s",pid)]] = 0
      }
    }
    res
  }

  write_delim(CS_common, sprintf("results/%s/%s/joint/CS_common_num.txt",
                                 Dir, pheno), delim = '\t')
  
  variant_0.9_UKBB = foreach(l = 1:nrow(results[[1]]), .combine = "rbind") %dopar%
  {
    i = results[[1]]$chromosome[l]
    block = results[[1]]$block[l]
    
    data = foreach(db = dbs) %do%
    {
      read.delim(sprintf("%s/%s/%s/chr%d/%d/variant_info.txt.gz",
                         Dir, pheno, db, i, block))
    }
    
    res = data.frame("pheno" = pheno, "chromosome" = i, "block" = block)
    for(pid in pids)
    {
      res[,sprintf("num_variants_0.9.UKBB.%s", pid)] = 
        sum(data[[1]][,sprintf("CL.%s", pid)] >= 0.9, na.rm = T)
      res[,sprintf("num_variants_0.9.UKBB_0.1.GLGC.%s", pid)] = sum(
        data[[1]][,sprintf("CL.%s", pid)] >= 0.9 & 
          data[[2]][,sprintf("CL.%s", pid)] < 0.1, na.rm = T)
    }
    
    res
  }
  
  print(sum(variant_0.9_UKBB$num_variants_0.9.UKBB.EUR))
  print(sum(variant_0.9_UKBB$num_variants_0.9.UKBB_0.1.GLGC.EUR))
  print(sum(variant_0.9_UKBB$num_variants_0.9.UKBB.AFR))
  print(sum(variant_0.9_UKBB$num_variants_0.9.UKBB_0.1.GLGC.AFR))
  print(sum(variant_0.9_UKBB$num_variants_0.9.UKBB.EAS))
  print(sum(variant_0.9_UKBB$num_variants_0.9.UKBB_0.1.GLGC.EAS))
  
  write_delim(variant_0.9_UKBB, sprintf("results/%s/%s/joint/variant_0.9_UKBB.txt",
                                        Dir, pheno), delim = '\t')
}
