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
library(MESuSiE)

setwd("real_data/glgc")

dbs = c("UKBB", "GLGC")
phenos = c("HDL", "LDL", "TG", "TC")
pids = c("EUR", "AFR", "EAS")

registerDoParallel(32)

Dir = "MESuSiE"

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
    
    foreach(l = 1:nrow(min_p), .combine = "c", .inorder = F) %dopar%
    {
      i = min_p$chromosome[l]
      block = min_p$block[l]
      
      dir.create(sprintf("%s/%s/%s/chr%d/%d", Dir, pheno, db, i, block), 
                 recursive = T)
      
      data = read.delim(gzfile(sprintf("data/%s/%s/merge/chr%d/block_%d.txt.gz", 
                                       db, pheno, i, block)))
      data = data[data$keep,-20]
      data$tags[is.na(data$tags)] = ""
      
      used_pids = foreach(pid = pids, .combine = "c") %do%
      {
        if(!all(is.na(data[,sprintf("neglog10_pval.%s", pid)]))) pid
      }
      
      R = foreach(pid = used_pids) %do%
      {
        LD = readMM(sprintf("LD/%s/%s/chr%d/%d/LD.mtx.gz", pheno, pid, i, block)) %>%
          as.matrix()
        SNPs = data$rsid[!is.na(data[,sprintf("neglog10_pval.%s", pid)])]
        rownames(LD) = colnames(LD) = SNPs
        diag(LD) = 1
        LD
      }
      
      if("EAS" %in% used_pids)
      {
        data %<>% filter(!is.na(neglog10_pval.EUR) & !is.na(neglog10_pval.AFR) &
                           !is.na(neglog10_pval.EAS))
      } else {
        data %<>% filter(!is.na(neglog10_pval.EUR) & !is.na(neglog10_pval.AFR))
      }
      
      M = nrow(data)
      
      if(db == "UKBB")
      {
        sumstat = foreach(pid = used_pids) %do%
        {
          tmp = data.frame("SNP" = data$rsid, "Beta" = 0)
          tmp[,2] = data[,sprintf("beta.%s", pid)] / data[,sprintf("se.%s", pid)] / 
            sqrt(N[pid, pheno])
          mutate(tmp, Se = sqrt(1 / N[pid, pheno]), Z = Beta / Se, 
                 N = N[pid, pheno])
        }
      } else {
        sumstat = foreach(pid = used_pids) %do%
        {
          tmp = data.frame("SNP" = data$rsid, "Beta" = 0)
          tmp[,2] = data[,sprintf("Z.%s", pid)] / sqrt(N[pid, pheno])
          mutate(tmp, Se = sqrt(1 / N[pid, pheno]), Z = Beta / Se, 
                 N = N[pid, pheno])
        }
      }
      rm(tmp)
      names(R) = names(sumstat) = used_pids
      
      for(pid in used_pids)
      {
        R[[pid]] = R[[pid]][data$rsid,data$rsid]
      }
      
      L = 5
      tic = Sys.time()
      MESuSiE_res = meSuSie_core(R, sumstat, L = L)
      while((length(MESuSiE_res$cs$cs) == L) & (L < 20))
      {
        L = L + 5
        MESuSiE_res = meSuSie_core(R, sumstat, L = L)
      }
      toc = Sys.time()
      
      saveRDS(MESuSiE_res, file = sprintf(
        "%s/%s/%s/chr%d/%d/MESuSiE_res.RData", Dir, pheno, db, i, block))
      
      data$PIP = MESuSiE_res$pip
      
      for(k in 1:3)
      {
        data[[sprintf("PIP.%s", pids[k])]] = 0
        data[[sprintf("mu.%s", pids[k])]] = 0
        data[sprintf("CS.%s", pids[k])] = NA
      }
      
      for(j in 1:ncol(MESuSiE_res$pip_config))
      {
        col_pids = strsplit(colnames(MESuSiE_res$pip_config)[j], "_")[[1]]
        for(k in 1:length(col_pids))
        {
          data[[sprintf("PIP.%s", col_pids[k])]] %<>% add(MESuSiE_res$pip_config[,j])
          data[[sprintf("mu.%s", col_pids[k])]] %<>% add(
            foreach(ll = 1:L, .combine = "+") %do%
            {
              if(length(col_pids) == 1)
              {
                MESuSiE_res$alpha[[ll]][,j] * MESuSiE_res$mu1[[ll]][[j]]
              } else {
                MESuSiE_res$alpha[[ll]][,j] * MESuSiE_res$mu1[[ll]][[j]][,k]
              }
            }
          )
        }
      }
      
      CS_info = data.frame("CS" = "",
                           "CS.size" = 0,
                           "CS.nlog10minp" = 0,
                           "CS.purity" = 0,
                           "CS.SNP" = "",
                           "CS.tags" = "")[-1,]
      
      n_CS = length(MESuSiE_res$cs$cs)
      if(n_CS > 0)
      {
        for(s in 1:n_CS)
        {
          set = MESuSiE_res$cs$cs[[s]]
          tags = foreach(j = set, .combine = "union") %do%
          {
            strsplit(data$tags[j], ",")[[1]]
          }
          for(pid in strsplit(MESuSiE_res$cs$cs_category[s], '_')[[1]])
          {
            if(length(set) == 1)
            {
              if(!is.na(data[set, sprintf("neglog10_pval.%s", pid)]))
              {
                purity = 1
              } else {
                purity = 0
              }
            } else {
              R_sub = R[[pid]][set, set]
              purity = min(abs(R_sub[upper.tri(R_sub)]))
            }
            
            if(purity >= 0.5)
            {
              CS_id = sprintf("CS:%s-%d-%d-%s", pid, i, block, 
                              names(MESuSiE_res$cs$cs)[s])
              data[set, sprintf("CS.%s", pid)] = CS_id
              CS_info %<>% rbind(data.frame(
                "CS" = CS_id,
                "CS.size" = length(set),
                "CS.nlog10minp" = max(data[set, sprintf("neglog10_pval.%s", pid)]),
                "CS.purity" = purity,
                "CS.SNP" = paste(data$rsid[set], collapse = ','),
                "CS.tags" = paste(tags, collapse = ',')))
            }
          }
        }
      }
      
      write_delim(data, gzfile(sprintf("%s/%s/%s/chr%d/%d/variant_info.txt.gz", 
                                       Dir, pheno, db, i, block)), 
                  delim = '\t')
      write_delim(CS_info, sprintf("%s/%s/%s/chr%d/%d/CS_info.txt", Dir, 
                                   pheno, db, i, block), delim = '\t')
      
      result = data.frame(
        "pheno" = pheno,
        "db" = db,
        "chromosome" = i,
        "block" = block,
        "num_variants" = M,
        "L" = L,
        "time" = difftime(toc, tic, units = "secs"),
        "num_variants_0.95" = sum(data$PIP >= 0.95),
        "num_variants_0.5" = sum(data$PIP >= 0.5),
        "num_CS95" = length(MESuSiE_res$cs$cs),
        "num_variants_CS95" = if(length(MESuSiE_res$cs$cs) > 0) sum(
          sapply(MESuSiE_res$cs$cs, length)) else 0)
      for(pid in pids)
      {
        result[[sprintf("num_variants_0.95.%s", pid)]] = sum(
          data[,sprintf("PIP.%s", pid)] >= 0.95 & 
            !is.na(data[,sprintf("neglog10_pval.%s", pid)]))
        result[[sprintf("num_variants_0.5.%s", pid)]] = sum(
          data[,sprintf("PIP.%s", pid)] >= 0.5 & 
            !is.na(data[,sprintf("neglog10_pval.%s", pid)]))
        result[[sprintf("num_CS95.%s", pid)]] = length(unique(na.omit(
          data[,sprintf("CS.%s", pid)])))
        result[[sprintf("num_variants_CS95.%s", pid)]] = 
          sum(!is.na(data[,sprintf("CS.%s", pid)]))
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

    CS_info_all = foreach(l = which(results$num_CS95 > 0),
                          .combine = "rbind") %dopar%
    {
      i = results$chromosome[l]
      block = results$block[l]

      if(file.exists(sprintf("%s/%s/%s/chr%d/%d/CS_info.txt",
                             Dir, pheno, db, i, block)))
        read.delim(sprintf("%s/%s/%s/chr%d/%d/CS_info.txt",
                           Dir, pheno, db, i, block))
    }
    CS_info_all$CS.tags[is.na(CS_info_all$CS.tags)] = ""
    write_delim(CS_info_all, sprintf("results/%s/%s/%s/CS_info.txt", Dir,
                                     pheno, db), delim = '\t')

    variants_CS95 = foreach(l = which(results$num_CS95 > 0),
                            .combine = "rbind") %dopar%
    {
      i = results$chromosome[l]
      block = results$block[l]

      data = read.delim(sprintf("data/%s/%s/merge/chr%d/block_%d.txt.gz",
                                db, pheno, i, block))
      used_pids = foreach(pid = pids, .combine = "c") %do%
      {
        if(!all(is.na(data[,sprintf("neglog10_pval.%s", pid)]))) pid
      }
      if("EAS" %in% used_pids)
      {
        data %<>% filter(!is.na(neglog10_pval.EUR) & !is.na(neglog10_pval.AFR) &
                           !is.na(neglog10_pval.EAS))
      } else {
        data %<>% filter(!is.na(neglog10_pval.EUR) & !is.na(neglog10_pval.AFR))
      }

      data_keep = read.delim(sprintf("%s/%s/%s/chr%d/%d/variant_info.txt.gz",
                                     Dir, pheno, db, i, block))

      data$PIP = NA
      for(pid in pids)
      {
        data[sprintf("PIP.%s", pid)] = NA
        data[sprintf("mu.%s", pid)] = NA
        data[sprintf("CS.%s", pid)] = NA
      }
      data[data$keep, 22:31] = data_keep[,21:30]

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

    CS_overlap = foreach(l = which(rowSums(results[,c(14,18,22)] > 0) > 1),
                         .combine = "rbind") %dopar%
    {
      i = results$chromosome[l]
      block = results$block[l]

      data = read.delim(sprintf("%s/%s/%s/chr%d/%d/variant_info.txt.gz",
                                Dir, pheno, db, i, block))

      foreach(s = 1:results$L[l], .combine = "rbind") %do%
      {
        Ks = foreach(k = 1:3, .combine = "c") %do%
        {
          if(is.element(sprintf("CS:%s-%d-%d-L%d", pids[k], i, block, s),
                        data[,sprintf("CS.%s",pids[k])]))
          {
            k
          }
        }

        if(length(Ks) == 2)
        {
          set = which(data[,sprintf("CS.%s",pids[Ks[1]])] ==
                        sprintf("CS:%s-%d-%d-L%d", pids[Ks[1]], i, block, s))

          tmp = data.frame("pheno" = pheno,
                           "db" = db,
                           "chromosome" = i,
                           "block" = block,
                           "CS.EUR" = "",
                           "CS.AFR" = "",
                           "CS.EAS" = "",
                           "size" = length(set),
                           "SNP" = paste(data$rsid[set], collapse = ','))
          tmp[,sprintf("CS.%s", pids[Ks[1]])] =
            sprintf("CS:%s-%d-%d-L%d", pids[Ks[1]], i, block, s)
          tmp[,sprintf("CS.%s", pids[Ks[2]])] =
            sprintf("CS:%s-%d-%d-L%d", pids[Ks[2]], i, block, s)
          tmp
        } else if(length(Ks) == 3){
          set = which(data[,sprintf("CS.%s",pids[1])] ==
                        sprintf("CS:%s-%d-%d-L%d", pids[1], i, block, s))

          data.frame("pheno" = pheno,
                     "db" = db,
                     "chromosome" = i,
                     "block" = block,
                     "CS.EUR" = sprintf("CS:%s-%d-%d-L%d", pids[1], i, block, s),
                     "CS.AFR" = sprintf("CS:%s-%d-%d-L%d", pids[2], i, block, s),
                     "CS.EAS" = sprintf("CS:%s-%d-%d-L%d", pids[3], i, block, s),
                     "size" = length(set),
                     "SNP" = paste(data$rsid[set], collapse = ','))
        }
      }
    }

    write_delim(CS_overlap, sprintf("results/%s/%s/%s/CS_overlap.txt", Dir,
                                    pheno, db), delim = '\t')

    variants_0.9 = foreach(l = which(results$num_variants_0.5 > 0),
                           .combine = "rbind") %dopar%
    {
      i = results$chromosome[l]
      block = results$block[l]

      data = read.delim(sprintf("%s/%s/%s/chr%d/%d/variant_info.txt.gz",
                                Dir, pheno, db, i, block))

      filter(data, PIP >= 0.9)
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

    if(results[[1]][l,10] > 0 & results[[2]][l,10] > 0)
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
        sum(data[[1]][,sprintf("PIP.%s", pid)] >= 0.9, na.rm = T)
      res[,sprintf("num_variants_0.9.UKBB_0.1.GLGC.%s", pid)] = sum(
        data[[1]][,sprintf("PIP.%s", pid)] >= 0.9 &
          data[[2]][,sprintf("PIP.%s", pid)] < 0.1, na.rm = T)
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
