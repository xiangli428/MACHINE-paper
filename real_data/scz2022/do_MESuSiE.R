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

setwd("real_data/scz2022")

pids = c("EUR", "EAS")
N = c("EUR" = 130644, "EAS" = 30761)

min_p = read.delim("data/min_p_sub.txt")

registerDoParallel(44)

Dir = "MESuSiE"

foreach(l = 1:nrow(min_p), .combine = "c", .inorder = F) %dopar%
{
  i = min_p$chromosome[l]
  block = min_p$block[l]
  
  data = read.delim(sprintf("data/merge/chr%d/block_%d.txt.gz", i, block))
  data = data[,-c(8,12)]
  
  M = sum(!is.na(data$PVAL.EUR) & !is.na(data$PVAL.EAS))
  if(M > 0)
  {
    dir.create(sprintf("%s/chr%d/%d", Dir, i, block), recursive = T)
    
    data %<>% filter(!is.na(data$PVAL.EUR) & !is.na(data$PVAL.EAS))
    
    R = foreach(pid = pids) %do%
    {
      LD = readMM(gzfile(sprintf("LD/%s/chr%d/%d/LD.mtx.gz", pid, i, block))) %>%
        as.matrix()
      SNPs = read.delim(sprintf(
        "LD/%s/chr%d/%d/variant.txt", pid, i, block), header = F)[,1]
      rownames(LD) = colnames(LD) = SNPs
      diag(LD) = 1
      LD[data$rsid,data$rsid]
    }
    rm(LD)
    
    sumstat = foreach(pid = pids) %do%
    {
      tmp = data.frame("SNP" = data$rsid, "Beta" = 0, "Se" = 0)
      tmp[,2]= data[,sprintf("BETA.%s", pid)] / data[,sprintf("SE.%s", pid)] /
        sqrt(N[pid])
      tmp[is.na(tmp[,2]),2] = 0
      
      mutate(tmp, 
             Se = sqrt(1 / N[pid]),
             Z = Beta / Se,
             N = N[pid])
    }
    rm(tmp)
    names(R) = names(sumstat) = pids
    
    L = 5
    tic = Sys.time()
    MESuSiE_res = meSuSie_core(R, sumstat, L = L)
    while((length(MESuSiE_res$cs$cs) == L) & (L < 20))
    {
      L = L + 5
      MESuSiE_res = meSuSie_core(R, sumstat, L = L)
    }
    toc = Sys.time()
    
    saveRDS(MESuSiE_res, file = sprintf("%s/chr%d/%d/MESuSiE_res.RData", 
                                        Dir, i, block))
    
    data$PIP = MESuSiE_res$pip
    
    for(k in 1:2)
    {
      data[[sprintf("PIP.%s", pids[k])]] = 0
      data[[sprintf("mu.%s", pids[k])]] = 0
      data[sprintf("CS.%s", pids[k])] = NA
    }
    
    data$PIP.EUR = rowSums(MESuSiE_res$pip_config[,c(1,3)])
    data$PIP.EAS = rowSums(MESuSiE_res$pip_config[,c(2,3)])
    
    data$mu.EUR = foreach(l = 1:L, .combine = "+") %do%
    {
      MESuSiE_res$alpha[[l]][,1] * MESuSiE_res$mu1[[l]][[1]] +
        MESuSiE_res$alpha[[l]][,3] * MESuSiE_res$mu1[[l]][[3]][,1]
    }
    
    data$mu.EAS = foreach(l = 1:L, .combine = "+") %do%
    {
      MESuSiE_res$alpha[[l]][,2] * MESuSiE_res$mu1[[l]][[2]] +
        MESuSiE_res$alpha[[l]][,3] * MESuSiE_res$mu1[[l]][[3]][,2]
    }
    
    CS_info = data.frame("CS" = "",
                         "CS.size" = 0,
                         "CS.minp" = 0,
                         "CS.purity" = 0,
                         "CS.SNP" = "")[-1,]
    
    if(length(MESuSiE_res$cs$cs) > 0)
    {
      for(s in 1:length(MESuSiE_res$cs$cs))
      {
        set = MESuSiE_res$cs$cs[[s]]
        for(pid in strsplit(MESuSiE_res$cs$cs_category[s], '_')[[1]])
        {
          if(length(set) == 1)
          {
            if(!is.na(data[set, sprintf("PVAL.%s", pid)]))
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
              "CS.minp" = min(data[set, sprintf("PVAL.%s", pid)]),
              "CS.purity" = purity,
              "CS.SNP" = paste(data$rsid[set], collapse = ',')))
          }
        }
      }
    }
    
    write_delim(data, gzfile(sprintf(
      "%s/chr%d/%d/variant_info.txt.gz", Dir, i, block)), delim = '\t')
    write_delim(CS_info, sprintf("%s/chr%d/%d/CS_info.txt", Dir, i, block),
                delim = '\t')
    
    result = data.frame(
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
        data[,sprintf("PIP.%s", pid)] >= 0.95 & !is.na(data[,sprintf("PVAL.%s", pid)]))
      result[[sprintf("num_variants_0.5.%s", pid)]] = sum(
        data[,sprintf("PIP.%s", pid)] >= 0.5 & !is.na(data[,sprintf("PVAL.%s", pid)]))
      result[[sprintf("num_CS95.%s", pid)]] = 
        length(unique(na.omit(data[[sprintf("CS.%s", pid)]])))
      result[[sprintf("num_variants_CS95.%s", pid)]] = 
        sum(!is.na(data[,sprintf("CS.%s", pid)]))
    }
    
    write_delim(result, sprintf("%s/chr%d/%d/result.txt", Dir, i, block),
                delim = '\t')
  }
  
  NULL
}

dir.create(sprintf("results/%s", Dir), recursive = T)

# results
results = foreach(l = 1:nrow(min_p), .combine = "rbind") %dopar%
{
  i = min_p$chromosome[l]
  block = min_p$block[l]

  read.delim(sprintf("%s/chr%d/%d/result.txt", Dir, i, block))
}

write_delim(results, sprintf("results/%s/results.txt", Dir), delim = '\t')

# summary
summary = foreach(pid = pids, .combine = "rbind") %do%
{
  data.frame(
    "pid" = pid,
    "num_variants_0.95" = sum(results[,sprintf("num_variants_0.95.%s", pid)]),
    "num_variants_0.5" = sum(results[,sprintf("num_variants_0.5.%s", pid)]),
    "num_CS95" = sum(results[,sprintf("num_CS95.%s", pid)]),
    "num_variants_CS95" = sum(results[,sprintf("num_variants_CS95.%s", pid)]))
}

write_delim(summary, sprintf("results/%s/summary.txt", Dir), delim = '\t')

# CS info
CS_info_all = foreach(l = which(rowSums(results[,c(12,16)]) > 0),
                      .combine = "rbind") %dopar%
{
  i = results$chromosome[l]
  block = results$block[l]

  CS_info = read.delim(sprintf("%s/chr%d/%d/CS_info.txt", Dir, i, block))
  CS_info$chromosome = i
  CS_info$block = block
  CS_info[,c(6:7,1:5)]
}

write_delim(CS_info_all, sprintf("results/%s/CS_info.txt", Dir), delim = '\t')

# CS overlap
Ls = which(results$num_CS95.EUR > 0 & results$num_CS95.EAS > 0)
sum(results$num_CS95.EUR[Ls])
sum(results$num_CS95.EAS[Ls])

CS_overlap = foreach(l = Ls, .combine = "rbind") %dopar%
{
  i = min_p$chromosome[l]
  block = min_p$block[l]
  
  data = read.delim(gzfile(sprintf(
    "%s/chr%d/%d/variant_info.txt.gz", Dir, i, block)))
  
  foreach(s = 1:results$L[l], .combine = "rbind") %do%
  {
    CS_id_1 = sprintf("CS:%s-%d-%d-L%d", "EUR", i, block, s)
    CS_id_2 = sprintf("CS:%s-%d-%d-L%d", "EAS", i, block, s)
    
    if(is.element(CS_id_1, data$CS.EUR) & is.element(CS_id_2, data$CS.EAS))
    {
      set = which(data$CS.EUR == CS_id_1)
      data.frame("chromosome" = i,
                 "block" = block,
                 "CS.EUR" = CS_id_1,
                 "CS.EUR.minp" = min(data$PVAL.EUR[set]),
                 "CS.EAS" = CS_id_2,
                 "CS.EAS.minp" = min(data$PVAL.EAS[set]),
                 "CS.size" = length(set),
                 "CS.SNP" = paste(data$rsid[set], collapse = ','))
    }
  }
}

write_delim(CS_overlap, sprintf("results/%s/CS_overlap.txt", Dir), delim = '\t')

# CL >= 0.5
variants_0.5 = foreach(l = which(results$num_variants_0.5 > 0),
                      .combine = "rbind") %dopar%
{
  i = results$chromosome[l]
  block = results$block[l]

  data = read.delim(sprintf(
    "%s/chr%d/%d/variant_info.txt.gz", Dir, i, block))

  filter(data, PIP >= 0.5)
}

write_delim(variants_0.5, gzfile(sprintf("results/%s/variants_0.5.txt.gz", Dir)),
            '\t')

# CS variants
variants_CS95 = foreach(l = which(rowSums(results[,c(12,16)]) > 0),
                        .combine = "rbind") %dopar%
{
  i = results$chromosome[l]
  block = results$block[l]
  
  data = read.delim(sprintf(
    "%s/chr%d/%d/variant_info.txt.gz", Dir, i, block))
  
  filter(data, !is.na(CS.EUR) | !is.na(CS.EAS))
}

write_delim(variants_CS95, gzfile(sprintf("results/%s/variants_CS95.txt.gz", Dir)),
            '\t')
