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

setwd("real_data/scz2022")

pids = c("EUR", "EAS")
N = c("EUR" = 130644, "EAS" = 30761)

min_p = read.delim("data/min_p_sub.txt")

gLDSC_var = foreach(pid = pids) %do%
{
  readRDS(sprintf("gLDSC_results/%s/gLDSC_var.RData", pid)) %>% pmax(max(.) / 20)
}
names(gLDSC_var) = pids

registerDoParallel(16)

Dir = "MACHINE-gLDSC"

foreach(l = 1:nrow(min_p), .combine = "c", .inorder = F) %dopar%
{
  i = min_p$chromosome[l]
  block = min_p$block[l]
  
  dir.create(sprintf("%s/chr%d/%d", Dir, i, block), recursive = T)
  
  data = read_delim(gzfile(sprintf("data/merge/chr%d/block_%d.txt.gz", i, block)), 
                    delim = '\t', progress = F, show_col_types = F)
  data = data[,-c(8,12)]
  
  R = foreach(pid = pids) %do%
  {
    readMM(gzfile(sprintf("LD/%s/chr%d/%d/LD.mtx.gz", pid, i, block)))
  }
  R_eig = foreach(pid = pids) %do%
  {
    readRDS(sprintf("LD/%s/chr%d/%d/LD_eig.RData", pid, i, block))
  }
  
  Z = foreach(pid = pids, .combine = "cbind") %do%
  {
    z = data[,sprintf("BETA.%s", pid)] / data[,sprintf("SE.%s", pid)]
  }
  
  gLDSC_var_sub = cbind(gLDSC_var[[1]][data$rsid], gLDSC_var[[2]][data$rsid])
  for(k in 1:2)
  {
    gLDSC_var_sub[,k] = gLDSC_var_sub[,k] / mean(gLDSC_var_sub[,k], na.rm = T) * 0.005
  }
  a = rowMeans(gLDSC_var_sub, na.rm = T)
  a = a / mean(a) * 0.005
  c = gLDSC_var_sub / rowSums(gLDSC_var_sub, na.rm = T) * 
    rowSums(!is.na(gLDSC_var_sub), na.rm = T) * 0.2
  
  MACHINE = CreateMACHINEObject(Z, 
                                R, 
                                N,
                                SNP_ID = data$rsid,
                                POP_ID = pids,
                                in_sample_LD = F,
                                R_eig = R_eig,
                                a = a,
                                c = c)
  rm(R, R_eig)
  tic = Sys.time()
  MACHINE = MACHINE_MCMC(MACHINE, mcmc_n = 5500, burn_in = 500, thin = 2,
                         get_CS = F)
  while(max(MACHINE@mcmc_samples$PSRF_beta) > 1.2 & 
        MACHINE@mcmc_samples$n_burnin < 5500)
  {
    MACHINE = MACHINE_MCMC(MACHINE, mcmc_n = 1000, burn_in = 1000, thin = 2,
                           get_CS = F)
  }
  MACHINE@CS = MACHINE_CS(MACHINE)
  toc = Sys.time()
  
  saveRDS(MACHINE, file = sprintf("%s/chr%d/%d/MACHINE.RData", Dir, i, block))
  
  for(pid in pids)
  {
    data[sprintf("sigma2_mean.%s", pid)] = NA
    data[sprintf("CL.%s", pid)] = NA
    data[sprintf("beta_mean.%s", pid)] = NA
    data[sprintf("CS.%s", pid)] = NA
  }
  
  CS_info = data.frame("CS" = "",
                       "CS.size" = 0,
                       "CS.minp" = 0,
                       "CS.purity" = 0,
                       "CS.SNP" = "")[-1,]
  
  for(pid in pids)
  {
    idx = which(!is.na(data[[sprintf("PVAL.%s", pid)]]))
    data[sprintf("sigma2_mean.%s", pid)] = NA
    data[idx, sprintf("sigma2_mean.%s", pid)] = 
      MACHINE@mcmc_mean$sigma2[idx, pid]
    data[sprintf("CL.%s", pid)] = NA
    data[idx, sprintf("CL.%s", pid)] = 
      MACHINE@CL[idx, pid]
    data[sprintf("beta_mean.%s", pid)] = NA
    data[idx, sprintf("beta_mean.%s", pid)] = 
      MACHINE@mcmc_mean$beta[idx, pid]
    
    data[sprintf("CS.%s", pid)] = NA
    if(length(MACHINE@CS[[pid]]$sets) > 0)
    {
      for(s in 1:length(MACHINE@CS[[pid]]$sets))
      {
        CS_id = sprintf("CS:%s-%d-%d-%s", pid, i, block, s)
        set = MACHINE@CS[[pid]]$sets[[s]]
        data[set, sprintf("CS.%s", pid)] = CS_id
        CS_info %<>% rbind(data.frame(
          "CS" = CS_id,
          "CS.size" = length(set),
          "CS.minp" = min(data[set, sprintf("PVAL.%s", pid)]),
          "CS.purity" = MACHINE@CS[[pid]]$purity[s,1],
          "CS.SNP" = paste(data$rsid[set], collapse = ',')))
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
    "num_variants" = MACHINE@M,
    "max_PSRF" = max(MACHINE@mcmc_samples$PSRF_beta),
    "max_beta" = lapply(MACHINE@mcmc_samples$beta, abs) %>% sapply(max) %>% max,
    "max_betamean" = max(abs(MACHINE@mcmc_mean$beta)),
    "mcmc_n" = MACHINE@mcmc_samples$n_samples + MACHINE@mcmc_samples$n_burnin,
    "time" = difftime(toc, tic, units = "secs"))
  
  for(pid in pids)
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
  }
  
  write_delim(result, sprintf("%s/chr%d/%d/result.txt", Dir, i, block),
              delim = '\t')
  
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
CS_info_all = foreach(l = which(rowSums(results[,c(13,19)]) > 0),
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
Ls = which(results[,13] > 0 & results[,19] > 0)
sum(results$num_CS95.EUR[Ls])
sum(results$num_CS95.EAS[Ls])

CS_overlap = foreach(l = Ls, .combine = "rbind") %dopar%
{
  i = min_p$chromosome[l]
  block = min_p$block[l]

  data = read.delim(gzfile(sprintf(
    "%s/chr%d/%d/variant_info.txt.gz", Dir, i, block)))

  CS1 = sort(unique(data$CS.EUR[!is.na(data$CS.EUR)]))

  foreach(set1 = CS1, .combine = "rbind") %do%
  {
    idx1 = which(data$CS.EUR == set1)
    SNP1 = data$rsid[idx1]
    if(any(!is.na(data$CS.EAS[idx1])))
    {
      foreach(set2 = unique(na.omit(data$CS.EAS[idx1])), .combine = "rbind") %do%
      {
        idx2 = which(data$CS.EAS == set2)
        SNP2 = data$rsid[idx2]
        data.frame("chromosome" = i,
                   "block" = block,
                   "CS.EUR" = set1,
                   "CS.EUR.size" = length(SNP1),
                   "CS.EUR.minp" = min(data$PVAL.EUR[idx1]),
                   "CS.EUR.SNP" = paste(SNP1, collapse = ','),
                   "CS.EAS" = set2,
                   "CS.EAS.size" = length(SNP2),
                   "CS.EAS.minp" = min(data$PVAL.EAS[idx2]),
                   "CS.EAS.SNP" = paste(SNP2, collapse = ','),
                   "intersect.size" = length(intersect(SNP1,SNP2)),
                   "intersect.SNP" = paste(intersect(SNP1,SNP2),
                                           collapse = ','))
      }
    }
  }
}

write_delim(CS_overlap, sprintf("results/%s/CS_overlap.txt", Dir), delim = '\t')

# CL >= 0.5
variants_0.5 = foreach(l = which(rowSums(results[,c(12,18)]) > 0),
                      .combine = "rbind") %dopar%
{
  i = results$chromosome[l]
  block = results$block[l]

  data = read.delim(gzfile(sprintf(
    "%s/chr%d/%d/variant_info.txt.gz", Dir, i, block)))

  filter(data, CL.EUR >= 0.5 | CL.EAS >= 0.5)
}

write_delim(variants_0.5, gzfile(sprintf("results/%s/variants_0.5.txt.gz", Dir)),
            '\t')

# CS variants
variants_CS95 = foreach(l = which(rowSums(results[,c(13,19)]) > 0),
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
