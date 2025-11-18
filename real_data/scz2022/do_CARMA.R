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
library(CARMA)

setwd("real_data/scz2022")

pids = c("EUR", "EAS")
N = c("EUR" = 130644, "EAS" = 30761)

min_p = read.delim("data/min_p_sub.txt")[,-3]
min_p %<>% pivot_longer(-c(chromosome,block), names_to = c(".value","pid"),
                        names_pattern = "([A-Za-z]+_?[A-Za-z]+)\\.(.*)")

registerDoParallel(22)

Dir = "CARMA"

foreach(l = 1:nrow(min_p), .combine = "c", .inorder = F) %dopar%
{
  pid = min_p$pid[l]
  i = min_p$chromosome[l]
  block = min_p$block[l]
  
  dir.create(sprintf("%s/%s/chr%d/%d", Dir, pid, i, block), recursive = T)
  
  data = read.delim(sprintf("data/%s/chr%d/block_%d.txt.gz", pid, i, block))
  data = data[, c(1:7,13:15)]
  
  R = readMM(sprintf("LD/%s/chr%d/%d/LD.mtx.gz", pid, i, block)) %>%
    as.matrix()
  diag(R) = 1
  
  M = nrow(data)
  n = N[pid]
  
  z = data$BETA / data$SE
  
  tic = Sys.time()
  CARMA_results <- CARMA(
    list(z), list(R), lambda.list = list(1),
    rho.index = 0.95, outlier.switch = T,
    output.labels = sprintf("%s/%s/chr%d/%d", Dir, pid, i, block))
  toc = Sys.time()
  
  saveRDS(CARMA_results, file = sprintf(
    "%s/%s/chr%d/%d/CARMA_results.RData", Dir, pid, i, block))
  
  data$outlier = F
  data$outlier[CARMA_results[[1]]$Outliers$Index] = T
  data$PIP = CARMA_results[[1]]$PIPs
  
  CS_info = data.frame("CS" = "",
                       "CS.size" = 0,
                       "CS.minp" = 0,
                       "CS.purity" = 0,
                       "CS.SNP" = "")[-1,]
  
  data$CS = NA
  n_CS = length(CARMA_results[[1]][[2]][[2]])
  if(n_CS > 0)
  {
    for(s in 1:n_CS)
    {
      CS_id = sprintf("CS:%s-%d-%d-%s", pid, i, block, s)
      set = CARMA_results[[1]][[2]][[2]][[s]]
      data$CS[set] = CS_id
      CS_info %<>% rbind(data.frame(
        "CS" = CS_id,
        "CS.size" = length(set),
        "CS.minp" = min(data$PVAL[set]),
        "CS.purity" = min(abs(R[set,set])),
        "CS.SNP" = paste(data$rsid[set], collapse = ',')))
    }
  }
  
  write_delim(data, gzfile(sprintf(
    "%s/%s/chr%d/%d/variant_info.txt.gz", Dir, pid, i, block)), delim = '\t')
  write_delim(CS_info, sprintf("%s/%s/chr%d/%d/CS_info.txt", Dir, pid, i, block),
              delim = '\t')
  
  result = data.frame(
    "pid" = pid,
    "chromosome" = i,
    "block" = block,
    "num_variants" = M,
    "time" = difftime(toc, tic, units = "secs"),
    "num_variants_0.95" = sum(data$PIP >= 0.95),
    "num_variants_0.5" = sum(data$PIP >= 0.5),
    "num_CS95" = n_CS,
    "num_variants_CS95" = sum(!is.na(data$CS)))
  
  write_delim(result, sprintf("%s/%s/chr%d/%d/result.txt", Dir, pid, i, block),
              delim = '\t')
  
  NULL
  
}

min_p = read.delim("data/min_p_sub.txt")[,-3]
dir.create(sprintf("results/%s/EUR", Dir), recursive = T)
dir.create(sprintf("results/%s/EAS", Dir), recursive = T)

for(pid in pids)
{
  # results
  results = foreach(l = 1:nrow(min_p), .combine = "rbind") %dopar%
  {
    i = min_p$chromosome[l]
    block = min_p$block[l]
    
    if(file.exists(sprintf("%s/%s/chr%d/%d/result.txt", Dir, pid, i, block)))
      read.delim(sprintf("%s/%s/chr%d/%d/result.txt", Dir, pid, i, block))
  }
  
  write_delim(results, sprintf("results/%s/%s/results.txt", Dir, pid), delim = '\t')
  
  # CS info
  CS_info_all = foreach(l = which(results$num_CS95 > 0), .combine = "rbind") %dopar%
  {
    i = results$chromosome[l]
    block = results$block[l]
    
    CS_info = read.delim(sprintf("%s/%s/chr%d/%d/CS_info.txt", Dir, pid, i, block))
    CS_info$pid = pid
    CS_info$chromosome = i
    CS_info$block = block
    CS_info[,c(6:8,1:5)]
  }
  
  write_delim(CS_info_all, sprintf("results/%s/%s/CS_info.txt", Dir, pid), 
              delim = '\t')
  
  # CL >= 0.5
  variants_0.5 = foreach(l = which(results$num_variants_0.5 > 0),
                         .combine = "rbind") %dopar%
  {
    i = results$chromosome[l]
    block = results$block[l]
    
    data = read.delim(sprintf(
      "%s/%s/chr%d/%d/variant_info.txt.gz", Dir, pid, i, block))
    
    filter(data, PIP >= 0.5)
  }
  
  if(!is.null(variants_0.5))
  {
    write_delim(variants_0.5, gzfile(sprintf(
      "results/%s/%s/variants_0.5.txt.gz", Dir, pid)), '\t')
  }
  
  # CS variants
  variants_CS95 = foreach(l = which(results$num_CS95 > 0),
                         .combine = "rbind") %dopar%
  {
    i = results$chromosome[l]
    block = results$block[l]
    
    data = read.delim(sprintf(
      "%s/%s/chr%d/%d/variant_info.txt.gz", Dir, pid, i, block))
    
    filter(data, !is.na(CS))
  }
  
  write_delim(variants_CS95, gzfile(sprintf(
    "results/%s/%s/variants_CS95.txt.gz", Dir, pid)), '\t')
}

results = foreach(pid = pids) %do%
{
  read.delim(sprintf("results/%s/%s/results.txt", Dir, pid))
}
names(results) = pids

# summary
summary = foreach(pid = pids, .combine = "rbind") %do%
{
  data.frame("pid" = pid,
             "num_variants_0.95" = sum(results[[pid]]$num_variants_0.95),
             "num_variants_0.5" = sum(results[[pid]]$num_variants_0.5),
             "num_CS95" = sum(results[[pid]]$num_CS95),
             "num_variants_CS95" = sum(results[[pid]]$num_variants_CS95))
}
write_delim(summary, sprintf("results/%s/summary.txt", Dir), delim = '\t')

# CS overlap
n_CS = results[[1]][,2:3]
n_CS$EUR = results[[1]]$num_CS95
n_CS$EAS = results[[2]]$num_CS95
Ls = which(rowSums(n_CS[,3:4] > 0) > 1)

CS_overlap = foreach(l = which(rowSums(n_CS[,3:4] > 0) > 1),
                     .combine = "rbind") %dopar%
{
  i = n_CS$chromosome[l]
  block = n_CS$block[l]
  
  data = foreach(pid = pids) %do%
  {
    read.delim(gzfile(sprintf(
      "%s/%s/chr%d/%d/variant_info.txt.gz", Dir, pid, i, block)))
  }
  
  CS1 = sort(unique(na.omit(data[[1]]$CS)))
  
  foreach(set1 = CS1, .combine = "rbind") %do%
  {
    SNP1 = data[[1]]$rsid[which(data[[1]]$CS == set1)]
    sub = data[[2]][is.element(data[[2]]$rsid, SNP1),]
    if(any(!is.na(sub$CS)))
    {
      foreach(set2 = unique(na.omit(sub$CS)), .combine = "rbind") %do%
      {
        SNP2 = data[[2]]$rsid[which(data[[2]]$CS == set2)]
        data.frame("chromosome" = i,
                   "block" = block,
                   "CS.EUR" = set1,
                   "CS.EUR.size" = length(SNP1),
                   "CS.EUR.minp" = min(data[[1]]$PVAL[
                     which(data[[1]]$CS == set1)]),
                   "CS.EUR.SNP" = paste(SNP1, collapse = ','),
                   "CS.EAS" = set2,
                   "CS.EAS.size" = length(SNP2),
                   "CS.EAS.minp" = min(data[[2]]$PVAL[
                     which(data[[2]]$CS == set2)]),
                   "CS.EAS.SNP" = paste(SNP2, collapse = ','),
                   "intersect.size" = length(intersect(SNP1,SNP2)),
                   "intersect.SNP" = paste(intersect(SNP1,SNP2),
                                           collapse = ','))
      }
    }
  }
}

write_delim(CS_overlap, sprintf("results/%s/CS_overlap.txt", Dir), delim = '\t')
