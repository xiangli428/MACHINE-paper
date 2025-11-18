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
library(susieR)

setwd("real_data/scz2022")

pids = c("EUR", "EAS")
N = c("EUR" = 130644, "EAS" = 30761)

min_p = read.delim("data/min_p_sub.txt")[,-3]
min_p %<>% pivot_longer(-c(chromosome,block), names_to = c(".value","pid"),
                        names_pattern = "([A-Za-z]+_?[A-Za-z]+)\\.(.*)")

registerDoParallel(42)

Dir = "SuSiE"

foreach(l = 1:nrow(min_p_sub), .combine = "c", .inorder = F) %dopar%
{
  pid = min_p_sub$pid[l]
  i = min_p_sub$chromosome[l]
  block = min_p_sub$block[l]
  
  dir.create(sprintf("%s/%s/chr%d/%d", Dir, pid, i, block), recursive = T)
  
  data = read.delim(gzfile(sprintf("data/%s/chr%d/block_%d.txt.gz", pid, i, block)))
  data = data[, c(1:7,22,15)]
  
  R = readMM(gzfile(sprintf("LD/%s/chr%d/%d/LD.mtx.gz", pid, i, block))) %>%
    as.matrix()
  diag(R) = 1
  
  tic = Sys.time()
  fitted_susie = susie_rss(data$BETA / data$SE, R, N[pid])
  toc = Sys.time()
  
  saveRDS(fitted_susie, file = sprintf("%s/%s/chr%d/%d/fitted_susie.RData", 
                                       Dir, pid, i, block))
  data$PIP = fitted_susie$pip
  data$beta_mean = colSums(fitted_susie$mu * fitted_susie$alpha)
  
  CS_info = data.frame("CS" = "",
                       "CS.size" = 0,
                       "CS.minp" = 0,
                       "CS.purity" = 0,
                       "CS.SNP" = "")[-1,]
  
  data$CS = NA
  if(length(fitted_susie$sets$cs) > 0)
  {
    for(s in 1:length(fitted_susie$sets$cs))
    {
      CS_id = sprintf(
        "CS:%s-%d-%d-%s", pid, i, block, names(fitted_susie$sets$cs)[s])
      set = fitted_susie$sets$cs[[s]]
      data$CS[set] = CS_id
      CS_info %<>% rbind(data.frame(
        "CS" = CS_id,
        "CS.size" = length(set),
        "CS.minp" = min(data$PVAL[set]),
        "CS.purity" = fitted_susie$sets$purity[s,1],
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
    "num_variants" = nrow(data),
    "time" = difftime(toc, tic, units = "secs"),
    "num_variants_0.95" = sum(data$PIP >= 0.95),
    "num_variants_0.5" = sum(data$PIP >= 0.5),
    "num_CS95" = length(fitted_susie$sets$cs),
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

  CS1 = sort(unique(data[[1]]$CS[!is.na(data[[1]]$CS)]))

  foreach(set1 = CS1, .combine = "rbind") %do%
  {
    SNP1 = data[[1]]$rsid[which(data[[1]]$CS == set1)]
    sub = data[[2]][is.element(data[[2]]$rsid, SNP1),]
    if(any(!is.na(sub$CS)))
    {
      foreach(set2 = unique(sub$CS[!is.na(sub$CS)]), .combine = "rbind") %do%
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
