options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)
library(mvtnorm)

setwd("simulation/data")

# gLDSC var
pops = c("EUR" = 1, "EAS" = 5)

registerDoParallel(10)

select_block = readRDS("select_block.RData")

variant = foreach(block = select_block, .combine = "rbind") %dopar%
{
  variant_pop = foreach(pid = names(pops)) %do%
  {
    read.delim(gzfile(sprintf(
      "%s/%s/%s/LD/chr%s/%s/variant_reblock.txt.gz",
      "~/Documents/GWAS/data/UKBB", pops[pid], 
      "imputed_genotype_unique_info-0.8_maf-0.001_hwe-1e-6", 
      i, block)))
  }

  merge(variant_pop[[1]], variant_pop[[2]], sort = F)
}
variant = variant[,-c(1,5)]

gLDSC_res = foreach(pid = names(pops)) %do%
{
  readRDS(sprintf("%s/gLDSC_result/%s/gldsc_res_pos.Rdata", 
                  "../../real_data/scz2022", pid))
}
names(gLDSC_res) = names(pops)

for(k in 1:2)
{
  names(gLDSC_res[[k]][[1]]) = rownames(gLDSC_res[[k]][[2]])
}

for(pid in names(pops))
{
  annos = rownames(gLDSC_res[[pid]][[2]])[-1]
  variant[[sprintf("var.%s", pid)]] = foreach(
    anno = annos, .combine = "+") %dopar%
    {
      bed = read_delim(sprintf(
        "../../annotations/baseline_bed_intersect/%s/%s.bed",
        pops[pid], anno), delim = '\t', col_names = F)
      snps = intersect(variant$rsid, bed$X4)
      s = rep(0, nrow(variant))
      names(s) = variant$rsid
      s[intersect(variant$rsid, bed$X4)] = 
        gLDSC_res[[pid]]$Taus[anno]
      s
    } + gLDSC_res[[pid]]$Taus[1]
}

N1 = 2e5
N2_seq = c(2e4,2e5)

n_causal = data.frame("n" = c(5,3,1),
                      "n1" = c(0,1,2),
                      "n2" = c(0,1,2))
rho = 0.8
sigma = 0.01 * sqrt(2)

for(s in 1:3)
{
  dir.create(sprintf("setting_%s", s), recursive = T)
}

foreach(block = select_block, .combine = "c") %dopar%
{
  data = variant[variant$block == block,]
  
  R = readRDS(sprintf("LD/%d.RData", block))
  
  data %<>% filter(is.element(rsid, colnames(R[[1]])))
  M = nrow(data)
  
  probs = cbind(data$var.EUR, data$var.EAS)
  
  foreach(s = 1:3, .combine = "c") %do%
  {
    set.seed(10*block + s)
    
    data$causal_1 = F
    data$causal_2 = F
    data$beta_1 = 0
    data$beta_2 = 0
    
    # Sample causal variants
    causal = causal_1 = causal_2 = numeric()
    candidate = candidate_1 = candidate_2 = 1:M
    
    if(n_causal$n[s] > 0)
    {
      causal = sample(candidate, n_causal$n[s], 
                      prob = rowSums(probs[candidate,]))
      
      candidate_1 %<>% setdiff(causal)
      candidate_2 %<>% setdiff(causal)
      
      data$causal_1[causal] = T
      data$causal_2[causal] = T
      data$beta_1[causal] = rnorm(n_causal$n[s], 0, sigma)
      data$beta_2[causal] = rnorm(n_causal$n[s], 
                                  rho*data$beta_1[causal], 
                                  sigma*sqrt(1 - rho^2))
    }
    
    if(n_causal$n1[s] > 0)
    {
      causal_1 = sample(candidate_1, n_causal$n1[s], 
                        prob = probs[candidate_1,1])
      
      candidate_2 %<>% setdiff(causal_1)
      
      data$causal_1[causal_1] = T
      data$beta_1[causal_1] = rnorm(n_causal$n1[s], 0, sigma)
    }
    
    if(n_causal$n2[s] > 0)
    {
      causal_2 = sample(candidate_2, n_causal$n2[s], 
                        prob = probs[candidate_2,2])
      
      data$causal_2[causal_2] = T
      data$beta_2[causal_2] = rnorm(n_causal$n2[s], 0, sigma)
    }
    
    LD = as.matrix(R[[1]])
    diag(LD) = 1
    data$z_1 = sqrt(N1) * (LD %*% data$beta_1)[,1] + 
      rmvnorm(1, sigma = LD)[1,]
    data$z_1 = sqrt(N1 - 1) * data$z_1 / sqrt(N1 - data$z_1^2)
    
    LD = as.matrix(R[[2]])
    diag(LD) = 1
    noise = rmvnorm(2, sigma = LD)
    for(m in 1:2)
    {
      z = sqrt(N2_seq[m]) * (LD %*% data$beta_2)[,1] + noise[m,]
      z = sqrt(N2_seq[m] - 1) * z / sqrt(N2_seq[m] - z^2)
      data[[sprintf("z_2.%d", N2_seq[m])]] = z
    }
    
    saveRDS(data, file = sprintf("setting_%s/%d.RData", s, block))
  }
}
