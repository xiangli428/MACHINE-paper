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
pops = c("EUR" = 1, "AFR" = 4, "EAS" = 5)

registerDoParallel(10)

select_block = readRDS("select_block.RData")

variant = readRDS("variant.RData")

gLDSC_res = foreach(pid = names(pops)) %do%
{
  readRDS(sprintf("%s/gLDSC_results/HDL/UKBB/%s/gLDSC_res.RData", 
                  "../../real_data/glgc", pid))
}
names(gLDSC_res) = names(pops)

for(pid in names(pops))
{
  Amatrix = read_delim(sprintf(
    "%s/baseline_bed_intersect/maf_g_input_53_%s/Amatrix.%s.annot", 
    "../../annotations", pops[pid], 1), delim = '\t', show_col_types = F)
  Amatrix = Amatrix[match(variant$rsid, Amatrix$SNP),]
  
  variant[[sprintf("var.%s", pid)]] = (as.matrix(Amatrix[,3:56]) %*% 
                                         gLDSC_res[[pid]]$Taus)[,1] %>% 
    pmax(max(.) / 20)
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

foreach(blk = select_block, .combine = "c") %dopar%
{
  data = filter(variant, block == blk)
  R = readRDS(sprintf("LD/%d/LD_UKBB.RData", blk))
  M = nrow(data)
  
  probs = cbind(data$var.EUR, data$var.AFR + data$var.EAS)
  
  foreach(s = 1:3, .combine = "c") %do%
  {
    set.seed(10*blk + s)
    
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
    
    LD = as.matrix(R[[2]])
    diag(LD) = 1
    noise = rmvnorm(2, sigma = LD)
    for(m in 1:2)
    {
      data[[sprintf("z_2.%d", N2_seq[m])]] = 
        sqrt(N2_seq[m]) * (LD %*% data$beta_2)[,1] + noise[m,]
    }
    
    LD = as.matrix(R[[3]])
    diag(LD) = 1
    noise = rmvnorm(2, sigma = LD)
    for(m in 1:2)
    {
      data[[sprintf("z_3.%d", N2_seq[m])]] = 
        sqrt(N2_seq[m]) * (LD %*% data$beta_2)[,1] + noise[m,]
    }
    
    saveRDS(data, file = sprintf("setting_%s/%d.RData", s, blk))
    NULL
  }
}
