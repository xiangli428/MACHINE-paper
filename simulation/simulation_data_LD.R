options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)
library(mvtnorm)

setwd("simulation/data")

# Common variants
pops = c("EUR" = 1, "EAS" = 5)
i = 1
Ks = 281

registerDoParallel(10)

variant = foreach(block = 1:Ks, .combine = "rbind") %dopar%
{
  variant_pop = foreach(pid = names(pops)) %do%
  {
    read.delim(gzfile(sprintf(
      "~/Documents/GWAS/data/UKBB/%s/%s/LD/chr%s/%s/variant_reblock.txt.gz",
      pops[pid], "imputed_genotype_unique_info-0.8_maf-0.001_hwe-1e-6", i, block)))
  }

  merge(variant_pop[[1]], variant_pop[[2]], sort = F)
}
variant = variant[,-c(1,5)]

nSNPs = table(variant$block)
set.seed(428)
select_block = sort(sample(which(nSNPs >= 1000 & nSNPs <= 6000), 200))
select_block = which(nSNPs >= 1000 & nSNPs <= 5390)
saveRDS(select_block, file = "select_block.RData")

variant %<>% filter(is.element(block, select_block))

# LD
dir.create("LD", recursive = T)

variant = foreach(block = select_block, .combine = "rbind") %dopar%
{
  data = variant[variant$block == block,]
  
  R = foreach(pid = names(pops)) %do%
  {
    LD = readMM(gzfile(sprintf(
      "~/Documents/GWAS/data/UKBB/%s/%s/LD/chr%s/%s/LD.mtx.gz",
      pops[pid], "imputed_genotype_unique_info-0.8_maf-0.001_hwe-1e-6", i, block)))
    rownames(LD) = colnames(LD) = read.delim(sprintf(
      "~/Documents/GWAS/data/UKBB/%s/%s/LD/chr%s/%s/variant.txt",
      pops[pid], "imputed_genotype_unique_info-0.8_maf-0.001_hwe-1e-6", i, block),
      header = F)[,1]
    as.matrix(LD[data$rsid, data$rsid])
  }

  idx = c(1)
  for(j in 2:nrow(data))
  {
    if(all(abs(R[[1]][idx,j]) < 0.99) | all(abs(R[[2]][idx,j]) < 0.99))
    {
      idx %<>% append(j)
    }
  }

  for(k in 1:2)
  {
    R[[k]] = as(R[[k]][idx,idx], "dsCMatrix")
  }

  saveRDS(R, file = sprintf("LD/%d.RData", block))
  
  data[idx,]
}

select_block_info = foreach(block = select_block, .combine = "rbind") %dopar%
{
  data = variant[variant$block == block,]

  data.frame("block" = block,
             "start" = min(variant$position[which(variant$block == block)]),
             "end" = max(variant$position[which(variant$block == block)]),
             "num_variants" = nrow(data))
}

write_delim(select_block_info, "select_block_info.txt", delim = '\t')