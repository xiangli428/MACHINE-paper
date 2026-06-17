options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)
library(mvtnorm)

setwd("simulation/data")

pops = c("EUR" = 1, "AFR" = 4, "EAS" = 5)

registerDoParallel(10)

# Common variants
## UKBB
i = 1
Ks = 281
variant = foreach(block = 1:Ks, .combine = "rbind") %dopar%
{
  variant_pop = foreach(pid = names(pops)) %do%
  {
    read.delim(gzfile(sprintf(
      "~/Documents/GWAS/data/UKBB/%s/%s/LD/chr%s/%s/variant_reblock.txt.gz",
      pops[pid], "imputed_genotype_unique_info-0.8_maf-0.001_hwe-1e-6", i, block)))
  }

  merge(variant_pop[[1]], variant_pop[[2]], sort = F) %>% 
    merge(variant_pop[[3]], sort = F)
}
variant = variant[,-c(1,5)]

## 1kg
kg_dir = "~/Documents/GWAS/data/1000GENOME/phase3"
kg_variant = read_delim(sprintf("%s/variants_info/chr1.txt.gz", kg_dir),
                        delim = '\t')

for(pid in names(pops))
{
  kg_freq = read_table(sprintf("%s/%s/chr1.frq.gz", kg_dir, pid))
  if(pid == "EUR") kg_variant$SNP = kg_freq$SNP
  kg_variant[,sprintf("MAF_%s", pid)] = kg_freq$MAF
  kg_variant[,sprintf("REFisA2_%s", pid)] = kg_variant$REF == kg_freq$A2
}
rm(kg_freq)

kg_variant %<>% filter(MAF_EUR != 0 & MAF_AFR != 0 & MAF_EAS != 0)


variant_merge = merge(variant, kg_variant, by.x = c(2:5), by.y = c(1:4),
                      sort = F)

nSNPs = table(variant_merge$block)
select_block = which(nSNPs >= 900 & nSNPs <= 5100)
saveRDS(select_block, file = "select_block.RData")

variant_merge %<>% filter(is.element(block, select_block))
rm(variant, kg_variant)

# LD
variant = foreach(blk = select_block, .combine = "rbind") %dopar%
{
  dir.create(sprintf("LD/%d", blk), recursive = T)
  
  data = filter(variant_merge, block == blk)
  
  R = foreach(pid = names(pops)) %do%
  {
    LD = readMM(gzfile(sprintf(
      "~/Documents/GWAS/data/UKBB/%s/%s/LD/chr%s/%s/LD.mtx.gz",
      pops[pid], "imputed_genotype_unique_info-0.8_maf-0.001_hwe-1e-6", i, blk)))
    rownames(LD) = colnames(LD) = read.delim(sprintf(
      "~/Documents/GWAS/data/UKBB/%s/%s/LD/chr%s/%s/variant.txt",
      pops[pid], "imputed_genotype_unique_info-0.8_maf-0.001_hwe-1e-6", i, blk),
      header = F)[,1]
    as.matrix(LD[data$rsid, data$rsid])
  }

  for(k in 1:3)
  {
    R[[k]] = as(R[[k]][idx,idx], "dsCMatrix")
  }

  saveRDS(R, file = sprintf("LD/%d/LD_UKBB.RData", blk))
  
  data[idx,]
}

saveRDS(variant, "variant.RData")

select_block_info = foreach(blk = select_block, .combine = "rbind") %dopar%
{
  data = filter(variant, block == blk)
  write.table(data$SNP, sprintf("LD/%d/SNP_1kg.txt", blk), sep = '\t',
              row.names = F, col.names = F, quote = F)

  data.frame("block" = blk,
             "start" = min(data$position),
             "end" = max(data$position),
             "num_variants" = nrow(data))
}

write_delim(select_block_info, "select_block_info.txt", delim = '\t')
