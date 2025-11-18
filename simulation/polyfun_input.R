options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)
library(arrow)

setwd("simulation")

dir.create("polyfun_input/EUR", recursive = T)
dir.create("polyfun_input/EAS", recursive = T)
dir.create("polyfun_input/EUR_1kg", recursive = T)
dir.create("polyfun_input/EAS_1kg", recursive = T)

select_block = readRDS("data/select_block.RData")

pops = c("EUR" = 1, "EAS" = 5)

N1 = 2e5
N2_seq = c(2e4,2e5)

suffix = c("N1-200000", sprintf("N2-%d", N2_seq))
idx = c(1,2,2)

registerDoParallel(10)

# All SNPs
data = foreach(block = select_block, .combine = "rbind") %dopar%
{
  readRDS(sprintf("data/setting_%s/%d.RData", 1, block))
}

# Annotation files
Amatrix = read_delim("../annotations/baseline_bed_intersect/g_input_52_1/Amatrix.1.annot",
                     delim = '\t')
Amatrix = Amatrix[match(data$rsid, Amatrix$SNP),]
anno_name = names(Amatrix)[3:55]

anno = data.frame("SNP" = data$rsid,
                  "CHR" = ceiling(data$block / 13),
                  "BP" = data$position,
                  "A1" = data$first_allele,
                  "A2" = data$alternative_alleles) %>% cbind(Amatrix[,3:55])

foreach(i = 1:22, .combine = "c") %dopar%
{
  for(k in 1:2)
  {
    write_parquet(anno[anno$CHR == i,], sprintf(
      "polyfun_input/%s/annotations.%d.annot.parquet", names(pops)[k], i))
    write.table(t(colSums(
      anno[anno$CHR == i, 6:ncol(anno)])),
      sprintf("polyfun_input/%s/annotations.%d.l2.M", names(pops)[k], i),
      sep = ' ', row.names = F, col.names = F, quote = F)
    
    write_parquet(anno[anno$CHR == i,], sprintf(
      "polyfun_input/%s_1kg/annotations.%d.annot.parquet", names(pops)[k], i))
    write.table(t(colSums(
      anno[anno$CHR == i, 6:ncol(anno)])),
      sprintf("polyfun_input/%s_1kg/annotations.%d.l2.M", names(pops)[k], i),
      sep = ' ', row.names = F, col.names = F, quote = F)
  }
  
  NULL
}

# LD scores
for(k in 1:2)
{
  ldscore = foreach(block = select_block, .combine = "rbind") %dopar%
  {
    data_sub = data[data$block == block,]
    anno_sub = anno[data$block == block,]
    
    R = readRDS(sprintf("data/LD/%d.RData", block))[[k]] %>% as.matrix()
    diag(R) = 1
    
    ldscore_sub = anno_sub[,1:5]
    cbind(
      ldscore_sub,
      foreach(an = anno_name, .combine = "cbind") %do%
      {
        colSums(R^2 * anno_sub[,an])
      }
    )
  }
  names(ldscore)[6:ncol(ldscore)] = anno_name
  
  w = ldscore[,1:6]
  names(w)[6] = "L2"
  
  foreach(i = 1:22, .combine = "c") %dopar%
  {
    write_parquet(filter(ldscore, CHR == i), sprintf(
      "polyfun_input/%s/annotations.%d.l2.ldscore.parquet", names(pops)[k], i))
    write_parquet(filter(w, CHR == i), sprintf(
      "polyfun_input/%s/weights.%d.l2.ldscore.parquet", names(pops)[k], i))
    NULL
  }
}

# LD scores 1kg
for(k in 1:2)
{
  ldscore = foreach(block = select_block, .combine = "rbind") %dopar%
  {
    data_sub = data[data$block == block,]
    anno_sub = anno[data$block == block,]
    
    R = readRDS(sprintf("data/LD_1kg/%d.RData", block))[[k]] %>% as.matrix()
    diag(R) = 1
    
    ldscore_sub = anno_sub[,1:5]
    cbind(
      ldscore_sub,
      foreach(an = anno_name, .combine = "cbind") %do%
      {
        colSums(R^2 * anno_sub[,an])
      }
    )
  }
  names(ldscore)[6:ncol(ldscore)] = anno_name
  
  w = ldscore[,1:6]
  names(w)[6] = "L2"
  
  foreach(i = 1:22, .combine = "c") %dopar%
    {
      write_parquet(filter(ldscore, CHR == i), sprintf(
        "polyfun_input/%s_1kg/annotations.%d.l2.ldscore.parquet", names(pops)[k], i))
      write_parquet(filter(w, CHR == i), sprintf(
        "polyfun_input/%s_1kg/weights.%d.l2.ldscore.parquet", names(pops)[k], i))
      NULL
    }
}

# sumstat

for(s in 1:3)
{
  dir.create(sprintf("polyfun_input/sumstat/setting_%d", s),
             recursive = T)
  
  data = foreach(block = select_block, .combine = "rbind") %dopar%
  {
    readRDS(sprintf("data/setting_%s/%d.RData", s, block))
  }
  
  for(k in 1:3)
  {
    sumstat = data.frame("SNP" = data$rsid,
                         "CHR" = ceiling(data$block / 13),
                         "BP" = data$position,
                         "A1" = data$first_allele,
                         "A2" = data$alternative_alleles,
                         "Z" = data[,12+k],
                         "N" = c(N1,N2_seq)[k])
    
    write_parquet(sumstat, sprintf(
      "polyfun_input/sumstat/setting_%d/sumstat.%s.parquet", s, suffix[k]))
  }
}
