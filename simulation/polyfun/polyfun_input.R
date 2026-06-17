options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)
library(arrow)

setwd("simulation")

select_block = readRDS("data/select_block.RData")

pops = c("EUR" = 1, "AFR" = 4, "EAS" = 5)
N1 = 2e5
N2_seq = c(2e4,2e5)

suffix = c("N1-200000", sprintf("N2-%d", N2_seq), sprintf("N3-%d", N2_seq))
idx = c(1,2,2,3,3)

registerDoParallel(20)

# All SNPs
data = readRDS("data/variant.RData")


for(k in 1:3)
{
  pid = names(pops)[k]
  # Annotation files
  Amatrix = read_delim(sprintf("%s/maf_g_input_53_%s/Amatrix.1.annot",
                               "../annotations/baseline_bed_intersect", pops[k]),
                       delim = '\t')
  Amatrix = Amatrix[match(data$rsid, Amatrix$SNP),3:56]
  anno_name = names(Amatrix)
  
  anno = data.frame("SNP" = data$rsid,
                    "CHR" = ceiling(data$block / 13),
                    "BP" = data$position,
                    "A1" = data$first_allele,
                    "A2" = data$alternative_alleles) %>% cbind(Amatrix)
  
  # LD scores
  for(t in c("UKBB","1kg"))
  {
    dir.create(sprintf("polyfun/input/annotations/%s/%s", t, pid), recursive = T)

    ldscore = foreach(block = select_block, .combine = "rbind") %dopar%
    {
      anno_sub = anno[data$block == block,]

      R = readRDS(sprintf("data/LD/%d/LD_%s.RData", block, t))[[k]] %>% as.matrix()
      diag(R) = 1

      cbind(
        anno_sub[,1:5],
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
      write_parquet(filter(anno, CHR == i), sprintf(
        "polyfun/input/annotations/%s/%s/annotations.%d.annot.parquet", t, pid, i))
      write.table(t(colSums(anno[anno$CHR == i, 6:ncol(anno)])),
                  sprintf("polyfun/input/annotations/%s/%s/annotations.%d.l2.M", t, pid, i),
                  sep = ' ', row.names = F, col.names = F, quote = F)
      write_parquet(filter(ldscore, CHR == i), sprintf(
        "polyfun/input/annotations/%s/%s/annotations.%d.l2.ldscore.parquet", t, pid, i))
      write_parquet(filter(w, CHR == i), sprintf(
        "polyfun/input/annotations/%s/%s/weights.%d.l2.ldscore.parquet", t, pid, i))
      NULL
    }
  }
}

# sumstat
for(s in 1:3)
{
  dir.create(sprintf("polyfun/input/sumstat/setting_%d", s), recursive = T)
  
  data = foreach(block = select_block, .combine = "rbind") %dopar%
  {
    readRDS(sprintf("data/setting_%s/%d.RData", s, block))
  }
  
  for(k in 1:5)
  {
    sumstat = data.frame("SNP" = data$rsid,
                         "CHR" = ceiling(data$block / 13),
                         "BP" = data$position,
                         "A1" = data$first_allele,
                         "A2" = data$alternative_alleles,
                         "Z" = data[,20+k],
                         "N" = c(N1,N2_seq,N2_seq)[k])
    
    write_parquet(sumstat, sprintf(
      "polyfun/input/sumstat/setting_%d/sumstat.%s.parquet", s, suffix[k]))
  }
}