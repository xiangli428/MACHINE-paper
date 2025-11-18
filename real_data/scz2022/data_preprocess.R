options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)

logit = function(x)
{
  return(log(x / (1-x)))
}

setwd("real_data/scz2022/data")

data_dir = "~/Documents/GWAS/data/Summarydata/scz2022"
pops = c("EUR" = 1, "EAS" = 5)

ncas = c("EUR" = 53386, "EAS" = 14004)
ncon = c("EUR" = 77258, "EAS" = 16757)

registerDoParallel(15)

for(pid in c("EUR", "EAS"))
{
  UKBB_dir = sprintf("%s/%d/imputed_genotype_unique_info-0.8_maf-0.001_hwe-1e-6",
                     "~/Documents/GWAS/data/UKBB", pops[pid])
  
  foreach(i = 1:22) %dopar%
  {
    dir.create(sprintf("%s/chr%s", pid, i), recursive = T)
    
    # GWAS summary data
    data = read_delim(gzfile(sprintf("%s/%s/chr%d.tsv.gz", data_dir, pid, i)),
                      delim = '\t')
    n_before = nrow(data)
    filter_sample_size = sum(data$NCAS != ncas[pid] | data$NCON != ncon[pid])
    data %<>% filter(NCAS == ncas[pid] & NCON == ncon[pid])
    
    filter_info = sum(data$IMPINFO < 0.8)
    data %<>% filter(IMPINFO >= 0.8)
    
    data %<>% mutate(FREQ = 1 - (FCAS * ncas[pid] + FCON * ncon[pid]) / 
                       sum(ncas[pid] + ncon[pid]))
    
    # SNP info
    variant = read_delim(gzfile(sprintf(
      "%s/LD/chr%d/variant_reblock.txt.gz", UKBB_dir, i)), delim = '\t')[,-5]
    afreq = read_delim(sprintf("%s/genotype_bgen-1.2/chr%d.afreq", UKBB_dir, i), 
                       delim = '\t')
    variant = merge(variant, afreq[,2:5],
                    by.x = c("rsid","first_allele","alternative_alleles"),
                    by.y = c("ID","REF","ALT"),
                    sort = F)
    rm(afreq)
    variant = variant[!is.na(variant$block),]
    
    variant_merge_1 = merge(
      variant, data, 
      by.x = c("chromosome", "position", "first_allele", "alternative_alleles"),
      by.y = c("CHROM", "POS", "A1", "A2"), sort = F)
    
    variant_merge_2 = merge(
      variant, data, 
      by.x = c("chromosome", "position", "first_allele", "alternative_alleles"),
      by.y = c("CHROM", "POS", "A2", "A1"), sort = F)
    
    if(nrow(variant_merge_2) > 0)
    {
      variant_merge_1$ALTisA2 = T
      variant_merge_2$ALTisA2 = F
      variant_merge_2$BETA = -variant_merge_2$BETA
      variant_merge = rbind(variant_merge_1, variant_merge_2)
    } else {
      variant_merge = variant_merge_1
      variant_merge$ALTisA2 = T
    }
    
    filter_merge = nrow(data) - nrow(variant_merge)
    
    variant_merge %<>% mutate(diff_logit_freq = logit(ALT_FREQS) - logit(FREQ))
    filter_freq = sum(abs(variant_merge$diff_logit_freq) > 0.5)
    variant_merge %<>% filter(abs(diff_logit_freq) <= 0.5)
    
    variant_merge %<>% arrange(position)
    variant_merge %<>% mutate(alternate_ids = sprintf(
      "%s:%s_%s_%s", chromosome, position, first_allele, alternative_alleles))
    
    for(block in unique(variant_merge$block))
    {
      write_delim(variant_merge[variant_merge$block == block,], 
                  gzfile(sprintf("%s/chr%d/block_%s.txt.gz", pid, i, block)), 
                  delim = '\t')
    }
    
    NULL
  }
}
