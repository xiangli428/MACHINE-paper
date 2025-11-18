options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)

high_LD = read.delim("~/Documents/GWAS/data/UKBB/high_LD_block.txt")
high_LD_block = split(high_LD$block, high_LD$chromosome)

setwd("real_data/glgc/data")

data_dir = c(
  "UKBB" = "~/Documents/GWAS/data/Summarydata/Pan-UKBB",
  "GLGC" = "~/Documents/GWAS/data/Summarydata/glgc-lipids2021/ancestry_specific")

phenos = list("UKBB" = c(
  "HDL" = "30760-HDL_cholesterol", "LDL" = "30780-LDL_direct", 
  "TG" = "30870-Triglycerides", "TC" = "30690-Cholesterol"), 
  "GLGC" = c("HDL" = "HDL", "LDL" = "LDL", "TG" = "logTG", "TC" = "TC"))

pops = c("EUR" = 1, "AFR" = 4, "EAS" = 5)

registerDoParallel(22)

sample_size = foreach(pheno = names(phenos$GLGC), .combine = "rbind") %do%
{
  data_UKBB = read_delim(gzfile(sprintf(
    "%s/%s/biomarkers-%s-both_sexes-irnt.tsv.bgz", data_dir["UKBB"], 
    phenos[["UKBB"]][pheno], gsub("-.*$", "", phenos[["UKBB"]][pheno]))), 
    delim = '\t')
  
  foreach(pid = names(pops), .combine = "rbind") %do%
  {
    UKBB_dir = sprintf("%s/%d/imputed_genotype_unique_info-0.8_maf-0.001_hwe-1e-6",
                       "~/Documents/GWAS/data/UKBB", pops[pid])
    
    idx = which(!is.na(data_UKBB[,sprintf("beta_%s", pid)]))
    
    data = list("UKBB" = data_UKBB[idx,1:4])
    data$UKBB[,"af"] = data_UKBB[idx, sprintf("af_%s", pid)]
    data$UKBB[,"beta"] = data_UKBB[idx, sprintf("beta_%s", pid)]
    data$UKBB[,"se"] = data_UKBB[idx, sprintf("se_%s", pid)]
    data$UKBB[,"neglog10_pval"] = data_UKBB[idx, sprintf("neglog10_pval_%s", pid)]
    
    if(pid != "EAS")
    {
      data[["GLGC"]] = read_table(gzfile(sprintf(
        "%s/%s_INV_%s_HRC_1KGP3_others_ALL.meta.singlevar.results.gz", 
        data_dir[["GLGC"]], phenos[["GLGC"]][pheno], pid)))
    } else {
      data[["GLGC"]] = read_table(gzfile(sprintf(
        "%s/%s_INV_%s_1KGP3_ALL.meta.singlevar.results.gz", 
        data_dir[["GLGC"]], phenos[["GLGC"]][pheno], pid)))
    }
    
    n = ceiling(max(data[["GLGC"]]$N) * 0.9)
    
    data[["GLGC"]] %<>% filter(N >= n)
    
    # SNP info
    foreach(i = 1:22, .combine = "c") %dopar%
    {
      dir.create(sprintf("UKBB/%s/%s/chr%d", pheno, pid, i), recursive = T)
      dir.create(sprintf("GLGC/%s/%s/chr%d", pheno, pid, i), recursive = T)
      
      variant = read_delim(gzfile(sprintf(
        "%s/LD/chr%d/variant_reblock.txt.gz", UKBB_dir, i)), delim = '\t')[,-5]
      afreq = read_delim(sprintf("%s/genotype_bgen-1.2/chr%d.afreq", UKBB_dir, i), 
                         delim = '\t')
      variant = merge(variant, afreq[,2:5],
                      by.x = c("rsid","first_allele","alternative_alleles"),
                      by.y = c("ID","REF","ALT"),
                      sort = F)
      rm(afreq)
      variant %<>% filter(!is.na(block))
      
      data_merge = list()
      
      #UKBB
      sub = filter(data$UKBB, chr == i)
      
      data_merge_1 = merge(
        variant, sub, 
        by.x = c("chromosome", "position", "first_allele", "alternative_alleles"),
        by.y = c("chr", "pos", "ref", "alt"), sort = F)
      
      data_merge_2 = merge(
        variant, sub, 
        by.x = c("chromosome", "position", "first_allele", "alternative_alleles"),
        by.y = c("chr", "pos", "alt", "ref"), sort = F)
      
      if(nrow(data_merge_2) > 0)
      {
        data_merge_1$ALTisALT = T
        data_merge_2$ALTisALT = F
        data_merge_2 %<>% mutate(af = 1-af, beta = -beta)
        data_merge$UKBB = rbind(data_merge_1, data_merge_2)
      } else {
        data_merge$UKBB = data_merge_1
        data_merge$UKBB$ALTisALT = T
      }
      
      data_merge$UKBB %<>% arrange(position)
      data_merge$UKBB %<>% mutate(alternate_ids = sprintf(
        "%s:%s_%s_%s", chromosome, position, first_allele, alternative_alleles))
      
      rm(sub, data_merge_1, data_merge_2)
      
      #GLGC
      sub = filter(data$GLGC, CHROM == i)
      
      data_merge_1 = merge(
        variant, sub, 
        by.x = c("chromosome", "position", "first_allele", "alternative_alleles"),
        by.y = c("CHROM", "POS_b37", "REF", "ALT"), sort = F)
      
      data_merge_2 = merge(
        variant, sub, 
        by.x = c("chromosome", "position", "first_allele", "alternative_alleles"),
        by.y = c("CHROM", "POS_b37", "ALT", "REF"), sort = F)
      
      if(nrow(data_merge_2) > 0)
      {
        data_merge_1$ALTisALT = T
        data_merge_2$ALTisALT = F
        data_merge_2 %<>% mutate(POOLED_ALT_AF = 1-POOLED_ALT_AF,
                                 EFFECT_SIZE = -EFFECT_SIZE)
        data_merge$GLGC = rbind(data_merge_1, data_merge_2)
      } else {
        data_merge$GLGC = data_merge_1
        data_merge$GLGC$ALTisALT = T
      }
      
      data_merge$GLGC %<>% arrange(position)
      data_merge$GLGC %<>% mutate(alternate_ids = sprintf(
        "%s:%s_%s_%s", chromosome, position, first_allele, alternative_alleles),
        Z = EFFECT_SIZE / SE, Z = Z * sqrt(n / (N+Z^2)), Z = Z * sqrt(n / (N-Z^2)),
        neglog10_pval = -(pnorm(-abs(Z), log.p = T) + log(2)) / log(10))
      
      rm(sub, data_merge_1, data_merge_2)
      
      # select common SNPs
      SNPs = intersect(data_merge$UKBB$alternate_ids, data_merge$GLGC$alternate_ids)
      
      data_merge$UKBB %<>% filter(alternate_ids %in% SNPs)
      data_merge$GLGC %<>% filter(alternate_ids %in% SNPs)
      
      for(block in setdiff(unique(data_merge$UKBB$block), 
                           high_LD_block[[as.character(i)]]))
      {
        if(sum(data_merge$UKBB$block == block) >= 50)
        {
          write_delim(data_merge$UKBB[data_merge$UKBB$block == block,], gzfile(
            sprintf("UKBB/%s/%s/chr%d/block_%d.txt.gz", pheno, pid, i, block)),
            delim = '\t')
          
          write_delim(data_merge$GLGC[data_merge$GLGC$block == block,], gzfile(
            sprintf("GLGC/%s/%s/chr%d/block_%d.txt.gz", pheno, pid, i, block)),
            delim = '\t')
        }
      }
    }
    
    data.frame("pheno" = pheno, "pid" = pid, "n" = n)
  }
}

write_delim(sample_size, "GLGC/sample_size.txt", delim = '\t')


# sample size info
sample_size = foreach(db = dbs, .combine = "cbind") %dopar%
{
  read.delim(sprintf("%s/sample_size.txt", db))
}
sample_size = sample_size[,-c(4,5)]
names(sample_size)[3:4] = c("UKBB.n", "GLGC.n")

write_delim(sample_size, "sample_size.txt", delim = '\t')
