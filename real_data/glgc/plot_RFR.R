options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(ggh4x)
library(latex2exp)
library(scales)
library(Hmisc)

setwd("real_data/glgc")

dbs = c("UKBB", "GLGC")
phenos = c("HDL", "LDL", "TG", "TC")
pids = c("EUR", "AFR", "EAS")

methods = c("MACHINE + g-LDSC" = "MACHINE-gLDSC",
            "MACHINE" = "MACHINE",
            "MESuSiE" = "MESuSiE",
            "h2-D2 + g-LDSC" = "h2D2-gLDSC",
            "h2-D2" = "h2D2",
            "SuSiE" = "SuSiE",
            "RSparsePro" = "RSparsePro",
            "CARMA" = "CARMA")

data = foreach(pheno = phenos, .combine = "rbind") %do%
{
  foreach(method = names(methods), .combine = "rbind") %do%
  {
    if(method %in% c("MACHINE + g-LDSC", "MACHINE", "MESuSiE"))
    {
      CS_common = read.delim(sprintf(
        "results/%s/%s/joint/CS_common_num.txt", methods[method], pheno))
      variant_0.9_UKBB = read.delim(sprintf(
        "results/%s/%s/joint/variant_0.9_UKBB.txt", methods[method], pheno),
        check.names = F)
      
      foreach(pid = pids, .combine = "rbind") %do%
      {
        data.frame(
          "pheno" = pheno, "pid" = pid, "method" = method,
          "num_variants_0.9-UKBB" = sum(variant_0.9_UKBB[,sprintf(
            "num_variants_0.9.UKBB.%s", pid)]),
          "num_variants_0.9-UKBB_0.1.GLGC" = sum(variant_0.9_UKBB[,sprintf(
            "num_variants_0.9.UKBB_0.1.GLGC.%s", pid)]),
          "RFR_0.9" = 0, "RFR_0.9_sd" = 0,
          "UKBB_nCS" = sum(CS_common[,sprintf("UKBB_nCS.%s",pid)]),
          "UKBB_in_GLGC" = sum(CS_common[,sprintf("UKBB_in_GLGC.%s",pid)]),
          "RFR_CS95" = 0, "RFR_CS95_sd" = 0)
      }
    } else {
      foreach(pid = pids, .combine = "rbind") %do%
      {
        CS_common = read.delim(sprintf("results/%s/%s/joint/%s/CS_common_num.txt",
                                       methods[method], pheno, pid))
        variant_0.9_UKBB = read.delim(sprintf(
          "results/%s/%s/joint/%s/variant_0.9_UKBB.txt", methods[method], pheno, 
          pid), check.names = F)
        
        data.frame(
          "pheno" = pheno, "pid" = pid, "method" = method,
          "num_variants_0.9.UKBB" = sum(variant_0.9_UKBB$num_variants_0.9.UKBB),
          "num_variants_0.9.UKBB_0.1.GLGC" = sum(
            variant_0.9_UKBB$num_variants_0.9.UKBB_0.1.GLGC),
          "RFR_0.9" = 0, "RFR_0.9_sd" = 0,
          "UKBB_nCS" = sum(CS_common$UKBB_nCS),
          "UKBB_in_GLGC" = sum(CS_common$UKBB_in_GLGC),
          "RFR_CS95" = 0, "RFR_CS95_sd" = 0, check.names = F)
      }
    }
  }
}

data %<>% mutate(
  RFR_0.9 = num_variants_0.9.UKBB_0.1.GLGC / num_variants_0.9.UKBB,
  RFR_0.9_sd = sqrt(RFR_0.9 * (1-RFR_0.9) / num_variants_0.9.UKBB),
  RFR_CS95 = 1 - UKBB_in_GLGC / UKBB_nCS,
  RFR_CS95_sd = sqrt(RFR_CS95 * (1-RFR_CS95) / UKBB_nCS))
data$pheno %<>% factor(levels = phenos)
data$pid %<>% factor(levels = pids)
data$method %<>% factor(levels = names(methods))
data %<>% arrange(pheno, pid, method)

saveRDS(data, "RFR.RData")
