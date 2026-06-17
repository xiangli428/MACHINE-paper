options(stringsAsFactors = F, check.names = F)

library(readr)
library(tidyr)
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

setwd("real_data/scz2022")

pids = c("EUR", "EAS")

methods = c("MACHINE + g-LDSC" = "MACHINE-gLDSC",
            "MACHINE" = "MACHINE",
            "MESuSiE" = "MESuSiE",
            "h2-D2 + g-LDSC" = "h2D2-gLDSC",
            "h2-D2" = "h2D2",
            "SuSiE" = "SuSiE",
            "RSparsePro" = "RSparsePro",
            "CARMA" = "CARMA")

data = foreach(method = names(methods), .combine = "rbind") %dopar%
{
  CS_overlap = read.delim(sprintf("results/%s/CS_overlap.txt", methods[method]))
  
  if(method %in% c("MACHINE + g-LDSC", "MACHINE", "MESuSiE"))
  {
    CS_info = read.delim(sprintf("results/%s/CS_info.txt", methods[method]))
  } else {
    CS_info = foreach(pid = pids, .combine = "rbind") %do%
    {
      read.delim(sprintf("results/%s/%s/CS_info.txt", methods[method], pid))
    }
  }
  
  data.frame("method" = method,
             "EUR_EAS" = nrow(CS_overlap),
             "EUR-specific" = length(setdiff(grep("EUR", CS_info$CS, value = T),
                                             CS_overlap$CS.EUR)), 
             "EAS-specific" = length(setdiff(grep("EAS", CS_info$CS, value = T),
                                             CS_overlap$CS.EAS)),
             check.names = F)
}

data$method %<>% factor(levels = names(methods))

saveRDS(data, "CS_overlap_num.RData")
