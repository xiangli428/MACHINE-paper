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
  foreach(method = names(methods), .combine = "rbind") %dopar%
  {
    output = data.frame("pheno" = pheno, "method" = method, "EUR_AFR_EAS" = 0, 
                        "EUR_AFR" = 0, "EUR_EAS" = 0, "AFR_EAS" = 0)
    
    if(file.exists(sprintf("results/%s/%s/%s/CS_overlap_union.txt", 
                           methods[method], pheno, "UKBB")) & 
       file.exists(sprintf("results/%s/%s/%s/CS_overlap_union.txt", 
                           methods[method], pheno, "GLGC")))
    {
      CS_overlap_union = foreach(db = dbs) %do%
      {
        tmp = read.delim(sprintf("results/%s/%s/%s/CS_overlap_union.txt", 
                                 methods[method], pheno, db))
        tmp[,5:7][is.na(tmp[,5:7])] = ""
        tmp
      }
      
      chr_block = merge(unique(CS_overlap_union[[1]][,3:4]),
                        unique(CS_overlap_union[[2]][,3:4]), sort = F)
      
      if(nrow(chr_block) > 0)
      {
        for(l in 1:nrow(chr_block))
        {
          i = chr_block[l,1]
          block = chr_block[l,2]
          
          sub = foreach(k = 1:2) %do%
          {
            filter(CS_overlap_union[[k]], chromosome == i & block == chr_block[l,2])
          }
          
          if(method %in% c("MACHINE + g-LDSC", "MACHINE", "MESuSiE"))
          {
            CS_info = read.delim(sprintf("%s/%s/%s/chr%d/%d/CS_info.txt",
                                         methods[method], pheno, "GLGC", i, block))
          } else {
            CS_info = foreach(pid = pids, .combine = "rbind") %do%
            {
              if(file.exists(sprintf("%s/%s/%s/%s/chr%d/%d/CS_info.txt",
                                     methods[method], pheno, "GLGC", pid, i, block)))
                read.delim(sprintf("%s/%s/%s/%s/chr%d/%d/CS_info.txt",
                                   methods[method], pheno, "GLGC", pid, i, block))
            }
          }
          CS_info$CS.tags[is.na(CS_info$CS.tags)] = ""
          
          for(j in 1:nrow(sub[[1]]))
          {
            Ks = which(sub[[1]][j,5:7] != '')
            Ms = which(apply(sub[[2]][,4+Ks] != '', 1, all))
            
            SNPs = foreach(m = Ms, .combine = "union") %do%
            {
              foreach(k = which(sub[[2]][m,5:7] != ''), .combine = "union") %do%
              {
                foreach(id = strsplit(sub[[2]][m, 4+k], ',')[[1]],
                        .combine = "union") %do%
                {
                  s = which(CS_info$CS == id)
                  c(strsplit(CS_info$CS.SNP[s], ',')[[1]],
                    strsplit(CS_info$CS.tags[s], ',')[[1]])
                }
              }
            }
            
            if(any(is.element(strsplit(sub[[1]]$union.SNP[j], ',')[[1]], SNPs)))
            {
              output[,paste(pids[Ks], collapse = '_')] %<>% add(1)
            }
          }
        }
      }
    }
    
    output
  }
}

data$pheno %<>% factor(levels = phenos)
data$method %<>% factor(levels = names(methods))

saveRDS(data, "CS_overlap_common_num.RData")

CS_overlap_num = readRDS("CS_overlap_num.RData")

data %<>% gather(pops, nCS_common, -pheno, -method)
CS_overlap_num = gather(CS_overlap_num[1:32,2:7], pops, nCS, -pheno, -method)
data$nCS = CS_overlap_num$nCS
data = data[,c(1:3,5,4)]
data %<>% mutate(RFR = 1 - nCS_common / nCS, RFR_sd = sqrt(RFR * (1-RFR) / nCS))

data$pops %<>% gsub("_", " & ", .)
data$pops %<>% factor(levels = c("EUR & AFR & EAS", "EUR & AFR", "EUR & EAS",
                                 "AFR & EAS"))

saveRDS(data, "CS_overlap_RFR.RData")

data = readRDS("CS_overlap_RFR.RData")

custom_theme = function()
{
  theme(
    aspect.ratio = 1,
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 8),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.25),
    axis.line = element_blank(),
    strip.text = element_text(size = 9),
    strip.background = element_rect(
      fill = "lightgray", color = "lightgray"),
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.spacing = unit(0.05, "in"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA))
}

p = ggplot(data, aes(x = method, y = RFR)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8,
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = RFR - RFR_sd, ymax = RFR + RFR_sd),
                width = 0.8, position = position_dodge2(),
                linewidth = 0.25, color = "black") +
  facet_grid(pops ~ pheno) +
  theme_classic() + custom_theme() +
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = NULL, y = "RFR (95% CS)", fill = "")

ggsave("CS_overlap_RFR.pdf", p,
       device = "pdf", width = 6.3, height = 6.8, units = "in", bg = "white")
