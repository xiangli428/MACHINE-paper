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
library(patchwork)
library(ggh4x)
library(latex2exp)
library(scales)

setwd("real_data/glgc")

dbs = c("GLGC","UKBB")
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

coding = c("NSF","NSM","NSN","ASS","DSS")

registerDoParallel(22)

data = foreach(pheno = phenos, .combine = "rbind") %do%
{
  eQTLs = foreach(i = 1:22, .combine = "rbind") %dopar%
  {
    unique(filter(read_delim(sprintf(
      "GTEx_eQTL_SuSiE_V10/%s/chr%d.txt.gz", pheno, i)), is.element(tissue, c(
        "Adipose_Subcutaneous","Adipose_Visceral_Omentum","Liver",
        "Whole_Blood")))[,c(4,10)])
  }
  eQTLs_keep = eQTLs$rsid[eQTLs$keep]
  eQTLs = eQTLs$rsid
  
  coding_SNPs = foreach(anno = coding, .combine = "rbind") %do%
  {
    foreach(i = 1:22, .combine = "rbind") %dopar%
    {
      read_delim(sprintf(
        "dbSNP151/%s/%s/chr%d.txt.gz", pheno, anno, i))[,c(4,10)]
    }
  }
  coding_SNPs %<>% unique()
  coding_SNPs_keep = coding_SNPs$rsid[coding_SNPs$keep]
  coding_SNPs = coding_SNPs$rsid
  
  min_p = read.delim(sprintf("data/%s/%s/min_p_sub.txt", "UKBB", pheno))
  
  foreach(db = dbs, .combine = "rbind") %do%
  {
    foreach(method = names(methods), .combine = "rbind") %dopar%
    {
      if(method %in% c("MACHINE + g-LDSC", "MACHINE", "MESuSiE"))
      {
        variants_0.9 = read_delim(sprintf(
          "results/%s/%s/%s/variants_0.9.txt.gz", methods[method], pheno, db))
        CS_info = read.delim(sprintf("results/%s/%s/%s/CS_info.txt", 
                                     methods[method], pheno, db))
      } else {
        variants_0.9 = foreach(pid = pids, .combine = "rbind") %do%
        {
          if(file.exists(sprintf(
            "results/%s/%s/%s/%s/variants_0.9.txt.gz", methods[method], pheno,
            db, pid)))
          read_delim(sprintf(
            "results/%s/%s/%s/%s/variants_0.9.txt.gz", methods[method], pheno,
            db, pid))
        }
        CS_info = foreach(pid = pids, .combine = "rbind") %do%
        {
          read.delim(sprintf("results/%s/%s/%s/%s/CS_info.txt", methods[method], 
                             pheno, db, pid))
        }
      }
      CS_info$CS.tags[is.na(CS_info$CS.tags)] = ""
      
      size = CS_info$CS.size + sapply(strsplit(CS_info$CS.tags, ','), length)
      
      output = data.frame(
        "pheno" = pheno, "db" = db, "method" = method,
        "num_variants" = sum(min_p$num_variants_keep),
        "num_variants_0.9" = length(unique(variants_0.9$rsid)),
        "num_CS95" = nrow(CS_info),
        "CS95_size_mean" = mean(size),
        "CS95_size_sd" = sd(size) / sqrt(nrow(CS_info)),
        "num_variants_coding" = length(coding_SNPs_keep),
        "num_variants_0.9_coding" = length(intersect(
          unique(variants_0.9$rsid), coding_SNPs_keep)),
        "num_CS95_coding" = 0,
        "num_variants_eQTL" = length(eQTLs_keep),
        "num_variants_0.9_eQTL" = length(intersect(
          unique(variants_0.9$rsid), eQTLs_keep)),
        "num_CS95_eQTL" = 0)
      if(method == "MESuSiE")
      {
        output$num_variants = sum(min_p$num_variants_keep.intersect)
      }
      
      for(l in 1:nrow(CS_info))
      {
        SNPs = c(strsplit(CS_info$CS.SNP[l], ",")[[1]], 
                 strsplit(CS_info$CS.tags[l], ",")[[1]])
        
        if(any(is.element(SNPs, coding_SNPs)))
        {
          output$num_CS95_coding %<>% add(1)
        }
        
        if(any(is.element(SNPs, eQTLs)))
        {
          output$num_CS95_eQTL %<>% add(1)
        }
      }
      
      output
    }
  }
}

data$pheno %<>% factor(levels = phenos)
data$db %<>% factor(levels = dbs)
data$method %<>% factor(levels = names(methods))

saveRDS(data, "eQTL_enrichment.RData")

data_RFR = readRDS("RFR.RData")

custom_theme = function()
{
  theme(
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
    plot.title = element_text(size = 10, hjust = 0.5),
    plot.margin = margin(l = 10, r = 5, b = 5, t = 5),
    panel.spacing = unit(0.05, "in"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA))
}


p1 = ggplot(data, aes(x = method, y = num_CS95)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8,
           position = position_dodge2()) +
  facet_grid(db ~ pheno, scales = "free_y") +
  theme_classic() + custom_theme() + 
  scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
  guides(fill = "none") +
  labs(x = NULL, y = "Number of 95% CSs", title = NULL)

p2 = ggplot(data, aes(x = method, y = CS95_size_mean)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8,
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = CS95_size_mean - CS95_size_sd, 
                    ymax = CS95_size_mean + CS95_size_sd),
                width = 0.8, position = position_dodge2(),
                linewidth = 0.25, color = "black") +
  facet_grid(db ~ pheno, scales = "free_y") +
  theme_classic() + custom_theme() + 
  scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
  guides(fill = "none") +
  labs(x = NULL, y = "Mean size of 95% CSs", title = NULL)

p3 = ggplot(data_RFR, aes(x = method, y = RFR_0.9)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8,
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = RFR_0.9 - RFR_0.9_sd, ymax = RFR_0.9 + RFR_0.9_sd),
                width = 0.8, position = position_dodge2(),
                linewidth = 0.25, color = "black") +
  facet_grid(pid ~ pheno) +
  theme_classic() + custom_theme() + 
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
  guides(fill = "none") +
  labs(x = NULL, y = "RFR", fill = "", title = TeX("CL or PIP $\\geq 0.9$"))

p4 = ggplot(data_RFR, aes(x = method, y = RFR_CS95)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8,
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = RFR_CS95 - RFR_CS95_sd, ymax = RFR_CS95 + RFR_CS95_sd),
                width = 0.8, position = position_dodge2(),
                linewidth = 0.25, color = "black") +
  facet_grid(pid ~ pheno) +
  theme_classic() + custom_theme() + 
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
  guides(fill = "none") +
  labs(x = NULL, y = "RFR", fill = "", title = "95% CS")

p5 = ggplot(data, aes(x = method, y = num_variants_0.9_coding / num_variants_0.9)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8,
           position = position_dodge2()) +
  facet_grid(db ~ pheno, scales = "free_y") +
  theme_classic() + custom_theme() + 
  scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Proportion", title = TeX("CL or PIP $\\geq 0.9$"))

p6 = ggplot(data, aes(x = method, y = num_CS95_coding / num_CS95)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8,
           position = position_dodge2()) +
  facet_grid(db ~ pheno, scales = "free_y") +
  theme_classic() + custom_theme() + 
  scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Proportion", title = "95% CS")

p7 = ggplot(data, aes(x = method, y = num_variants_0.9_eQTL / num_variants_0.9)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8,
           position = position_dodge2()) +
  facet_grid(db ~ pheno, scales = "free_y") +
  theme_classic() + custom_theme() + 
  scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Proportion", title = TeX("CL or PIP $\\geq 0.9$"))

p8 = ggplot(data, aes(x = method, y = num_CS95_eQTL / num_CS95)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8,
           position = position_dodge2()) +
  facet_grid(db ~ pheno, scales = "free_y") +
  theme_classic() + custom_theme() + 
  scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Proportion", title = "95% CS")

ggsave("CS_num_RFR_eQTL_enrichment.pdf",
       p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + 
         plot_layout(nrow = 4, guides = "collect", heights = c(1,1.5,1,1)) & 
         plot_annotation(tag_levels = "a") &
         theme(legend.position = "bottom",
               plot.margin = margin(l = 2, r = 2, b = 2, t = 2),
               plot.tag = element_text(size = 12, face = "bold")),
       device = "pdf", width = 9, height = 11.5, units = "in")
