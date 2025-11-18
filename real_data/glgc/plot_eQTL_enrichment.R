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

annos = list("Coding" = c("NSF","NSM","NSN"),
             "Regulatory" = c("R3","R5","SYN","U3","U5"))

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
  
  anno_SNPs = foreach(type = names(annos)) %do%
  {
    foreach(anno = annos[[type]], .combine = "rbind") %do%
    {
      foreach(i = 1:22, .combine = "rbind") %dopar%
      {
        read_delim(sprintf(
          "dbSNP151/%s/%s/chr%d.txt.gz", pheno, anno, i))[,c(4,10)]
      }
    }
  }
  names(anno_SNPs) = names(annos)
  anno_SNPs[[1]] %<>% unique()
  anno_SNPs[[2]] %<>% unique()
  anno_SNPs[[2]] %<>% filter(!is.element(rsid, anno_SNPs[[1]]$rsid))
  
  anno_SNPs_keep = list()
  for(type in names(annos))
  {
    anno_SNPs_keep[[type]] = anno_SNPs[[type]]$rsid[anno_SNPs[[type]]$keep]
    anno_SNPs[[type]] = anno_SNPs[[type]]$rsid
  }
  
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
      
      output = data.frame(
        "pheno" = pheno, "db" = db, "method" = method,
        "num_variants" = sum(min_p$num_variants_keep),
        "num_variants_0.9" = length(unique(variants_0.9$rsid)),
        "num_CS95" = nrow(CS_info),
        "num_variants_eQTL" = length(eQTLs_keep),
        "num_variants_0.9_eQTL" = length(intersect(
          unique(variants_0.9$rsid), eQTLs_keep)),
        "num_CS95_eQTL" = 0)
      if(method == "MESuSiE")
      {
        output$num_variants = sum(min_p$num_variants_keep.intersect)
      }
      for(type in names(annos))
      {
        output[[sprintf("num_variants_%s", type)]] = length(anno_SNPs_keep[[type]])
        output[[sprintf("num_variants_0.9_%s", type)]] = length(intersect(
          unique(variants_0.9$rsid), anno_SNPs_keep[[type]]))
        output[[sprintf("num_CS95_%s", type)]] = 0
      }
      
      for(l in 1:nrow(CS_info))
      {
        SNPs = c(strsplit(CS_info$CS.SNP[l], ",")[[1]], 
                 strsplit(CS_info$CS.tags[l], ",")[[1]])
        
        if(any(is.element(SNPs, eQTLs)))
        {
          output$num_CS95_eQTL %<>% add(1)
        }
        
        if(any(is.element(SNPs, anno_SNPs[["Coding"]])))
        {
          output$num_CS95_Coding %<>% add(1)
        } else if(any(is.element(SNPs, anno_SNPs[["Regulatory"]]))) {
          output$num_CS95_Regulatory %<>% add(1)
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

data = readRDS("eQTL_enrichment.RData")
data_RFR = readRDS("RFR.RData")

data_anno_0.9 = data[,1:3]
data_anno_0.9$Coding = data$num_variants_0.9_Coding / data$num_variants_0.9
data_anno_0.9[,"Putative regulatory"] = data$num_variants_0.9_Regulatory / 
  data$num_variants_0.9
data_anno_0.9 %<>% gather(type, prop, -pheno, -db, -method)
data_anno_0.9$type %<>% factor(levels = c("Putative regulatory", "Coding"))

data_anno_CS95 = data[,1:3]
data_anno_CS95$Coding = data$num_CS95_Coding / data$num_CS95
data_anno_CS95[,"Putative regulatory"] = data$num_CS95_Regulatory / data$num_CS95
data_anno_CS95 %<>% gather(type, prop, -pheno, -db, -method)
data_anno_CS95$type %<>% factor(levels = c("Putative regulatory", "Coding"))

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

p1 = ggplot(data_RFR, aes(x = method, y = RFR_0.9)) +
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

p2 = ggplot(data_RFR, aes(x = method, y = RFR_CS95)) +
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

p3 = ggplot(data_anno_0.9, aes(
  x = method, y = prop, color = NULL, fill = method, alpha = type)) +
  geom_bar(stat = "identity") +
  facet_grid(db ~ pheno, scales = "fixed") +
  theme_classic() + custom_theme() + 
  scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
  scale_alpha_manual(breaks = c("Coding","Putative regulatory"), 
                     values = c(1,0.6)) +
  guides(fill = "none", alpha = "none") +
  labs(x = NULL, y = "Proportion", title = TeX("CL or PIP $\\geq 0.9$"))

p4 = ggplot(data_anno_CS95, aes(
  x = method, y = prop, color = NULL, fill = method, alpha = type)) +
  geom_bar(stat = "identity") +
  facet_grid(db ~ pheno, scales = "fixed") +
  theme_classic() + custom_theme() + 
  scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
  scale_alpha_manual(breaks = c("Coding","Putative regulatory"), 
                     values = c(1,0.6)) +
  guides(fill = guide_legend(nrow = 2, order = 1), 
         alpha = guide_legend(nrow = 2, order = 2)) +
  labs(x = NULL, y = "Proportion", title = "95% CS")

p5 = ggplot(data, aes(x = method, y = num_variants_0.9_eQTL / num_variants_0.9)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8,
           position = position_dodge2()) +
  facet_grid(db ~ pheno, scales = "free_y") +
  theme_classic() + custom_theme() + 
  scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
  guides(fill = "none") +
  labs(x = NULL, y = "Proportion", title = TeX("CL or PIP $\\geq 0.9$"))

p6 = ggplot(data, aes(x = method, y = num_CS95_eQTL / num_CS95)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8,
           position = position_dodge2()) +
  facet_grid(db ~ pheno, scales = "free_y") +
  theme_classic() + custom_theme() + 
  scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
  guides(fill = "none") +
  labs(x = NULL, y = "Proportion", title = "95% CS")

p12 = ggarrange(p1, p2,
                nrow = 1, ncol = 2, common.legend = T, legend = "bottom",
                widths = c(1,1), labels = c("a","b"),
                font.label = list(size = 12, color = "black", face = "bold",
                                  family = NULL), align = "hv")

p3456 = ggarrange(p3, p4, p5, p6, 
                  nrow = 2, ncol = 2, common.legend = T, legend = "bottom",
                  widths = c(1,1), heights = c(1,1), 
                  labels = c("c","d","e","f"),
                  font.label = list(size = 12, color = "black", face = "bold",
                                    family = NULL), align = "hv")

ggsave("RFR_eQTL_enrichment.pdf",
       ggarrange(p12, p3456,
                 nrow = 2, ncol = 1, common.legend = T, legend = "bottom",
                 heights = c(1,1.66),
                 font.label = list(size = 12, color = "black", face = "bold",
                                   family = NULL), align = "hv"),
       device = "pdf", width = 9, height = 9.08, units = "in", bg = "white")
