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

annos = list("Coding" = c("NSF","NSM","NSN"),
             "Regulatory" = c("R3","R5","SYN","U3","U5"))

registerDoParallel(22)

eQTLs = foreach(i = 1:22, .combine = "c") %dopar%
{
  unique(filter(read_delim(sprintf(
    "GTEx_eQTL_SuSiE_V10/chr%d.txt.gz", i)), startsWith(tissue, "Brain"))$rsid)
}

anno_SNPs = foreach(type = names(annos)) %do%
{
  foreach(anno = annos[[type]], .combine = "c") %do%
  {
    foreach(i = 1:22, .combine = "c") %dopar%
    {
      read_delim(sprintf(
        "dbSNP151/%s/chr%d.txt.gz", anno, i))$rsid
    }
  }
}

names(anno_SNPs) = names(annos)
anno_SNPs[[1]] %<>% unique()
anno_SNPs[[2]] %<>% unique()
anno_SNPs[[2]] %<>% setdiff(anno_SNPs[[1]])

min_p = read.delim("data/min_p_sub.txt")

data = foreach(method = names(methods), .combine = "rbind") %dopar%
{
  if(method %in% c("MACHINE + g-LDSC", "MACHINE", "MESuSiE"))
  {
    variants_0.5 = read.delim(sprintf(
      "results/%s/variants_0.5.txt.gz", methods[method]))
    
    CS_info = read.delim(sprintf("results/%s/CS_info.txt", methods[method]))
  } else {
    variants_0.5 = foreach(pid = pids, .combine = "rbind") %do%
    {
      if(file.exists(sprintf(
        "results/%s/%s/variants_0.5.txt.gz", methods[method], pid)))
        read.delim(sprintf(
          "results/%s/%s/variants_0.5.txt.gz", methods[method], pid))
    }
    
    CS_info = foreach(pid = pids, .combine = "rbind") %do%
    {
      read.delim(sprintf("results/%s/%s/CS_info.txt", methods[method], pid))
    }
  }
  
  output = data.frame(
    "method" = method,
    "num_variants" = sum(min_p$num_variants),
    "num_variants_0.5" = length(unique(variants_0.5$rsid)),
    "num_CS95" = nrow(CS_info),
    "num_variants_eQTL" = length(eQTLs),
    "num_variants_0.5_eQTL" = length(intersect(
      unique(variants_0.5$rsid), eQTLs)),
    "num_CS95_eQTL" = 0)
  for(type in names(annos))
  {
    output[[sprintf("num_variants_%s", type)]] = length(anno_SNPs[[type]])
    output[[sprintf("num_variants_0.5_%s", type)]] = length(intersect(
      unique(variants_0.5$rsid), anno_SNPs[[type]]))
    output[[sprintf("num_CS95_%s", type)]] = 0
  }
  
  for(l in 1:nrow(CS_info))
  {
    SNPs = strsplit(CS_info$CS.SNP[l], ",")[[1]]
    
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

data$method %<>% factor(levels = names(methods))

saveRDS(data, "eQTL_enrichment.RData")

data = readRDS("eQTL_enrichment.RData")

data_anno_0.5 = data["method"]
data_anno_0.5$Coding = data$num_variants_0.5_Coding / data$num_variants_0.5
data_anno_0.5[,"Putative regulatory"] = data$num_variants_0.5_Regulatory / 
  data$num_variants_0.5
data_anno_0.5 %<>% gather(type, prop, -method)
data_anno_0.5$type %<>% factor(levels = c("Putative regulatory", "Coding"))

data_anno_CS95 = data["method"]
data_anno_CS95$Coding = data$num_CS95_Coding / data$num_CS95
data_anno_CS95[,"Putative regulatory"] = data$num_CS95_Regulatory / data$num_CS95
data_anno_CS95 %<>% gather(type, prop, -method)
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

p1 = readRDS("plot_CS_overlap_num.RData")

p2 = ggplot(data_anno_0.5, aes(
  x = method, y = prop, color = NULL, fill = method, alpha = type)) +
  geom_bar(stat = "identity") +
  theme_classic() + custom_theme() + 
  scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
  scale_alpha_manual(breaks = c("Coding","Putative regulatory"), 
                     values = c(1,0.6)) +
  guides(fill = "none", alpha = "none") +
  labs(x = NULL, y = "Proportion", title = TeX("CL or PIP $\\geq 0.5$"))

p3 = ggplot(data_anno_CS95, aes(
  x = method, y = prop, color = NULL, fill = method, alpha = type)) +
  geom_bar(stat = "identity") +
  theme_classic() + custom_theme() + 
  scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
  scale_alpha_manual(breaks = c("Coding","Putative regulatory"), 
                     values = c(1,0.6)) +
  guides(fill = guide_legend(nrow = 2, order = 1), 
         alpha = guide_legend(nrow = 2, order = 2)) +
  labs(x = NULL, y = "Proportion", title = "95% CS")

p4 = ggplot(data, aes(x = method, y = num_variants_0.5_eQTL / num_variants_0.5)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8,
           position = position_dodge2()) +
  theme_classic() + custom_theme() + 
  scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
  guides(fill = "none") +
  labs(x = NULL, y = "Proportion", title = TeX("CL or PIP $\\geq 0.5$"))

p5 = ggplot(data, aes(x = method, y = num_CS95_eQTL / num_CS95)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8,
           position = position_dodge2()) +
  theme_classic() + custom_theme() + 
  scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
  guides(fill = "none") +
  labs(x = NULL, y = "Proportion", title = "95% CS")

p2345 = ggarrange(p2, p3, p4, p5, 
                  nrow = 2, ncol = 2, common.legend = T, legend = "none",
                  widths = c(1,1), heights = c(1,1), 
                  labels = c("b","c","d","e"),
                  font.label = list(size = 12, color = "black", face = "bold",
                                    family = NULL), align = "hv")

ggsave("eQTL_enrichment.pdf",
       ggarrange(p1, p2345, 
                 nrow = 1, ncol = 2, common.legend = T, legend = "bottom",
                 legend.grob = ggpubr::get_legend(p3),
                 widths = c(1,1.575), labels = c("a"),
                 font.label = list(size = 12, color = "black", face = "bold",
                                   family = NULL)),
       device = "pdf", width = 9, height = 5.92, units = "in", bg = "white")
