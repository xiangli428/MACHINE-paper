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
library(patchwork)
library(ggpubr)
library(ggh4x)
library(latex2exp)
library(scales)

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

coding = c("NSF","NSM","NSN","ASS","DSS")

registerDoParallel(22)

eQTLs = foreach(i = 1:22, .combine = "c") %dopar%
{
  unique(filter(read_delim(sprintf(
    "GTEx_eQTL_SuSiE_V10/chr%d.txt.gz", i)), startsWith(tissue, "Brain"))$rsid)
}

coding_SNPs = foreach(anno = coding, .combine = "c") %do%
{
  foreach(i = 1:22, .combine = "c") %dopar%
  {
    read_delim(sprintf(
      "dbSNP151/%s/chr%d.txt.gz", anno, i))$rsid
  }
}

coding_SNPs %<>% unique()

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
  
  idx1 = grep("EUR", CS_info$CS)
  idx2 = grep("EAS", CS_info$CS)
  
  output = data.frame(
    "method" = method,
    "num_variants" = sum(min_p$num_variants),
    "num_variants_0.5" = length(unique(variants_0.5$rsid)),
    "num_CS95" = nrow(CS_info),
    "CS95_size_mean" = mean(CS_info$CS.size),
    "CS95_size_sd" = sd(CS_info$CS.size) / sqrt(nrow(CS_info)),
    "num_CS95_EUR" = length(idx1),
    "CS95_size_mean_EUR" = mean(CS_info$CS.size[idx1]),
    "CS95_size_sd_EUR" = sd(CS_info$CS.size[idx1]) / sqrt(length(idx1)),
    "num_CS95_EAS" = length(idx2),
    "CS95_size_mean_EAS" = mean(CS_info$CS.size[idx2]),
    "CS95_size_sd_EAS" = sd(CS_info$CS.size[idx1]) / sqrt(length(idx2)),
    "num_variants_coding" = length(coding_SNPs),
    "num_variants_0.5_coding" = length(intersect(
      unique(variants_0.5$rsid), coding_SNPs)),
    "num_CS95_coding" = 0,
    "num_variants_eQTL" = length(eQTLs),
    "num_variants_0.5_eQTL" = length(intersect(
      unique(variants_0.5$rsid), eQTLs)),
    "num_CS95_eQTL" = 0)
  
  for(l in 1:nrow(CS_info))
  {
    SNPs = strsplit(CS_info$CS.SNP[l], ",")[[1]]
    
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

data$method %<>% factor(levels = names(methods))

data %<>% mutate(
  prop_variants_0.5_coding = num_variants_0.5_coding / num_variants_0.5,
  prop_variants_0.5_coding_sd = sqrt(prop_variants_0.5_coding * (
    1-prop_variants_0.5_coding) / num_variants_0.5),
  prop_CS95_coding = num_CS95_coding / num_CS95,
  prop_CS95_coding_sd = sqrt(prop_CS95_coding * (1-prop_CS95_coding) / num_CS95),
  prop_variants_0.5_eQTL = num_variants_0.5_eQTL / num_variants_0.5,
  prop_variants_0.5_eQTL_sd = sqrt(prop_variants_0.5_eQTL * (
    1-prop_variants_0.5_eQTL) / num_variants_0.5),
  prop_CS95_eQTL = num_CS95_eQTL / num_CS95,
  prop_CS95_eQTL_sd = sqrt(prop_CS95_eQTL * (1-prop_CS95_eQTL) / num_CS95))

saveRDS(data, "eQTL_enrichment.RData")

CS_overlap_num = readRDS("CS_overlap_num.RData")
CS_overlap_num %<>% gather(pops, nCSs, -method)
CS_overlap_num$pops %<>% gsub("_", " & ", .)
CS_overlap_num$pops %<>% factor(levels = c("EUR & EAS", "EUR-specific", "EAS-specific"))

CS_size = pivot_longer(data[,c(1,8,9,11,12)], -c(method), names_to = c(".value","pid"),
                       names_pattern = "(CS95_size_[a-z]+)_(.*)")
CS_size$pid %<>% factor(levels = c("EUR","EAS"))

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
    panel.spacing = unit(0.5, "in"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA))
}

colors = c("#521A13","#996035FF","#DA2222",
           "#FF9933FF","#FFDD00FF","#00AD00FF","#5D7A2BFF",
           "#009393FF","#69D2E7FF","#0066CCFF",
           "#1E2085FF","#8785B2FF","#953272FF")

p1 = ggplot(CS_overlap_num, aes(x = method, y = nCSs)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8,
           position = position_dodge2()) +
  geom_text(aes(label = nCSs), size = 2, vjust = -0.3) +
  facet_grid(pops ~ ., scales = "free_y") +
  theme_classic() + custom_theme() +
  theme(plot.background = element_rect(color = "black", linetype = "solid",
                                       linewidth = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.05,0.15))) +
  scale_fill_manual(values = colors[c(1,3,4,8,10:13)]) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Number of CSs", fill = "", title = NULL)

p2 = ggplot(CS_size, aes(x = method, y = CS95_size_mean)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8,
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = CS95_size_mean - CS95_size_sd, 
                    ymax = CS95_size_mean + CS95_size_sd),
                width = 0.8, position = position_dodge2(),
                linewidth = 0.25, color = "black") +
  facet_grid(. ~ pid, scales = "free_y") +
  theme_classic() + custom_theme() +
  theme(plot.background = element_rect(color = "black", linetype = "solid",
                                       linewidth = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.05,0.15))) +
  scale_fill_manual(values = colors[c(1,3,4,8,10:13)]) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Mean size of 95% CSs", fill = "", title = NULL)

p3 = ggplot(data, aes(x = method, y = prop_variants_0.5_coding, fill = method)) +
  geom_bar(stat = "identity", width = 0.8, position = position_dodge2()) +
  geom_errorbar(aes(ymin = prop_variants_0.5_coding - prop_variants_0.5_coding_sd, 
                    ymax = prop_variants_0.5_coding + prop_variants_0.5_coding_sd),
                width = 0.8, position = position_dodge2(),
                linewidth = 0.25, color = "black") +
  theme_classic() + custom_theme() + 
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  scale_fill_manual(values = colors[c(1,3,4,8,10:13)]) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Proportion", title = TeX("CL or PIP $\\geq 0.5$"))

p4 = ggplot(data, aes(x = method, y = prop_CS95_coding, fill = method)) +
  geom_bar(stat = "identity", width = 0.8, position = position_dodge2()) +
  geom_errorbar(aes(ymin = prop_CS95_coding - prop_CS95_coding_sd, 
                    ymax = prop_CS95_coding + prop_CS95_coding_sd),
                width = 0.8, position = position_dodge2(),
                linewidth = 0.25, color = "black") +
  theme_classic() + custom_theme() + 
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  scale_fill_manual(values = colors[c(1,3,4,8,10:13)]) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Proportion", title = "95% CS")

p5 = ggplot(data, aes(x = method, y = prop_variants_0.5_eQTL, fill = method)) +
  geom_bar(stat = "identity", width = 0.8, position = position_dodge2()) +
  geom_errorbar(aes(ymin = prop_variants_0.5_eQTL - prop_variants_0.5_eQTL_sd, 
                    ymax = prop_variants_0.5_eQTL + prop_variants_0.5_eQTL_sd),
                width = 0.8, position = position_dodge2(),
                linewidth = 0.25, color = "black") +
  theme_classic() + custom_theme() + 
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  scale_fill_manual(values = colors[c(1,3,4,8,10:13)]) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Proportion", title = TeX("CL or PIP $\\geq 0.5$"))

p6 = ggplot(data, aes(x = method, y = prop_CS95_eQTL, fill = method)) +
  geom_bar(stat = "identity", width = 0.8, position = position_dodge2()) +
  geom_errorbar(aes(ymin = prop_CS95_eQTL - prop_CS95_eQTL_sd, 
                    ymax = prop_CS95_eQTL + prop_CS95_eQTL_sd),
                width = 0.8, position = position_dodge2(),
                linewidth = 0.25, color = "black") +
  theme_classic() + custom_theme() + 
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  scale_fill_manual(values = colors[c(1,3,4,8,10:13)]) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Proportion", title = "95% CS")

layout = "
ABB
ACD
AEF
"

pdf("CS_num_eQTL_enrichment.pdf", width = 9, height = 9.2, useDingbats = F)
p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(design = layout, guides = "collect") &
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom",
        plot.margin = margin(l = 2, r = 2, b = 2, t = 2),
        plot.tag = element_text(size = 12, face = "bold"))
dev.off()
