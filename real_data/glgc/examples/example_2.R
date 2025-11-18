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
library(gggenes)
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

db = "GLGC"
pheno = "LDL"
i = 19
block = 38

data_all = foreach(method = names(methods)[c(1,4)]) %do%
{
  if(method %in% c("MACHINE + g-LDSC", "MACHINE"))
  {
    read.delim(sprintf("%s/%s/%s/chr%d/%d/variant_info.txt.gz", 
                       methods[method], pheno, db, i, block))
  } else {
    foreach(pid = pids) %do%
    {
      read.delim(sprintf("%s/%s/%s/%s/chr%d/%d/variant_info.txt.gz", 
                         methods[method], pheno, db, pid, i, block))
    }
  }
}
names(data_all[[2]]) = pids

gene_position = read.delim(
  "~/Documents/GWAS/data/reference/human/v42/GRCh37/gene_position.hg19.v42.bed",
  header = F) %>% filter(V1 == sprintf("chr%d",i) & !(
    (V3 < min(data_all[[1]]$position)) | (V2 > max(data_all[[1]]$position))) &
      V7 == "protein_coding")

gene_table = data.frame("type" = "",
                        "start" = gene_position$V2,
                        "end" = gene_position$V3,
                        "name" = gene_position$V4,
                        "orientation" = c("+" = T, "-" = F)[gene_position$V6],
                        stringsAsFactors = F)
gene_table$start %<>% pmax(min(data_all[[1]]$position))
gene_table$end %<>% pmin(max(data_all[[1]]$position))

data_all[[1]]$CS.EUR[which(data_all[[1]]$CS.EUR == "CS:EUR-19-38-1")] = "CS-1"
data_all[[1]]$CS.EUR[which(data_all[[1]]$CS.EUR == "CS:EUR-19-38-2")] = "CS-2"
data_all[[1]]$CS.EUR[which(data_all[[1]]$CS.EUR == "CS:EUR-19-38-3")] = "CS-3"
data_all[[1]]$CS.AFR[which(data_all[[1]]$CS.AFR == "CS:AFR-19-38-2")] = "CS-3"
data_all[[1]]$CS.AFR[which(data_all[[1]]$CS.AFR == "CS:AFR-19-38-1")] = "CS-4"
data_all[[1]]$CS.EAS[which(data_all[[1]]$CS.EAS == "CS:EAS-19-38-1")] = "CS-3"
data_all[[2]][[1]]$CS[which(data_all[[2]][[1]]$CS == "CS:EUR-19-38-1")] = "CS-1"
data_all[[2]][[1]]$CS[which(data_all[[2]][[1]]$CS == "CS:EUR-19-38-2")] = "CS-2"
data_all[[2]][[2]]$CS[which(data_all[[2]][[2]]$CS == "CS:AFR-19-38-2")] = "CS-3"
data_all[[2]][[2]]$CS[which(data_all[[2]][[2]]$CS == "CS:AFR-19-38-1")] = "CS-4"
data_all[[2]][[3]]$CS[which(data_all[[2]][[3]]$CS == "CS:EAS-19-38-1")] = "CS-3"

data_all[[1]]$CS.EUR %<>% factor(levels = sprintf("CS-%s",1:4))
data_all[[1]]$CS.AFR %<>% factor(levels = sprintf("CS-%s",1:4))
data_all[[1]]$CS.EAS %<>% factor(levels = sprintf("CS-%s",1:4))
data_all[[2]][[1]]$CS %<>% factor(levels = sprintf("CS-%s",1:4))
data_all[[2]][[2]]$CS %<>% factor(levels = sprintf("CS-%s",1:4))
data_all[[2]][[3]]$CS %<>% factor(levels = sprintf("CS-%s",1:4))

save(data_all, gene_table, file = "examples/example_2.RData")

load("examples/example_2.RData")

custom_theme = function()
{
  theme(
    aspect.ratio = 1/4,
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.25),
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

data_1 = pivot_longer(data_all[[1]][,c(2,5,11,15,19,22,24,26,28,30,32)],
                      -c(position, rsid), names_to = c(".value","pid"),
                      names_pattern = "([A-Za-z]+[0-9]*_?[A-Za-z]+)\\.(.*)") %>%
  filter(!is.na(neglog10_pval))
data_1$pid %<>% factor(levels = pids)
data_1 %<>% arrange(!is.na(CS), CS, pid, position)

data_2 = foreach(pid = pids, .combine = "rbind") %do%
{
  mutate(data_all[[2]][[pid]][,c(2,5,23,25)], pid = pid)
}
data_2$pid %<>% factor(levels = pids)
data_2 %<>% arrange(!is.na(CS), CS, pid, position)


p1 = ggplot(data_1, aes(x = position / 1e6, y = neglog10_pval, color = CS)) +
  geom_point(size = 0.5) +
  facet_grid(pid ~ "Marginal association", scales = "free_y") +
  theme_classic() + custom_theme() + 
  scale_x_continuous(limits = range(data_all[[1]]$position) / 1e6) +
  scale_color_manual(breaks = sprintf("CS-%s",1:4), values = hue_pal()(4)) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x = "Position on chr19 (Mb)", y = TeX("$-log_{10}(\\italic(P))$"), 
       color = "")

p2 = ggplot(data_1, aes(x = position / 1e6, y = CL, color = CS)) +
  geom_point(size = 0.5) +
  facet_grid(pid ~ "MACHINE + g-LDSC") +
  theme_classic() + custom_theme() + 
  scale_x_continuous(limits = range(data_all[[1]]$position) / 1e6) +
  scale_y_continuous(limits = c(0,1)) +
  scale_color_manual(breaks = sprintf("CS-%s",1:4), values = hue_pal()(4)) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x = "Position on chr19 (Mb)", y = "CL", color = "")

p3 = ggplot(data_2, aes(x = position / 1e6, y = CL, color = CS)) +
  geom_point(size = 0.5) +
  facet_grid(pid ~ "h2-D2 + g-LDSC") +
  theme_classic() + custom_theme() + 
  scale_x_continuous(limits = range(data_all[[1]]$position) / 1e6) +
  scale_y_continuous(limits = c(0,1)) +
  scale_color_manual(breaks = sprintf("CS-%s",1:4), values = hue_pal()(4)) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x = "Position on chr19 (Mb)", y = "CL", color = "")

gene_table$y_adjust = rep_len(c(1.2,1.2,1.2,0.8,0.8,0.8), 
                              length.out = nrow(gene_table))
gene_table$group = rep_len(1:3, length.out = nrow(gene_table))

p4 = ggplot(gene_table, aes(
  xmin = start / 1e6, xmax = end / 1e6, y = type, forward = orientation)) +
  geom_gene_arrow(
    arrowhead_height = unit(1, "mm"), arrowhead_width = unit(1, "mm"), 
    arrow_body_height = unit(1, "mm"), fill = "lightgreen") +
  geom_text(aes(x = (start+end)/2e6, y = y_adjust, label = name), size = 2.5) +
  facet_grid(group ~ "Gene") +
  theme_genes() + 
  theme(aspect.ratio = 1/4,
        panel.spacing = unit(0.05, "in"),
        axis.text.x = element_text(size = 7),
        axis.title.x = element_text(size = 8),
        axis.ticks.x = element_line(linewidth = 0.25),
        axis.line.x = element_line(linewidth = 0.25)) +
  scale_x_continuous(limits = range(data_all[[1]]$position) / 1e6) +
  labs(x = "Position on chr19 (Mb)", y = NULL)

ggsave("examples/example_2.pdf",
       ggarrange(p1, p2, p3, p4,
                 nrow = 2, ncol = 2, 
                 align = "hv", widths = c(1,1), heights = c(1,1),
                 common.legend = T, legend = "bottom", 
                 labels = c("a","b","c","d"),
                 font.label = list(size = 12, color = "black", face = "bold",
                                   family = NULL)),
       device = "pdf", width = 9, height = 7.5, units = "in", bg = "white")
