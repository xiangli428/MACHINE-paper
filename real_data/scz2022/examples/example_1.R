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

i = 12
block = 139

data_all = foreach(method = names(methods)[c(1,4)]) %do%
{
  if(method %in% c("MACHINE + g-LDSC", "MACHINE"))
  {
    read.delim(sprintf("%s/chr%d/%d/variant_info.txt.gz", 
                       methods[method], i, block))
  } else {
    foreach(pid = pids) %do%
    {
      read.delim(sprintf("%s/%s/chr%d/%d/variant_info.txt.gz", 
                         methods[method], pid, i, block))
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
gene_table$name[9] = "F5H7X1"
gene_table$name[16] = "LOC128125816"
gene_table$name[18] = "H3BMM5"

data_all[[1]]$CS.EUR[which(data_all[[1]]$CS.EUR == "CS:EUR-12-139-1")] = "CS-1"
data_all[[1]]$CS.EUR[which(data_all[[1]]$CS.EUR == "CS:EUR-12-139-2")] = "CS-2"
data_all[[1]]$CS.EUR[which(data_all[[1]]$CS.EUR == "CS:EUR-12-139-3")] = "CS-3"
data_all[[1]]$CS.EAS[which(data_all[[1]]$CS.EAS == "CS:EAS-12-139-1")] = "CS-1"
data_all[[1]]$CS.EAS[which(data_all[[1]]$CS.EAS == "CS:EAS-12-139-2")] = "CS-2"

data_all[[2]][[1]]$CS[which(data_all[[2]][[1]]$CS == "CS:EUR-12-139-1")] = "CS-1"
data_all[[2]][[1]]$CS[which(data_all[[2]][[1]]$CS == "CS:EUR-12-139-2")] = "CS-3"
data_all[[2]][[2]]$CS[which(data_all[[2]][[2]]$CS == "CS:EAS-12-139-1")] = "CS-1"

data_all[[1]]$CS.EUR %<>% factor(levels = sprintf("CS-%s",1:3))
data_all[[1]]$CS.EAS %<>% factor(levels = sprintf("CS-%s",1:3))
data_all[[2]][[1]]$CS %<>% factor(levels = sprintf("CS-%s",1:3))
data_all[[2]][[2]]$CS %<>% factor(levels = sprintf("CS-%s",1:3))

save(data_all, gene_table, file = "examples/example_1.RData")
load("examples/example_1.RData")

custom_theme = function()
{
  theme(
    aspect.ratio = 1/3,
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

data_1 = pivot_longer(data_all[[1]][,c(2,5,10,13,15,17,19,21)],
                      -c(position, rsid), names_to = c(".value","pid"),
                      names_pattern = "([A-Za-z]+_?[A-Za-z]+)\\.(.*)") %>%
  filter(!is.na(PVAL))
data_1$pid %<>% factor(levels = pids)
data_1 %<>% arrange(!is.na(CS), CS, pid, position)

data_2 = foreach(pid = pids, .combine = "rbind") %do%
{
  mutate(data_all[[2]][[pid]][,c(2,5,10,11,13)], pid = pid)
}
data_2$pid %<>% factor(levels = pids)
data_2 %<>% arrange(!is.na(CS), CS, pid, position)

p1 = ggplot(data_1, aes(x = position / 1e6, y = -log10(PVAL), color = CS)) +
  geom_point(size = 0.5) +
  facet_grid(pid ~ "Marginal association", scales = "free_y") +
  theme_classic() + custom_theme() + 
  scale_x_continuous(limits = range(data_all[[1]]$position) / 1e6) +
  scale_color_manual(breaks = sprintf("CS-%s",1:3), values = hue_pal()(3)) +
  guides(color = guide_legend(nrow = 1)) +
  labs(x = "Position on chr12 (Mb)", y = TeX("$-log_{10}(\\italic(P))$"), 
       color = "")

p2 = ggplot(data_1, aes(x = position / 1e6, y = CL, color = CS)) +
  geom_point(size = 0.5) +
  facet_grid(pid ~ "MACHINE + g-LDSC") +
  theme_classic() + custom_theme() + 
  scale_x_continuous(limits = range(data_all[[1]]$position) / 1e6) +
  scale_y_continuous(limits = c(0,0.5)) +
  scale_color_manual(breaks = sprintf("CS-%s",1:3), values = hue_pal()(3)) +
  guides(color = "none") +
  labs(x = "Position on chr12 (Mb)", y = "CL", color = "")

p3 = ggplot(data_2, aes(x = position / 1e6, y = CL, color = CS)) +
  geom_point(size = 0.5) +
  facet_grid(pid ~ "h2-D2 + g-LDSC") +
  theme_classic() + custom_theme() + 
  scale_x_continuous(limits = range(data_all[[1]]$position) / 1e6) +
  scale_y_continuous(limits = c(0,0.5)) +
  scale_color_manual(breaks = sprintf("CS-%s",1:3), values = hue_pal()(3)) +
  guides(color = "none") +
  labs(x = "Position on chr12 (Mb)", y = "CL", color = "")

gene_table$y_adjust = rep_len(c(rep(1.25,5),rep(0.75,5)), 
                              length.out = nrow(gene_table))
gene_table$group = rep_len(1:5, length.out = nrow(gene_table))

p4 = ggplot(gene_table, aes(
  xmin = start / 1e6, xmax = end / 1e6, y = type, forward = orientation)) +
  geom_gene_arrow(
    arrowhead_height = unit(1, "mm"), arrowhead_width = unit(1, "mm"), 
    arrow_body_height = unit(1, "mm"), fill = "lightgreen") +
  geom_text(aes(x = (start+end)/2e6, y = y_adjust, label = name), size = 2) +
  facet_grid(group ~ "Gene") +
  theme_genes() + 
  theme(axis.text.x = element_text(size = 7),
        axis.title.x = element_text(size = 8),
        axis.ticks.x = element_line(linewidth = 0.25),
        axis.line.x = element_line(linewidth = 0.25)) +
  scale_x_continuous(limits = range(data_all[[1]]$position) / 1e6) +
  labs(x = "Position on chr12 (Mb)", y = NULL)

a1234 = align_plots(p1 + theme(legend.position = "none"), p2, p3, p4, 
                    align = "hv", axis = "tblr")

ggsave("examples/example_1.pdf",
       ggarrange(plot_grid(plotlist = a1234,
                           nrow = 2, ncol = 2,
                           labels = c('a','b','c','d'), label_size = 12),
                 common.legend = T, legend = "bottom",
                 legend.grob = ggpubr::get_legend(p1)),
       device = "pdf", width = 9, height = 6.9, units = "in", bg = "white")
