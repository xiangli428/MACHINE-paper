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
library(ggrepel)
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
pheno = "TC"
i = 8
block = 125

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
  (V3 < min(data_all[[1]]$position)) | (V2 > max(data_all[[1]]$position))))

gene_table = data.frame("type" = "",
                        "start" = gene_position$V2,
                        "end" = gene_position$V3,
                        "name" = gene_position$V4,
                        "orientation" = c("+" = T, "-" = F)[gene_position$V6],
                        stringsAsFactors = F)
gene_table$start %<>% pmax(min(data_all[[1]]$position))
gene_table$end %<>% pmin(max(data_all[[1]]$position))
gene_table$name[3] = "TRIB1AL"

data_all[[1]]$CS.EUR[which(data_all[[1]]$CS.EUR == "CS:EUR-8-125-1")] = "CS-1"
data_all[[1]]$CS.EUR[which(data_all[[1]]$CS.EUR == "CS:EUR-8-125-2")] = "CS-2"
data_all[[1]]$CS.AFR[which(data_all[[1]]$CS.AFR == "CS:AFR-8-125-1")] = "CS-1"
data_all[[1]]$CS.AFR[which(data_all[[1]]$CS.AFR == "CS:AFR-8-125-2")] = "CS-2"
data_all[[1]]$CS.EAS[which(data_all[[1]]$CS.EAS == "CS:EAS-8-125-1")] = "CS-1"
data_all[[1]]$CS.EAS[which(data_all[[1]]$CS.EAS == "CS:EAS-8-125-2")] = "CS-2"
data_all[[2]][[1]]$CS[which(data_all[[2]][[1]]$CS == "CS:EUR-8-125-1")] = "CS-3"
data_all[[2]][[1]]$CS[which(data_all[[2]][[1]]$CS == "CS:EUR-8-125-2")] = "CS-4"
data_all[[2]][[2]]$CS[which(data_all[[2]][[2]]$CS == "CS:AFR-8-125-1")] = "CS-1"
data_all[[2]][[3]]$CS[which(data_all[[2]][[3]]$CS == "CS:EAS-8-125-1")] = "CS-1"
data_all[[2]][[3]]$CS[which(data_all[[2]][[3]]$CS == "CS:EAS-8-125-2")] = "CS-2"

data_all[[1]]$CS.EUR %<>% factor(levels = sprintf("CS-%s",1:4))
data_all[[1]]$CS.AFR %<>% factor(levels = sprintf("CS-%s",1:4))
data_all[[1]]$CS.EAS %<>% factor(levels = sprintf("CS-%s",1:4))
data_all[[2]][[1]]$CS %<>% factor(levels = sprintf("CS-%s",1:4))
data_all[[2]][[2]]$CS %<>% factor(levels = sprintf("CS-%s",1:4))
data_all[[2]][[3]]$CS %<>% factor(levels = sprintf("CS-%s",1:4))

label_SNPs = c("rs2980888","rs8180991","rs28601761","rs112875651")

R = foreach(pid = pids) %do%
{
  readMM(sprintf("LD/%s/%s/chr%d/%d/LD.mtx.gz", pheno, pid, i, block))
}

R_df = foreach(k = 1:3, .combine = "rbind") %do%
{
  diag(R[[k]]) = 1
  pid = pids[k]
  sub = data.frame("pid" = pid, "i" = R[[k]]@i+1, "j" = R[[k]]@j+1, 
                   "r" = R[[k]]@x)
  sub$i = data_all[[1]]$rsid[!is.na(data_all[[1]][,sprintf("Z.%s",pid)])][sub$i]
  sub$j = data_all[[1]]$rsid[!is.na(data_all[[1]][,sprintf("Z.%s",pid)])][sub$j]
  filter(sub, i %in% label_SNPs & j %in% label_SNPs)
}
idx = which(match(R_df$i, label_SNPs) > match(R_df$j, label_SNPs))
R_df[idx,c(2,3)] = R_df[idx,c(3,2)]

R_df$pid %<>% factor(levels = pids)
R_df$i %<>% factor(levels = label_SNPs)
R_df$j %<>% factor(levels = label_SNPs)

save(data_all, gene_table, R_df, label_SNPs, file = "examples/example_1.RData")

load("examples/example_1.RData")

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

data_1 %<>% mutate(alpha = if_else(!is.na(CS) | rsid %in% label_SNPs, 1, 0.25))
label_df = filter(data_1, rsid %in% label_SNPs)
label_df %<>% arrange(pid, position)
label_df$nudge_x = c(0,-0.05,0,0.02,-0.05,0,0.05,0.02,-0.05,0,0.02,0.02)
label_df$nudge_y = c(-100,50,-100,0,0,2,2,-1,-2,5,5,5)

p1 = ggplot(data_1, aes(x = position / 1e6, y = neglog10_pval, color = CS)) +
  geom_point(aes(alpha = alpha), size = 0.5) +
  geom_text_repel(data = label_df, aes(
    x = position / 1e6, y = neglog10_pval, label = rsid), color = "black",
    size = 2, force = 50, force_pull = 0, max.time = 5, max.iter = Inf,
    nudge_x = label_df$nudge_x, nudge_y = label_df$nudge_y, 
    box.padding = 0.5, point.size = 0.5, min.segment.length = 0,
    segment.alpha = 0.5, na.rm = T, seed = 2, show.legend = F, verbose = T) +
  facet_grid(pid ~ "Marginal association", scales = "free_y") +
  theme_classic() + custom_theme() + 
  scale_x_continuous(limits = range(data_all[[1]]$position) / 1e6) +
  scale_alpha_continuous(limits = c(0,1)) +
  scale_color_manual(breaks = sprintf("CS-%s",1:4), values = hue_pal()(4)) +
  guides(color = "none", alpha = "none") +
  labs(x = "Position on chr8 (Mb)", y = TeX("$-log_{10}(\\italic(P))$"), 
       color = "")

label_df$nudge_x = rep(c(-0.05,-0.05,0.05,0.05),3)
label_df$nudge_y = rep(c(0.5,0,0.5,0),3)

p2 = ggplot(data_1, aes(x = position / 1e6, y = CL, color = CS)) +
  geom_point(aes(alpha = alpha), size = 0.5) +
  geom_text_repel(data = label_df, aes(
    x = position / 1e6, y = CL, label = rsid), color = "black",
    size = 2, force = 50, force_pull = 0.1, max.time = 5, max.iter = Inf,
    nudge_x = label_df$nudge_x, nudge_y = label_df$nudge_y, 
    box.padding = 0.5, point.size = 0.5, min.segment.length = 0, 
    segment.alpha = 0.5, na.rm = T, seed = 2, show.legend = F, verbose = T) +
  facet_grid(pid ~ "MACHINE + g-LDSC") +
  theme_classic() + custom_theme() + 
  scale_x_continuous(limits = range(data_all[[1]]$position) / 1e6) +
  scale_y_continuous(limits = c(0,1)) +
  scale_alpha_continuous(limits = c(0,1)) +
  scale_color_manual(breaks = sprintf("CS-%s",1:4), values = hue_pal()(4)) +
  guides(color = "none", alpha = "none") +
  labs(x = "Position on chr8 (Mb)", y = "CL", color = "")

data_2 %<>% mutate(alpha = if_else(!is.na(CS) | rsid %in% label_SNPs, 1, 0.25))
label_df = filter(data_2, rsid %in% label_SNPs)
label_df %<>% arrange(pid, position)
label_df$nudge_x = rep(c(-0.05,0.05,0.05,0.05),3)
label_df$nudge_y = c(0,0.6,0.6,0.6,rep(0.5,4),0.6,0.6,0.6,0.1)

p3 = ggplot(data_2, aes(x = position / 1e6, y = CL, color = CS)) +
  geom_point(aes(alpha = alpha), size = 0.5) +
  geom_text_repel(data = label_df, aes(
    x = position / 1e6, y = CL, label = rsid), color = "black",
    size = 2, force = 50, force_pull = 0.1, max.time = 5, max.iter = Inf,
    nudge_x = label_df$nudge_x, nudge_y = label_df$nudge_y,
    box.padding = 0.5, point.size = 0.5, min.segment.length = 0, 
    segment.alpha = 0.5, na.rm = T, seed = 2, show.legend = F, verbose = T) +
  facet_grid(pid ~ "h2-D2 + g-LDSC") +
  theme_classic() + custom_theme() + 
  scale_x_continuous(limits = range(data_all[[1]]$position) / 1e6) +
  scale_y_continuous(limits = c(0,1)) +
  scale_alpha_continuous(limits = c(0,1)) +
  scale_color_manual(breaks = sprintf("CS-%s",1:4), values = hue_pal()(4)) +
  guides(color = guide_legend(nrow = 1), alpha = "none") +
  labs(x = "Position on chr8 (Mb)", y = "CL", color = "")

gene_table$y_adjust = rep_len(c(1.2,0.8), length.out = nrow(gene_table))

p4 = ggplot(gene_table, aes(
  xmin = start / 1e6, xmax = end / 1e6, y = type, forward = orientation)) +
  geom_gene_arrow(
    arrowhead_height = unit(1, "mm"), arrowhead_width = unit(1, "mm"), 
    arrow_body_height = unit(1, "mm"), fill = "lightgreen") +
  geom_text(aes(x = (start+end)/2e6, y = y_adjust, label = name), size = 2.5) +
  theme_genes() + 
  theme(aspect.ratio = 1/4,
        panel.spacing = unit(0.05, "in"),
        axis.text = element_text(size = 7),
        axis.title.x = element_text(size = 8),
        axis.ticks.x = element_line(linewidth = 0.25),
        axis.line.x = element_line(linewidth = 0.25)) +
  scale_x_continuous(limits = range(data_all[[1]]$position) / 1e6) +
  labs(x = "Position on chr8 (Mb)", y = "")

p5 = ggplot(R_df, aes(x = i, y = j)) + 
  geom_tile(aes(fill = r)) + coord_equal() +
  geom_text(aes(label = sprintf("%0.2f", r)), size = 2) +
  facet_grid(. ~ pid) + theme_classic() +
  theme(
    aspect.ratio = 1,
    panel.spacing = unit(0.05, "in"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.text.x = element_text(size = 9),
    strip.text.y = element_blank(),
    strip.background = element_rect(
      fill = "lightgray", color = "lightgray"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  scale_y_discrete(limits = rev(label_SNPs)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1,1)) +
  guides(fill = "none")

a123 = align_plots(p1, p2, p3 + theme(legend.position = "none"), 
                   align = "hv", axis = "tblr")

a34 = align_plots(a123[[3]], p4, 
                  align = "v", axis = "lr")

p45 = plot_grid(a34[[2]], p5,
                labels = c('d','e'), label_size = 12, nrow = 2, ncol = 1)

ggsave("examples/example_1.pdf",
       ggarrange(plot_grid(a123[[1]], a123[[2]], a34[[1]], p45, 
                           nrow = 2, ncol = 2,
                           labels = c('a','b','c',''), label_size = 12),
                 common.legend = T, legend = "bottom",
                 legend.grob = ggpubr::get_legend(p3)),
       device = "pdf", width = 9, height = 7.5, units = "in", bg = "white")
