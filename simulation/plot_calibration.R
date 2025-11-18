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
library(latex2exp)
library(scales)

setwd("simulation/results")

source("../simulation_setting.R")

groups = c("[0,0.1)","[0.1,0.5)","[0.5,0.9)","[0.9,1]")

calibration_table = foreach(s = 1:3, .combine = "rbind") %do% 
{
  foreach(m = 1:2, .combine = "rbind") %do% 
  {
    foreach(ld = lds, .combine = "rbind") %do% 
    {
      foreach(method = names(ld_methods[[ld]]), .combine = "rbind") %do% 
      {
        calibration = readRDS(sprintf(
          "setting_%s/%s/calibration_N2-%d.RData",
          s, ld_methods[[ld]][method], N2_seq[m]))
        
        foreach(pop = names(pops), .combine = "rbind") %do%
        {
          data.frame("LD" = ld,
                     "scenario" = s,
                     "N2" = N2_seq[m],
                     "POP" = pops[pop],
                     "method" = method,
                     "group" = groups,
                     "nVar" = calibration[[pop]]$n,
                     "Expected" = calibration[[pop]]$Expected,
                     "Prop" = calibration[[pop]]$Prop)
        }
      }
    }
  }
}

calibration_table$LD %<>% factor(levels = lds)
calibration_table$POP %<>% factor(levels = pops)
calibration_table$method %<>% factor(levels = methods)
calibration_table$group %<>% factor(levels = groups)
calibration_table %<>% arrange(LD, scenario, N2, POP, group, method)

calibration_table %<>% mutate(Prop_sd = sqrt(Prop * (1-Prop) / nSNPs))

write_delim(calibration_table, "calibration_table.txt", delim = '\t')
saveRDS(calibration_table, "calibration_table.RData")

calibration_table = readRDS("calibration_table.RData")

custom_theme = function()
{
  theme(
    aspect.ratio = 1/3.2,
    axis.text = element_text(size = 7),  
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

calibration_table$N2 = sprintf("N^{(2)}==%d*k", calibration_table$N2 / 1e3)
calibration_table$N2 %<>% factor(levels = sprintf("N^{(2)}==%d*k", N2_seq / 1e3))


nmed = length(ld_methods[[1]])
data = filter(calibration_table, LD == lds[1])
data$x = as.numeric(data$group) + seq(
  from = -0.4, by = 0.8/nmed, length.out = nmed)
data$xend = as.numeric(data$group) + seq(
  to = 0.4, by = 0.8/nmed, length.out = nmed)

pdf("calibration_UKBB.pdf", width = 9, height = 8.4, bg = "white")
ggplot(data, aes(x = group, y = Prop, color = method)) +
  geom_point(position = position_dodge2(width = 0.8), size = 0.5) +
  geom_errorbar(aes(ymin = Prop - Prop_sd, ymax = Prop + Prop_sd),
                position = position_dodge2(width = 0.8), 
                width = 0.8, linewidth = 0.1) +
  geom_segment(aes(x = x, xend = xend, y = Expected, yend = Expected),
               color = "black", alpha = 0.5) +
  facet_grid(POP + N2 ~ scenario, scales = "free", labeller = 
               labeller(scenario = scenarios, N2 = label_parsed)) +
  theme_classic() + custom_theme() + 
  scale_y_continuous(limits = c(0,1)) +
  scale_color_manual(values = hue_pal()(13)[c(1:11,13)]) +
  guides(color = guide_legend(nrow = 2)) +
  labs(x = "CL or PIP bins", y = "Proportion of causal variants", colour = "")
dev.off()


nmed = length(ld_methods[[2]])
data = filter(calibration_table, LD == lds[2])
data$x = as.numeric(data$group) + seq(
  from = -0.4, by = 0.8/nmed, length.out = nmed)
data$xend = as.numeric(data$group) + seq(
  to = 0.4, by = 0.8/nmed, length.out = nmed)

pdf("calibration_1kG.pdf", width = 9, height = 8.4, bg = "white")
ggplot(data, aes(x = group, y = Prop, color = method)) +
  geom_point(position = position_dodge2(width = 0.8), size = 0.5) +
  geom_errorbar(aes(ymin = Prop - Prop_sd, ymax = Prop + Prop_sd),
                position = position_dodge2(width = 0.8), 
                width = 0.8, linewidth = 0.1) +
  geom_segment(aes(x = x, xend = xend, y = Expected, yend = Expected),
               color = "black", alpha = 0.5) +
  facet_grid(POP + N2 ~ scenario, scales = "free", labeller = 
               labeller(scenario = scenarios, N2 = label_parsed)) +
  theme_classic() + custom_theme() + 
  scale_y_continuous(limits = c(0,1)) +
  scale_color_manual(values = hue_pal()(13)[c(1:3,8:10,12,13)]) +
  guides(color = guide_legend(nrow = 2)) +
  labs(x = "CL or PIP bins", y = "Proportion of causal variants", colour = "")
dev.off()
