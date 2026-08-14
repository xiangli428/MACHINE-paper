options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)
library(ggplot2)
library(cowplot)
library(ggpattern)
library(ggpubr)
library(ggh4x)
library(latex2exp)
library(scales)

setwd("simulation/results")

source("../simulation_setting.R")

data = readRDS("FDR_power.RData")

custom_theme = function()
{
  theme(
    aspect.ratio = 1,
    panel.spacing = unit(0.05, "in"),
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
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA))
}

colors = c("#521A13","#996035FF","#DA2222",
           "#FF9933FF","#FFDD00FF","#00AD00FF","#5D7A2BFF",
           "#009393FF","#69D2E7FF","#0066CCFF",
           "#1E2085FF","#8785B2FF","#953272FF")

p_list = list()

for(k in 1:2)
{
  sub = filter(data, LD == lds[k] & N2 == 2e5 & method %in% methods[8:13] &
                 task %in% pids & !(pops == "EUR+EAS" & task == "EUR"))
  
  p1 = ggplot(sub, aes(x = method, y = power_0.9)) +
    geom_bar_pattern(aes(fill = method, pattern = task), 
                     pattern_fill = "gray", pattern_color = "gray",
                     pattern_size = 0, pattern_density = 1/3,
                     stat = "identity", width = 0.8, 
                     position = position_dodge2(padding = 0.2), show.legend = T) +
    geom_errorbar(aes(ymin = power_0.9 - power_0.9_sd, 
                      ymax = power_0.9 + power_0.9_sd),
                  width = 0.8, position = position_dodge2(padding = 0.2), 
                  linewidth = 0.25, color = "black") +
    facet_grid(. ~ scenario, scales = "free",
               labeller = labeller(scenario = scenarios2)) +
    theme_classic() + custom_theme() +
    scale_y_continuous(limits = c(0,0.35)) +
    scale_fill_manual(values = colors[8:13], limits = methods[8:13],
                      drop = F) +
    scale_pattern_manual(values = c("none", "circle", "stripe")) +
    guides(fill = guide_legend(nrow = 1, override.aes = list(pattern = "none"))) +
    labs(x = NULL, y = TeX("Power at CL or PIP \\geq 0.9"), fill = "")
  
  p2 = ggplot(sub, aes(x = method, y = size_CS95)) +
    geom_bar_pattern(aes(fill = method, pattern = task), 
                     pattern_fill = "gray", pattern_color = "gray",
                     pattern_size = 0, pattern_density = 1/3,
                     stat = "identity", width = 0.8, 
                     position = position_dodge2(padding = 0.2), show.legend = T) +
    geom_errorbar(aes(ymin = size_CS95 - size_CS95_sd, 
                      ymax = size_CS95 + size_CS95_sd),
                  width = 0.8, position = position_dodge2(padding = 0.2), 
                  linewidth = 0.25, color = "black") +
    facet_grid(. ~ scenario, scales = "free",
               labeller = labeller(scenario = scenarios2)) +
    theme_classic() + custom_theme() +
    scale_y_continuous(limits = c(0,14)) +
    scale_fill_manual(values = colors[8:13], limits = methods[8:13],
                      drop = F) +
    scale_pattern_manual(values = c("none", "circle", "stripe")) +
    guides(fill = guide_legend(nrow = 1, override.aes = list(pattern = "none"))) +
    labs(x = NULL, y = "Size of 95% CSs", fill = "")
  
  p3 = ggplot(sub, aes(x = method, y = power_CS95)) +
    geom_bar_pattern(aes(fill = method, pattern = task), 
                     pattern_fill = "gray", pattern_color = "gray",
                     pattern_size = 0, pattern_density = 1/3,
                     stat = "identity", width = 0.8, 
                     position = position_dodge2(padding = 0.2), show.legend = T) +
    geom_errorbar(aes(ymin = power_CS95 - power_CS95_sd, 
                      ymax = power_CS95 + power_CS95_sd),
                  width = 0.8, position = position_dodge2(padding = 0.2), 
                  linewidth = 0.25, color = "black") +
    facet_grid(. ~ scenario, scales = "free",
               labeller = labeller(scenario = scenarios2)) +
    theme_classic() + custom_theme() +
    scale_y_continuous(limits = c(0,0.5)) +
    scale_fill_manual(values = colors[8:13], limits = methods[8:13],
                      drop = F) +
    scale_pattern_manual(values = c("none", "circle", "stripe")) +
    guides(fill = guide_legend(nrow = 1, override.aes = list(pattern = "none"))) +
    labs(x = NULL, y = "Power of 95% CSs", fill = "")
  
  p_list[[k]] = annotate_figure(ggarrange(
    p1, p2, p3, nrow = 3, ncol = 1, common.legend = T, legend = "none",
    heights = c(1,1,1), labels = list(c("a","b","c"),c("d","e","f"))[[k]],
    font.label = list(size = 12, color = "black", face = "bold", 
                      family = NULL), align = "v"),
    top = text_grob(lds[k], face = "bold", size = 12))
}

pdf("single_ancestry.pdf", width = 9, height = 6, useDingbats = F, bg = "white")
ggarrange(plotlist = p_list, nrow = 1, ncol = 2, 
          common.legend = T, legend = "bottom", legend.grob = get_legend(p1),
  heights = c(1,1,1), align = "hv")
dev.off()
