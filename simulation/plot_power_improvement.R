options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(ggpattern)
library(latex2exp)
library(scales)

setwd("simulation/results")

source("../simulation_setting.R")

data = readRDS("FDR_power.RData")
data %<>% filter(is.element(method, methods[c(1:3,8:10)]))

power = foreach(ld = lds, .combine = "rbind") %do%
{
  foreach(s = 1:3, .combine = "rbind") %do%
  {
    foreach(k = 2:3, .combine = "rbind") %do%
    {
      foreach(m = 1:2, .combine = "rbind") %do%
      {
        foreach(pid = names(tasks[[k-1]]), .combine = "rbind") %do%
        {
          sub = filter(data, LD == ld & scenario == s & 
                         pops == sprintf("EUR+%s", pids[k]) &
                         N2 == N2_seq[m] & task == tasks[[k-1]][pid])
          data.frame("LD" = ld,
                     "scenario" = s,
                     "pops" = sprintf("EUR+%s", pids[k]),
                     "N2" = N2_seq[m],
                     "task" = tasks[[k-1]][pid],
                     "anno_method" = c("g-LDSC", "PolyFun", "none"),
                     "h2-D2_0.9" = sub$power_0.9[4:6],
                     "MACHINE_0.9" = sub$power_0.9[1:3],
                     "improve_0.9" = 0,
                     "h2-D2_CS95" = sub$power_CS95[4:6],
                     "MACHINE_CS95" = sub$power_CS95[1:3],
                     "improve_CS95" = 0,
                     check.names = F)
        }
      }
    }
  }
}

power %<>% mutate(improve_0.9 = (MACHINE_0.9 - `h2-D2_0.9`) / `h2-D2_0.9`,
                  improve_CS95 = (MACHINE_CS95 - `h2-D2_CS95`) / `h2-D2_CS95`)
power %<>% mutate(improve_0.9 = MACHINE_0.9 - `h2-D2_0.9`,
                  improve_CS95 = MACHINE_CS95 - `h2-D2_CS95`)
power$pops %<>% factor(levels = c("EUR+AFR","EUR+EAS"))
power$task %<>% factor(levels = c("Cross", "EUR", "AFR", "EAS", "Shared"))
power$anno_method %<>% factor(levels = c("g-LDSC", "PolyFun", "none"))

saveRDS(power, "power_improvement.RData")

power$task2 = as.character(power$task)
power$task2[power$task2 %in% c("AFR","EAS")] = "AFR / EAS"
power$task2 %<>% factor(levels = c("Cross", "EUR", "AFR / EAS", "Shared"))

custom_theme = function()
{
  theme(
    aspect.ratio = 1/3,
    panel.spacing = unit(0.05, "in"),
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
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA))
}

pdf("power_improvement.pdf", width = 9, height = 9.06, bg = "white")
for(k in 1:2)
{
  p1 = ggplot(filter(power, LD == lds[k]), aes(
    x = factor(sprintf("%dk", N2/1e3), levels = N2_label), y = improve_0.9)) +
    geom_bar_pattern(aes(fill = anno_method, group = pops, pattern = pops), 
                     pattern_fill = "gray", pattern_color = "gray",
                     pattern_size = 0, pattern_density = 1/3,
                     stat = "identity", width = 0.75, 
                     position = position_dodge2(padding = 0.2), show.legend = T) +
    facet_grid(task2 ~ scenario, scales = "free",
               labeller = labeller(scenario = scenarios)) +
    theme_classic() + custom_theme() +
    scale_x_discrete(expand = expansion(mult = c(0.5,0.5)), 
                     breaks = N2_label, labels = c("20k" = TeX(
                       "$N^{(2)}=20k$"), "200k" = TeX("$N^{(2)}=200k$"))) +
    scale_pattern_manual(values = c("circle", "stripe")) +
    guides(fill = guide_legend(override.aes = list(pattern = "none"), nrow = 1)) +
    labs(x = NULL, y = TeX("Power improvement at CL or PIP \\geq 0.9"), 
         fill = "", title = NULL)
  
  p2 = ggplot(filter(power, LD == lds[k]), aes(
    x = factor(sprintf("%dk", N2/1e3), levels = N2_label), y = improve_CS95)) +
    geom_bar_pattern(aes(fill = anno_method, group = pops, pattern = pops), 
                     pattern_fill = "gray", pattern_color = "gray",
                     pattern_size = 0, pattern_density = 1/3,
                     stat = "identity", width = 0.75, 
                     position = position_dodge2(padding = 0.2), show.legend = T) +
    facet_grid(task2 ~ scenario, scales = "free",
               labeller = labeller(scenario = scenarios)) +
    theme_classic() + custom_theme() +
    scale_x_discrete(expand = expansion(mult = c(0.5,0.5)), 
                     breaks = N2_label, labels = c("20k" = TeX(
                       "$N^{(2)}=20k$"), "200k" = TeX("$N^{(2)}=200k$"))) +
    scale_pattern_manual(values = c("circle", "stripe")) +
    guides(fill = guide_legend(override.aes = list(pattern = "none"), nrow = 1)) +
    labs(x = NULL, y = TeX("Power improvement of 95% CSs"), 
         fill = "", title = NULL)
  
  print(ggarrange(p1, p2,
                  nrow = 2, ncol = 1, common.legend = T, legend = "bottom",
                  heights = c(1,1), labels = c("a","b"),
                  font.label = list(size = 12, color = "black", face = "bold", 
                                    family = NULL), align = "v"))
}
dev.off()
