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
library(ggh4x)
library(latex2exp)
library(scales)
library(Hmisc)

setwd("simulation/results")

source("../simulation_setting.R")

data = readRDS("FDR_power.RData")

custom_theme = function()
{
  theme(
    aspect.ratio = 1/5,
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

ylimits = foreach(pop = pops, .combine = "rbind") %do%
{
  sub = filter(data, LD == lds[1] & POP == pop)
  
  data.frame("POP" = pop,
             "purity_l" = floor(min(sub$purity_CS95 - sub$purity_CS95_sd, 
                                    na.rm = T) * 20) / 20)
}

p_n_CS = ggplot(filter(data, LD == lds[1]), aes(
  x = factor(sprintf("%dk", N2/1e3), levels = c("20k","200k")), y = n_CS95 / 200)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8, 
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = n_CS95 / 200 - n_CS95_sd, ymax = n_CS95 / 200 + n_CS95_sd),
                width = 0.8, position = position_dodge2(), 
                linewidth = 0.25, color = "black") +
  facet_grid(POP ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(expand = expansion(mult = c(0.5,0.5)), 
                   breaks = N2_label, labels = c("20k" = TeX(
                     "$N^{(2)}=20k$"), "200k" = TeX("$N^{(2)}=200k$"))) +
  scale_fill_manual(values = hue_pal()(13)[c(1:11,13)]) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = NULL, y = "Number of 95% CSs", fill = "")

p_size_CS = ggplot(filter(data, LD == lds[1]), aes(
  x = factor(sprintf("%dk", N2/1e3), levels = c("20k","200k")), y = size_CS95)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8, 
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = size_CS95 - size_CS95_sd, ymax = size_CS95 + size_CS95_sd),
                width = 0.8, position = position_dodge2(), 
                linewidth = 0.25, color = "black") +
  facet_grid(POP ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(expand = expansion(mult = c(0.5,0.5)), 
                   breaks = N2_label, labels = c("20k" = TeX(
                     "$N^{(2)}=20k$"), "200k" = TeX("$N^{(2)}=200k$"))) +
  scale_fill_manual(values = hue_pal()(13)[c(1:11,13)]) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = NULL, y = "Size of 95% CSs", fill = "")

p_purity_CS = ggplot(filter(data, LD == lds[1]), aes(
  x = factor(sprintf("%dk", N2/1e3), levels = c("20k","200k")), y = purity_CS95)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8, 
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = purity_CS95 - purity_CS95_sd, ymax = purity_CS95 + purity_CS95_sd),
                width = 0.8, position = position_dodge2(), 
                linewidth = 0.25, color = "black") +
  facet_grid(POP ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(expand = expansion(mult = c(0.5,0.5)), 
                   breaks = N2_label, labels = c("20k" = TeX(
                     "$N^{(2)}=20k$"), "200k" = TeX("$N^{(2)}=200k$"))) +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(ylimits[1,2], 1)),
    scale_y_continuous(limits = c(ylimits[2,2], 1)),
    scale_y_continuous(limits = c(ylimits[3,2], 1)),
    scale_y_continuous(limits = c(ylimits[4,2], 1)))) +
  scale_fill_manual(values = hue_pal()(13)[c(1:11,13)]) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = NULL, y = "Purity of 95% CSs", fill = "")

ggsave("size_purity_UKBB.pdf",
       ggarrange(p_n_CS, p_size_CS, p_purity_CS,
                 nrow = 3, ncol = 1, common.legend = T, legend = "bottom",
                 heights = c(1,1,1), labels = c("a","b","c"),
                 font.label = list(size = 12, color = "black", face = "bold", 
                                   family = NULL), align = "v"),
       device = "pdf", width = 9, height = 9.4, units = "in", bg = "white")


ylimits = foreach(pop = pops, .combine = "rbind") %do%
{
  sub = filter(data, LD == lds[2] & POP == pop)
  
  data.frame("POP" = pop,
             "purity_l" = floor(min(sub$purity_CS95 - sub$purity_CS95_sd, 
                                    na.rm = T) * 20) / 20)
}

p_n_CS = ggplot(filter(data, LD == lds[2]), aes(
  x = factor(sprintf("%dk", N2/1e3), levels = c("20k","200k")), y = n_CS95 / 200)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8, 
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = n_CS95 / 200 - n_CS95_sd, ymax = n_CS95 / 200 + n_CS95_sd),
                width = 0.8, position = position_dodge2(), 
                linewidth = 0.25, color = "black") +
  facet_grid(POP ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(expand = expansion(mult = c(0.5,0.5)), 
                   breaks = N2_label, labels = c("20k" = TeX(
                     "$N^{(2)}=20k$"), "200k" = TeX("$N^{(2)}=200k$"))) +
  scale_fill_manual(values = hue_pal()(13)[c(1:3,8:10,12,13)]) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Number of 95% CSs", fill = "")

p_size_CS = ggplot(filter(data, LD == lds[2]), aes(
  x = factor(sprintf("%dk", N2/1e3), levels = c("20k","200k")), y = size_CS95)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8, 
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = size_CS95 - size_CS95_sd, ymax = size_CS95 + size_CS95_sd),
                width = 0.8, position = position_dodge2(), 
                linewidth = 0.25, color = "black") +
  facet_grid(POP ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(expand = expansion(mult = c(0.5,0.5)), 
                   breaks = N2_label, labels = c("20k" = TeX(
                     "$N^{(2)}=20k$"), "200k" = TeX("$N^{(2)}=200k$"))) +
  scale_fill_manual(values = hue_pal()(13)[c(1:3,8:10,12,13)]) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Size of 95% CSs", fill = "")

p_purity_CS = ggplot(filter(data, LD == lds[2]), aes(
  x = factor(sprintf("%dk", N2/1e3), levels = c("20k","200k")), y = purity_CS95)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8, 
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = purity_CS95 - purity_CS95_sd, ymax = purity_CS95 + purity_CS95_sd),
                width = 0.8, position = position_dodge2(), 
                linewidth = 0.25, color = "black") +
  facet_grid(POP ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(expand = expansion(mult = c(0.5,0.5)), 
                   breaks = N2_label, labels = c("20k" = TeX(
                     "$N^{(2)}=20k$"), "200k" = TeX("$N^{(2)}=200k$"))) +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(ylimits[1,2], 1)),
    scale_y_continuous(limits = c(ylimits[2,2], 1)),
    scale_y_continuous(limits = c(ylimits[3,2], 1)),
    scale_y_continuous(limits = c(ylimits[4,2], 1)))) +
  scale_fill_manual(values = hue_pal()(13)[c(1:3,8:10,12,13)]) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Purity of 95% CSs", fill = "")

ggsave("size_purity_1kG.pdf",
       ggarrange(p_n_CS, p_size_CS, p_purity_CS,
                 nrow = 3, ncol = 1, common.legend = T, legend = "bottom",
                 heights = c(1,1,1), labels = c("a","b","c"),
                 font.label = list(size = 12, color = "black", face = "bold", 
                                   family = NULL), align = "v"),
       device = "pdf", width = 9, height = 9.1, units = "in", bg = "white")
