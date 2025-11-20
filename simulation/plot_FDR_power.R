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

select_block = readRDS("../data/select_block.RData")

#####
# CL or PIP >= 0.9
data_0.9 = foreach(s = 1:3, .combine = "rbind") %do%
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
          n = calibration[[pop]][4,2]
          x = ifelse(n == 0, 0, 
                     round(calibration[[pop]][4,2] * calibration[[pop]][4,4]))
          fdr = ifelse(n == 0, 0, 1 - x / n)
          fdr_sd = sqrt(fdr * (1-fdr) / max(n,1))
          p = x / (n_causal[s,pops[pop]] * 200)
          p_sd = sqrt(p * (1-p) / (n_causal[s,pops[pop]] * 200))
          
          data.frame("LD" = ld,
                     "scenario" = s,
                     "N2" = N2_seq[m],
                     "POP" = pops[pop],
                     "method" = method,
                     "nVar_0.9" = n,
                     "FDR_0.9" = fdr,
                     "FDR_0.9_sd" = fdr_sd,
                     "power_0.9" = p,
                     "power_0.9_sd" = p_sd)
        }
      }
    }
  }
}

data_CS = foreach(s = 1:3, .combine = "rbind") %do%
{
  foreach(m = 1:2, .combine = "rbind") %do%
  {
    foreach(ld = lds, .combine = "rbind") %do%
    {
      foreach(method = names(ld_methods[[ld]]), .combine = "rbind") %do%
      {
        results = read.delim(sprintf(
          "setting_%s/%s/results_N2-%d.txt",
          s, ld_methods[[ld]][method], N2_seq[m]))
        
        CS_causal = c("Cross" = sum(results$CS95_causal),
                      "EUR" = sum(results$CS95_causal_1),
                      "EAS" = sum(results$CS95_causal_2),
                      "Shared" = sum(results$CS95_causal_0))
        
        results_CS95 = list()
        
        if(method %in% methods[1:7])
        {
          results_CS95[["Cross"]] = read.delim(sprintf(
            "setting_%s/%s/results_CS95_N2-%d_cross.txt",
            s, ld_methods[[ld]][method], N2_seq[m]))
          
          results_CS95[["EUR"]] = read.delim(sprintf(
            "setting_%s/%s/results_CS95_%d_N2-%d.txt",
            s, ld_methods[[ld]][method], 1, N2_seq[m]))
          
          results_CS95[["EAS"]] = read.delim(sprintf(
            "setting_%s/%s/results_CS95_%d_N2-%d.txt",
            s, ld_methods[[ld]][method], 2, N2_seq[m]))
          
          results_CS95[["Shared"]] = read.delim(sprintf(
            "setting_%s/%s/results_CS95_N2-%d_shared.txt",
            s, ld_methods[[ld]][method], N2_seq[m]))
        } else {
          results_CS95[["Cross"]] = read.delim(sprintf(
            "setting_%s/%s/results_CS95_N2-%d_cross.txt",
            s, ld_methods[[ld]][method], N2_seq[m]))
          
          results_CS95[["EUR"]] = read.delim(sprintf(
            "setting_%s/%s/results_CS95_N1-200000.txt",
            s, ld_methods[[ld]][method]))
          
          results_CS95[["EAS"]] = read.delim(sprintf(
            "setting_%s/%s/results_CS95_N2-%d.txt",
            s, ld_methods[[ld]][method], N2_seq[m]))
          
          results_CS95[["Shared"]] = read.delim(sprintf(
            "setting_%s/%s/results_CS95_N2-%d_shared.txt",
            s, ld_methods[[ld]][method], N2_seq[m]))
        }
        
        foreach(pop = pops, .combine = "rbind") %do%
        {
          n = nrow(results_CS95[[pop]])
          c = sum(results_CS95[[pop]]$coverage > 0) / n
          p = CS_causal[pop] / (n_causal[s,pop] * 200)
          
          n_CS = rep(0,200)
          names(n_CS) = select_block
          n_CS_table = table(results_CS95[[pop]]$block)
          n_CS[names(n_CS_table)] = n_CS_table
          
          data.frame("LD" = ld,
                     "scenario" = s,
                     "N2" = N2_seq[m],
                     "POP" = pop,
                     "method" = method,
                     "n_CS95" = n,
                     "n_CS95_sd" = sd(n_CS) / sqrt(200),
                     "coverage_CS95" = c,
                     "coverage_CS95_sd" = sqrt(c * (1-c) / n),
                     "power_CS95" = p,
                     "power_CS95_sd" = sqrt(p * (1-p) / (n_causal[s,pop] * 200)),
                     "size_CS95" = mean(results_CS95[[pop]]$size),
                     "size_sd" = sd(results_CS95[[pop]]$size) / sqrt(n),
                     "purity_CS95" = mean(results_CS95[[pop]]$min.abs.corr),
                     "purity_sd" = sd(results_CS95[[pop]]$min.abs.corr) / sqrt(n))
        }
      }
    }
  }
}

data = cbind(data_0.9, data_CS[,6:15])

data$LD %<>% factor(levels = lds)
data$POP %<>% factor(levels = pops)
data$method %<>% factor(levels = methods)
data %<>% arrange(LD, scenario, N2, POP, method)

write_delim(data, "FDR_power.txt", delim = '\t')
saveRDS(data, "FDR_power.RData")

data = readRDS("FDR_power.RData")

custom_theme = function()
{
  theme(
    aspect.ratio = 1/5,
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

# UKBB LD

ylimits = foreach(pop = pops, .combine = "rbind") %do%
{
  sub = filter(data, LD == lds[1] & POP == pop)
  
  data.frame("POP" = pop,
             "FDR_l" = floor(min(sub$FDR_0.9 - sub$FDR_0.9_sd, 
                                 na.rm = T) * 20) / 20,
             "FDR_u" = ceiling(max(sub$FDR_0.9 + sub$FDR_0.9_sd, 
                                   na.rm = T) * 20) / 20,
             "power_l" = floor(min(sub$power_0.9 - sub$power_0.9_sd, 
                                   na.rm = T) * 20) / 20,
             "power_u" = ceiling(max(sub$power_0.9 + sub$power_0.9_sd, 
                                     na.rm = T) * 20) / 20,
             "coverage_l" = floor(min(sub$coverage_CS95 - sub$coverage_CS95_sd, 
                                      na.rm = T) * 20) / 20,
             "coverage_u" = ceiling(max(sub$coverage_CS95 + sub$coverage_CS95_sd, 
                                        na.rm = T) * 20) / 20,
             "powerCS_l" = floor(min(sub$power_CS95 - sub$power_CS95_sd, 
                                     na.rm = T) * 20) / 20,
             "powerCS_u" = ceiling(max(sub$power_CS95 + sub$power_CS95_sd, 
                                       na.rm = T) * 20) / 20)
}

p_FDR_0.9 = ggplot(filter(data, LD == lds[1]), aes(
  x = factor(sprintf("%dk", N2/1e3), levels = N2_label), y = FDR_0.9)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8, 
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = FDR_0.9 - FDR_0.9_sd, ymax = FDR_0.9 + FDR_0.9_sd),
                width = 0.8, position = position_dodge2(), 
                linewidth = 0.25, color = "black") +
  facet_grid(POP ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(expand = expansion(mult = c(0.5,0.5)), 
                   breaks = N2_label, labels = c("20k" = TeX(
                     "$N^{(2)}=20k$"), "200k" = TeX("$N^{(2)}=200k$"))) +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(ylimits[1,2], ylimits[1,3])),
    scale_y_continuous(limits = c(ylimits[2,2], ylimits[2,3])),
    scale_y_continuous(limits = c(ylimits[3,2], ylimits[3,3])),
    scale_y_continuous(limits = c(ylimits[4,2], ylimits[4,3])))) +
  scale_fill_manual(values = hue_pal()(13)[c(1:11,13)]) +
  guides(fill = guide_legend(nrow = 2)) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey", 
             linewidth = 0.25) +
  labs(x = NULL, y = TeX("FDR at CL or PIP \\geq 0.9"), fill = "")

p_power_0.9 = ggplot(filter(data, LD == lds[1]), aes(
  x = factor(sprintf("%dk", N2/1e3), levels = c("20k","200k")), y = power_0.9)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8, 
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = power_0.9 - power_0.9_sd, 
                    ymax = power_0.9 + power_0.9_sd),
                width = 0.8, position = position_dodge2(), 
                linewidth = 0.25, color = "black") +
  facet_grid(POP ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(expand = expansion(mult = c(0.5,0.5)), 
                   breaks = N2_label, labels = c("20k" = TeX(
                     "$N^{(2)}=20k$"), "200k" = TeX("$N^{(2)}=200k$"))) +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(ylimits[1,4], ylimits[1,5])),
    scale_y_continuous(limits = c(ylimits[2,4], ylimits[2,5])),
    scale_y_continuous(limits = c(ylimits[3,4], ylimits[3,5])),
    scale_y_continuous(limits = c(ylimits[4,4], ylimits[4,5])))) +
  scale_fill_manual(values = hue_pal()(13)[c(1:11,13)]) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = NULL, y = TeX("Power at CL or PIP \\geq 0.9"), fill = "")

p_coverage_CS = ggplot(filter(data, LD == lds[1]), aes(
  x = factor(sprintf("%dk", N2/1e3), levels = N2_label), y = coverage_CS95)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8, 
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = coverage_CS95 - coverage_CS95_sd, 
                    ymax = coverage_CS95 + coverage_CS95_sd),
                width = 0.8, position = position_dodge2(), 
                linewidth = 0.25, color = "black") +
  facet_grid(POP ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(expand = expansion(mult = c(0.5,0.5)), 
                   breaks = N2_label, labels = c("20k" = TeX(
                     "$N^{(2)}=20k$"), "200k" = TeX("$N^{(2)}=200k$"))) +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(ylimits[1,6], ylimits[1,7])),
    scale_y_continuous(limits = c(ylimits[2,6], ylimits[2,7])),
    scale_y_continuous(limits = c(ylimits[3,6], ylimits[3,7])),
    scale_y_continuous(limits = c(ylimits[4,6], ylimits[4,7])))) +
  scale_fill_manual(values = hue_pal()(13)[c(1:11,13)]) +
  guides(fill = guide_legend(nrow = 2)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey", 
             linewidth = 0.25) +
  labs(x = NULL, y = "Coverage of 95% CSs", fill = "")

p_power_CS = ggplot(filter(data, LD == lds[1]), aes(
  x = factor(sprintf("%dk", N2/1e3), levels = c("20k","200k")), y = power_CS95)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8, 
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = power_CS95 - power_CS95_sd, 
                    ymax = power_CS95 + power_CS95_sd),
                width = 0.8, position = position_dodge2(), 
                linewidth = 0.25, color = "black") +
  facet_grid(POP ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(expand = expansion(mult = c(0.5,0.5)), 
                   breaks = N2_label, labels = c("20k" = TeX(
                     "$N^{(2)}=20k$"), "200k" = TeX("$N^{(2)}=200k$"))) +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(ylimits[1,8], ylimits[1,9])),
    scale_y_continuous(limits = c(ylimits[2,8], ylimits[2,9])),
    scale_y_continuous(limits = c(ylimits[3,8], ylimits[3,9])),
    scale_y_continuous(limits = c(ylimits[4,8], ylimits[4,9])))) +
  scale_fill_manual(values = hue_pal()(13)[c(1:11,13)]) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = NULL, y = "Power of 95% CSs", fill = "")

ggsave("FDR_power_UKBB.pdf",
       ggarrange(p_FDR_0.9, p_power_0.9, p_coverage_CS, p_power_CS, 
                 nrow = 4, ncol = 1, common.legend = T, legend = "bottom",
                 heights = c(1,1,1,1), labels = c("a","b","c","d"),
                 font.label = list(size = 12, color = "black", face = "bold", 
                                   family = NULL), align = "v"),
       device = "pdf", width = 9, height = 12.31, units = "in", bg = "white")

# 1kg LD
ylimits = foreach(pop = pops, .combine = "rbind") %do%
{
  sub = filter(data, LD == lds[2] & POP == pop)
  
  data.frame("POP" = pop,
             "FDR_l" = floor(min(sub$FDR_0.9 - sub$FDR_0.9_sd, 
                                 na.rm = T) * 20) / 20,
             "FDR_u" = ceiling(max(sub$FDR_0.9 + sub$FDR_0.9_sd, 
                                   na.rm = T) * 20) / 20,
             "power_l" = floor(min(sub$power_0.9 - sub$power_0.9_sd, 
                                   na.rm = T) * 20) / 20,
             "power_u" = ceiling(max(sub$power_0.9 + sub$power_0.9_sd, 
                                     na.rm = T) * 20) / 20,
             "coverage_l" = floor(min(sub$coverage_CS95 - sub$coverage_CS95_sd, 
                                      na.rm = T) * 20) / 20,
             "coverage_u" = ceiling(max(sub$coverage_CS95 + sub$coverage_CS95_sd, 
                                        na.rm = T) * 20) / 20,
             "powerCS_l" = floor(min(sub$power_CS95 - sub$power_CS95_sd, 
                                     na.rm = T) * 20) / 20,
             "powerCS_u" = ceiling(max(sub$power_CS95 + sub$power_CS95_sd, 
                                       na.rm = T) * 20) / 20)
}

p_FDR_0.9 = ggplot(filter(data, LD == lds[2]), aes(
  x = factor(sprintf("%dk", N2/1e3), levels = N2_label), y = FDR_0.9)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8, 
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = FDR_0.9 - FDR_0.9_sd, ymax = FDR_0.9 + FDR_0.9_sd),
                width = 0.8, position = position_dodge2(), 
                linewidth = 0.25, color = "black") +
  facet_grid(POP ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(expand = expansion(mult = c(0.5,0.5)), 
                   breaks = N2_label, labels = c("20k" = TeX(
                     "$N^{(2)}=20k$"), "200k" = TeX("$N^{(2)}=200k$"))) +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(ylimits[1,2], ylimits[1,3])),
    scale_y_continuous(limits = c(ylimits[2,2], ylimits[2,3])),
    scale_y_continuous(limits = c(ylimits[3,2], ylimits[3,3])),
    scale_y_continuous(limits = c(ylimits[4,2], ylimits[4,3])))) +
  scale_fill_manual(values = hue_pal()(13)[c(1:3,8:10,12,13)]) +
  guides(fill = guide_legend(nrow = 1)) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey", 
             linewidth = 0.25) +
  labs(x = NULL, y = TeX("FDR at CL or PIP \\geq 0.9"), fill = "")

p_power_0.9 = ggplot(filter(data, LD == lds[2]), aes(
  x = factor(sprintf("%dk", N2/1e3), levels = c("20k","200k")), y = power_0.9)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8, 
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = power_0.9 - power_0.9_sd, 
                    ymax = power_0.9 + power_0.9_sd),
                width = 0.8, position = position_dodge2(), 
                linewidth = 0.25, color = "black") +
  facet_grid(POP ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(expand = expansion(mult = c(0.5,0.5)), 
                   breaks = N2_label, labels = c("20k" = TeX(
                     "$N^{(2)}=20k$"), "200k" = TeX("$N^{(2)}=200k$"))) +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(ylimits[1,4], ylimits[1,5])),
    scale_y_continuous(limits = c(ylimits[2,4], ylimits[2,5])),
    scale_y_continuous(limits = c(ylimits[3,4], ylimits[3,5])),
    scale_y_continuous(limits = c(ylimits[4,4], ylimits[4,5])))) +
  scale_fill_manual(values = hue_pal()(13)[c(1:3,8:10,12,13)]) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = TeX("Power at CL or PIP \\geq 0.9"), fill = "")

p_coverage_CS = ggplot(filter(data, LD == lds[2]), aes(
  x = factor(sprintf("%dk", N2/1e3), levels = N2_label), y = coverage_CS95)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8, 
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = coverage_CS95 - coverage_CS95_sd, 
                    ymax = coverage_CS95 + coverage_CS95_sd),
                width = 0.8, position = position_dodge2(), 
                linewidth = 0.25, color = "black") +
  facet_grid(POP ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(expand = expansion(mult = c(0.5,0.5)), 
                   breaks = N2_label, labels = c("20k" = TeX(
                     "$N^{(2)}=20k$"), "200k" = TeX("$N^{(2)}=200k$"))) +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(ylimits[1,6], ylimits[1,7])),
    scale_y_continuous(limits = c(ylimits[2,6], ylimits[2,7])),
    scale_y_continuous(limits = c(ylimits[3,6], ylimits[3,7])),
    scale_y_continuous(limits = c(ylimits[4,6], ylimits[4,7])))) +
  scale_fill_manual(values = hue_pal()(13)[c(1:3,8:10,12,13)]) +
  guides(fill = guide_legend(nrow = 1)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey", 
             linewidth = 0.25) +
  labs(x = NULL, y = "Coverage of 95% CSs", fill = "")

p_power_CS = ggplot(filter(data, LD == lds[2]), aes(
  x = factor(sprintf("%dk", N2/1e3), levels = c("20k","200k")), y = power_CS95)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8, 
           position = position_dodge2()) +
  geom_errorbar(aes(ymin = power_CS95 - power_CS95_sd, 
                    ymax = power_CS95 + power_CS95_sd),
                width = 0.8, position = position_dodge2(), 
                linewidth = 0.25, color = "black") +
  facet_grid(POP ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(expand = expansion(mult = c(0.5,0.5)), 
                   breaks = N2_label, labels = c("20k" = TeX(
                     "$N^{(2)}=20k$"), "200k" = TeX("$N^{(2)}=200k$"))) +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(ylimits[1,8], ylimits[1,9])),
    scale_y_continuous(limits = c(ylimits[2,8], ylimits[2,9])),
    scale_y_continuous(limits = c(ylimits[3,8], ylimits[3,9])),
    scale_y_continuous(limits = c(ylimits[4,8], ylimits[4,9])))) +
  scale_fill_manual(values = hue_pal()(13)[c(1:3,8:10,12,13)]) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Power of 95% CSs", fill = "")

ggsave("FDR_power_1kG.pdf",
       ggarrange(p_FDR_0.9, p_power_0.9, p_coverage_CS, p_power_CS, 
                 nrow = 4, ncol = 1, common.legend = T, legend = "bottom",
                 heights = c(1,1,1,1), labels = c("a","b","c","d"),
                 font.label = list(size = 12, color = "black", face = "bold", 
                                   family = NULL), align = "v"),
       device = "pdf", width = 9, height = 12, units = "in", bg = "white")
