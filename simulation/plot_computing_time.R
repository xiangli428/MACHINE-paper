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

select_block_info = read.delim("../data/select_block_info.txt")

methods = c("MACHINE", "MESuSiE", "XMAP", "SuSiEx", "MultiSuSiE", 
            "h2-D2", "SuSiE", "RSparsePro", "CARMA")

ld_methods = list(
  "In-sample LD" = c(
    "MACHINE" = "MACHINE",
    "MESuSiE" = "MESuSiE",
    "SuSiEx" = "SuSiEx",
    "XMAP" = "XMAP",
    "MultiSuSiE" = "MultiSuSiE",
    "h2-D2" = "h2D2-polyfun",
    "SuSiE" = "SuSiE",
    "CARMA" = "CARMA"),
  "1kG LD" = c(
    "MACHINE" = "MACHINE_1kg",
    "h2-D2" = "h2D2_1kg",
    "RSparsePro" = "RSparsePro_1kg",
    "CARMA" = "CARMA_1kg"))

data = foreach(s = 1:3, .combine = "rbind") %do%
{
  foreach(m = 1:2, .combine = "rbind") %do%
  {
    foreach(ld = lds, .combine = "rbind") %do%
    {
      foreach(method = names(ld_methods[[ld]]), .combine = "rbind") %do%
      {
        if(method %in% methods[1:5])
        {
          results = read.delim(sprintf(
            "setting_%s/%s/results_N2-%d.txt",
            s, ld_methods[[ld]][method], N2_seq[m]))
        } else {
          results_1 = read.delim(sprintf(
            "setting_%s/%s/results_N1-200000.txt",
            s, ld_methods[[ld]][method]))
          results_2 = read.delim(sprintf(
            "setting_%s/%s/results_N2-%d.txt",
            s, ld_methods[[ld]][method], N2_seq[m]))
        }
        
        data.frame("LD" = ld,
                   "scenario" = s,
                   "N2" = N2_seq[m],
                   "method" = method,
                   "block" = select_block_info$block,
                   "num_variants" = select_block_info$num_variants,
                   "time" = time)
      }
    }
  }
}

data$LD %<>% factor(levels = lds)
data$method %<>% factor(levels = methods)

saveRDS(data, "computing_time.RData")
data = readRDS("computing_time.RData")

data_mean = group_by(data, LD, method, block, num_variants) %>% 
  summarise(time_mean = mean(time))

model_info = foreach(ld = lds, .combine = "rbind") %do%
{
  foreach(med = names(ld_methods[[ld]]), .combine = "rbind") %do%
  {
    data_sub = filter(data_mean, LD == ld & method == med)
    mod = summary(lm(log10(time_mean) ~ log10(num_variants), data = data_sub))
    data.frame("LD" = ld,
               "method" = med,
               "intercept" = mod$coefficients[1,1],
               "intercept_sd" = mod$coefficients[1,2],
               "slope" = mod$coefficients[2,1],
               "slope_sd" = mod$coefficients[2,2],
               "adj.r.squared" = mod$adj.r.squared)
  }
}

model_info$LD %<>% factor(levels = lds)
model_info$method %<>% factor(levels = methods)

write_delim(model_info, "computing_time_model_info.txt", delim = '\t')
saveRDS(model_info, "computing_time_model_info.RData")
model_info = readRDS("computing_time_model_info.RData")

custom_theme = function()
{
  theme(
    aspect.ratio = 1,
    panel.spacing = unit(0.05, "in"),
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
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA))
}

model_info %<>% mutate(label = if_else(intercept > 0, sprintf(
  "y==%.2f*x+%.2f~~R^2==%.2f", slope, intercept, adj.r.squared), sprintf(
    "y==%.2f*x%.2f~~R^2==%.2f", slope, intercept, adj.r.squared)))
model_info$x = 3.05
model_info$y = 0.5
model_info$y[is.element(model_info$method, c("MESuSiE", "XMAP", "MultiSuSiE", "SuSiE"))] = 3

p1 = ggplot(filter(data_mean, LD == lds[1] & is.element(method, names(
  ld_methods[[1]])[1:4])), aes(x = log10(num_variants), y = log10(time_mean))) +
  geom_point(aes(color = method), size = 0.5) +
  geom_text(data = filter(model_info, LD == lds[1] & is.element(method, names(
    ld_methods[[1]])[1:4])), aes(x = x, y = y, label = label), size = 2.5, 
    color = "black", parse = T) +
  facet_grid(. ~ method, scales = "free") +
  theme_classic() + custom_theme() +
  scale_y_continuous(limits = c(-0.5,4)) +
  scale_color_manual(values = hue_pal()(13)[c(3:6)]) +
  guides(color = "none") +
  labs(x = NULL, y = NULL)

p2 = ggplot(filter(data_mean, LD == lds[1] & is.element(method, names(
  ld_methods[[1]])[5:8])), aes(x = log10(num_variants), y = log10(time_mean))) +
  geom_point(aes(color = method), size = 0.5) +
  geom_text(data = filter(model_info, LD == lds[1] & is.element(method, names(
    ld_methods[[1]])[5:8])), aes(x = x, y = y, label = label), size = 2.5, 
    color = "black", parse = T) +
  facet_grid(. ~ method, scales = "free") +
  theme_classic() + custom_theme() +
  scale_y_continuous(limits = c(-0.5,4)) +
  scale_color_manual(values = hue_pal()(13)[c(7,10,11,13)]) +
  guides(color = "none") +
  labs(x = NULL, y = NULL)

p3 = ggplot(filter(data_mean, LD == lds[2] & is.element(method, names(
  ld_methods[[2]]))), aes(x = log10(num_variants), y = log10(time_mean))) +
  geom_point(aes(color = method), size = 0.5) +
  geom_text(data = filter(model_info, LD == lds[2] & is.element(method, names(
    ld_methods[[2]]))), aes(x = x, y = y, label = label), size = 2.5, 
    color = "black", parse = T) +
  facet_grid(. ~ method, scales = "free") +
  theme_classic() + custom_theme() +
  scale_y_continuous(limits = c(-0.5,4)) +
  scale_color_manual(values = hue_pal()(13)[c(3,10,12,13)]) +
  guides(color = "none") +
  labs(x = NULL, y = NULL)

p12 = annotate_figure(ggarrange(
  p1, p2, nrow = 2, ncol = 1, heights = c(1,1), align = "hv"),
  top = text_grob(lds[1], face = "bold", size = 12))

p3 = annotate_figure(p3, top = text_grob(lds[2], face = "bold", size = 12))

p123 = ggarrange(p12, p3, nrow = 2, ncol = 1, heights = c(2,1.05), align = "hv")

p123 = annotate_figure(p123, left = text_grob(TeX(
  "$log_{10}$ (CPU time/s)"), rot = 90, size = 8),
  bottom = text_grob(TeX("$log_{10}$(Num variants)"), size = 8))

ggsave("computing_time.pdf", p123,
       device = "pdf", width = 9, height = 8.6, units = "in", bg = "white")
