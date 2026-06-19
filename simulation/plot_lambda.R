options(stringsAsFactors = F, check.names = F)

library(readr)
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

setwd("simulation/results")

select_block = readRDS("../data/select_block.RData")

source("../simulation_setting.R")

registerDoParallel(40)


data = foreach(s = 1:3, .combine = "rbind") %do%
{
  foreach(suffix = c("N1-200000", "N2-20000", "N2-200000", "N3-20000","N3-200000"), 
          .combine = "rbind") %do%
  {
    foreach(block = select_block, .combine = "rbind") %dopar%
    {
      vare = system(sprintf(
        "cat setting_%s/RSparsePro_1kg/%s/%s/result.rsparsepro.log %s",
        s, suffix, block, "|grep vare |tail -n 1 |cut -d ' ' -f 6"), intern = T)
      h2D2 = readRDS(sprintf("setting_%s/h2D2_1kg/%s/%s/h2D2.RData",
                             s, suffix, block))
      data.frame("scenario" = s,
                 "suffix" = suffix,
                 "block" = block,
                 "h2D2" = h2D2@lambda,
                 "RSparsePro" = as.numeric(vare))
    }
  }
}

saveRDS(data, "lambda.RData")

custom_theme1 = function()
{
  theme(
    axis.text = element_text(size = 7),  
    axis.title = element_text(size = 8),
    axis.ticks.x = element_line(linewidth = 0.25),
    axis.ticks.y = element_line(linewidth = 0.25),
    axis.line = element_blank(),
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.spacing = unit(0.05, "in"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA))
}

custom_theme2 = function()
{
  theme(
    axis.text.x = element_blank(),  
    axis.text.y = element_text(size = 7),  
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.ticks.x = element_line(linewidth = 0.25),
    axis.ticks.y = element_line(linewidth = 0.25),
    axis.line = element_blank(),
    strip.text.x = element_text(size = 9),
    strip.text.y = element_blank(),
    strip.background.x = element_rect(
      fill = "lightgray", color = "lightgray"),
    strip.background.y = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.spacing = unit(0.05, "in"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA))
}

custom_theme3 = function()
{
  theme(
    axis.text.x = element_text(size = 7),  
    axis.text.y = element_blank(),  
    axis.title.x = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.25),
    axis.ticks.y = element_line(linewidth = 0.25),
    axis.line = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_text(size = 9),
    strip.background.x = element_blank(),
    strip.background.y = element_rect(
      fill = "lightgray", color = "lightgray"),
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.spacing = unit(0.05, "in"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA))
}

suffixes = c("N1-200000" = "N^{(EUR)}==200*k",
             "N2-20000" = "N^{(AFR)}==20*k",
             "N2-200000" = "N^{(AFR)}==200*k",
             "N3-20000" = "N^{(EAS)}==20*k",
             "N3-200000" = "N^{(EAS)}==200*k")
data$suffix = factor(suffixes[data$suffix], levels = suffixes)

p1 = ggplot(data, aes(x = h2D2, y = RSparsePro)) + 
  geom_point(color = "tomato", size = 0.5, alpha = 0.5) +
  facet_grid(suffix ~ scenario, scales = "free") +
  theme_classic() + custom_theme1() + 
  scale_x_continuous(trans = log1p_trans(), limits = c(0,140),
                     breaks = c(0,2^(0:7))) +
  scale_y_continuous(trans = log1p_trans(), limits = c(0,140),
                     breaks = c(0,2^(0:7))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.25) + 
  labs(x = "MACHINE / h2-D2", y = "RSparsePro")

p2 = ggplot(data, aes(x = h2D2)) + 
  geom_histogram(fill = "tomato", color = "black", binwidth = 0.1, alpha = 0.5,
                 linewidth = 0.1) +
  facet_grid(suffix ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios)) +
  theme_classic() + custom_theme2() + 
  scale_x_continuous(trans = log1p_trans(), limits = c(0,140), oob = oob_keep,
                     breaks = c(0,2^(0:7))) +
  facetted_pos_scales(y = list(
    scale_y_continuous(breaks = c(0,10,20)),
    scale_y_continuous(breaks = c(0,40,80,120)),
    scale_y_continuous(breaks = c(0,10,20,30,40)),
    scale_y_continuous(breaks = c(0,40,80,120)),
    scale_y_continuous(breaks = c(0,10,20,30))))

p3 = ggplot(data, aes(x = RSparsePro)) + 
  geom_histogram(fill = "tomato", color = "black", binwidth = 0.1, alpha = 0.5,
                 linewidth = 0.1) +
  facet_grid(suffix ~ scenario, scales = "free",
             labeller = labeller(suffix = label_parsed)) +
  theme_classic() + custom_theme3() + coord_flip() +
  scale_x_continuous(trans = log1p_trans(), limits = c(0,140), oob = oob_keep,
                     breaks = c(0,2^(0:7))) +
  scale_y_continuous(breaks = c(0,100))

pdf("lambda.pdf", width = 9, height = 13.5, useDingbats = F, bg = "white")
p2 + plot_spacer() + p1 + p3 + 
  plot_layout(nrow = 2, ncol = 2, widths = c(4,1), heights = c(1,4))
dev.off()
