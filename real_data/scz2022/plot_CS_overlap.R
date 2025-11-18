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

data = foreach(method = names(methods), .combine = "rbind") %dopar%
{
  CS_overlap = read.delim(sprintf("results/%s/CS_overlap.txt", methods[method]))
  
  if(method %in% c("MACHINE + g-LDSC", "MACHINE", "MESuSiE"))
  {
    CS_info = read.delim(sprintf("results/%s/CS_info.txt", methods[method]))
  } else {
    CS_info = foreach(pid = pids, .combine = "rbind") %do%
    {
      read.delim(sprintf("results/%s/%s/CS_info.txt", methods[method], pid))
    }
  }
  
  data.frame("method" = method,
             "EUR_EAS" = nrow(CS_overlap),
             "EUR-specific" = length(setdiff(grep("EUR", CS_info$CS, value = T),
                                             CS_overlap$CS.EUR)), 
             "EAS-specific" = length(setdiff(grep("EAS", CS_info$CS, value = T),
                                             CS_overlap$CS.EAS)),
             check.names = F)
}

data$method %<>% factor(levels = names(methods))

saveRDS(data, "CS_overlap_num.RData")

data = readRDS("CS_overlap_num.RData")

data %<>% gather(pops, nCSs, -method)
data$pops %<>% gsub("_", " & ", .)
data$pops %<>% factor(levels = c("EUR & EAS", "EUR-specific", "EAS-specific"))

custom_theme = function()
{
  theme(
    aspect.ratio = 3/5,
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
    panel.spacing = unit(0.05, "in"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA))
}

p = ggplot(data, aes(x = method, y = nCSs)) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.8,
           position = position_dodge2()) +
  geom_text(aes(label = nCSs), size = 2, vjust = -0.3) +
  facet_grid(pops ~ ., scales = "free_y") +
  theme_classic() + custom_theme() +
  scale_y_continuous(expand = expansion(mult = c(0,0.15))) +
  scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
  guides(fill = "none") +
  labs(x = NULL, y = "Number of CSs", fill = "", title = NULL)

saveRDS(p, "plot_CS_overlap_num.RData")
