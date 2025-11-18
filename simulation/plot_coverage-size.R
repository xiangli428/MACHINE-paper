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
library(Hmisc)

setwd("simulation/results")

source("../simulation_setting.R")

size_range = list("1" = c(1,1),
                  "2" = c(2,2),
                  "3-5" = c(3,5),
                  "6-10" = c(6,10),
                  ">10" = c(11,Inf))

data = foreach(s = 1:3, .combine = "rbind") %do%
{
  foreach(m = 1:2, .combine = "rbind") %do%
  {
    foreach(ld = lds, .combine = "rbind") %do%
    {
      foreach(method = names(ld_methods[[ld]]), .combine = "rbind") %do%
      {
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
          foreach(g = names(size_range), .combine = "rbind") %do%
          {
            idx = which(results_CS95[[pop]]$size >= size_range[[g]][1] &
                          results_CS95[[pop]]$size <= size_range[[g]][2])
            n = length(idx)
            if(n > 0)
            {
              c = mean(results_CS95[[pop]]$coverage[idx] > 0)
              data.frame("LD" = ld,
                         "scenario" = s,
                         "N2" = N2_seq[m],
                         "POP" = pop,
                         "method" = method,
                         "size" = g,
                         "n_CS95" = n,
                         "coverage" = c,
                         "coverage_sd" = sqrt(c * (1-c) / n))
            } else {
              data.frame("LD" = ld,
                         "scenario" = s,
                         "N2" = N2_seq[m],
                         "POP" = pop,
                         "method" = method,
                         "size" = g,
                         "n_CS95" = 0,
                         "coverage" = NaN,
                         "coverage_sd" = NaN)
            }
          }
        }
      }
    }
  }
}

data$LD %<>% factor(levels = lds)
data$POP %<>% factor(levels = pops)
data$method %<>% factor(levels = methods)
data$size %<>% factor(levels = names(size_range))
data %<>% arrange(LD, scenario, N2, POP, method, size)

write_delim(data, "coverage-size.txt", delim = '\t')
saveRDS(data, "coverage-size.RData")
data = readRDS("coverage-size.RData")

custom_theme = function()
{
  theme(
    aspect.ratio = 1/10,
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

data$N2 = sprintf("N^{(2)}==%d*k", data$N2 / 1e3)
data$N2 %<>% factor(levels = sprintf("N^{(2)}==%d*k", N2_seq / 1e3))

pdf("coverage-size_UKBB.pdf", width = 9, height = 8.2, onefile = T, bg = "white")
for(s in 1:3)
{
  ylimits = foreach(pop = pops, .combine = "rbind") %do%
  {
    foreach(m = 1:2, .combine = "rbind") %do%
    {
      sub = filter(data, LD == lds[1] & scenario == s & POP == pop & 
                     N2 == sprintf("N^{(2)}==%d*k", N2_seq[m] / 1e3))
      data.frame("POP" = pop,
                 "N2" = N2_seq[m],
                 "coverage_l" = floor(min(sub$coverage - sub$coverage_sd, 
                                          na.rm = T) * 20) / 20,
                 "coverage_u" = ceiling(max(sub$coverage + sub$coverage_sd, 
                                            na.rm = T) * 20) / 20)
      
    }
  }
  
  print(ggplot(filter(data, LD == lds[1] & scenario == s),
               aes(x = size, y = coverage)) +
          geom_bar(aes(fill = method), stat = "identity", width = 0.8, 
                   position = position_dodge2()) +
          geom_errorbar(aes(ymin = coverage - coverage_sd, 
                            ymax = coverage + coverage_sd),
                        width = 0.8, position = position_dodge2(), 
                        linewidth = 0.25, color = "black") +
          facet_grid(POP + N2 ~ ., scales = "free", labeller = label_parsed) +
          theme_classic() + custom_theme() + 
          scale_x_discrete(expand = expansion(mult = c(0.125,0.125))) +
          facetted_pos_scales(y = list(
            scale_y_continuous(limits = c(ylimits[1,3], ylimits[1,4])),
            scale_y_continuous(limits = c(ylimits[2,3], ylimits[2,4])),
            scale_y_continuous(limits = c(ylimits[3,3], ylimits[3,4])),
            scale_y_continuous(limits = c(ylimits[4,3], ylimits[4,4])),
            scale_y_continuous(limits = c(ylimits[5,3], ylimits[5,4])),
            scale_y_continuous(limits = c(ylimits[6,3], ylimits[6,4])),
            scale_y_continuous(limits = c(ylimits[7,3], ylimits[7,4])),
            scale_y_continuous(limits = c(ylimits[8,3], ylimits[8,4])))) +
          scale_fill_manual(values = hue_pal()(13)[c(1:11,13)]) +
          guides(fill = guide_legend(nrow = 2)) +
          geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey", 
                     linewidth = 0.25) +
          geom_vline(xintercept = 1.5, color = "grey", linewidth = 0.25) +
          geom_vline(xintercept = 2.5, color = "grey", linewidth = 0.25) +
          geom_vline(xintercept = 3.5, color = "grey", linewidth = 0.25) +
          geom_vline(xintercept = 4.5, color = "grey", linewidth = 0.25) +
          labs(x = "Size of 95% CSs", y = "Coverage of 95% CSs", fill = "", 
               title = scenarios[s]))
}
dev.off()

pdf("coverage-size_1kG.pdf", width = 9, height = 8.2, onefile = T, bg = "white")
for(s in 1:3)
{
  ylimits = foreach(pop = pops, .combine = "rbind") %do%
  {
    foreach(m = 1:2, .combine = "rbind") %do%
    {
      sub = filter(data, LD == lds[2] & scenario == s & POP == pop & 
                     N2 == sprintf("N^{(2)}==%d*k", N2_seq[m] / 1e3))
      data.frame("POP" = pop,
                 "N2" = N2_seq[m],
                 "coverage_l" = floor(min(sub$coverage - sub$coverage_sd, 
                                          na.rm = T) * 20) / 20,
                 "coverage_u" = ceiling(max(sub$coverage + sub$coverage_sd, 
                                            na.rm = T) * 20) / 20)
      
    }
  }
  
  print(ggplot(filter(data, LD == lds[2] & scenario == s),
               aes(x = size, y = coverage)) +
          geom_bar(aes(fill = method), stat = "identity", width = 0.8, 
                   position = position_dodge2()) +
          geom_errorbar(aes(ymin = coverage - coverage_sd, 
                            ymax = coverage + coverage_sd),
                        width = 0.8, position = position_dodge2(), 
                        linewidth = 0.25, color = "black") +
          facet_grid(POP + N2 ~ ., scales = "free", labeller = label_parsed) +
          theme_classic() + custom_theme() + 
          scale_x_discrete(expand = expansion(mult = c(0.125,0.125))) +
          facetted_pos_scales(y = list(
            scale_y_continuous(limits = c(ylimits[1,3], ylimits[1,4])),
            scale_y_continuous(limits = c(ylimits[2,3], ylimits[2,4])),
            scale_y_continuous(limits = c(ylimits[3,3], ylimits[3,4])),
            scale_y_continuous(limits = c(ylimits[4,3], ylimits[4,4])),
            scale_y_continuous(limits = c(ylimits[5,3], ylimits[5,4])),
            scale_y_continuous(limits = c(ylimits[6,3], ylimits[6,4])),
            scale_y_continuous(limits = c(ylimits[7,3], ylimits[7,4])),
            scale_y_continuous(limits = c(ylimits[8,3], ylimits[8,4])))) +
          scale_fill_manual(values = hue_pal()(13)[c(1:3,8:10,12,13)]) +
          guides(fill = guide_legend(nrow = 2)) +
          geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey", 
                     linewidth = 0.25) +
          geom_vline(xintercept = 1.5, color = "grey", linewidth = 0.25) +
          geom_vline(xintercept = 2.5, color = "grey", linewidth = 0.25) +
          geom_vline(xintercept = 3.5, color = "grey", linewidth = 0.25) +
          geom_vline(xintercept = 4.5, color = "grey", linewidth = 0.25) +
          labs(x = "Size of 95% CSs", y = "Coverage of 95% CSs", fill = "", 
               title = scenarios[s]))
}
dev.off()
