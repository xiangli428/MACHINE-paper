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

size_range = list("1" = c(1,1),
                  "2" = c(2,2),
                  "3-5" = c(3,5),
                  "6-10" = c(6,10),
                  ">10" = c(11,Inf))

data = foreach(s = 1:3, .combine = "rbind") %do%
{
  foreach(k = 2:3, .combine = "rbind") %do%
  {
    foreach(m = 1:2, .combine = "rbind") %do%
    {
      foreach(ld = lds, .combine = "rbind") %do%
      {
        foreach(method = names(ld_methods[[ld]]), .combine = "rbind") %do%
        {
          results_CS95 = list()
          
          results_CS95[["cross"]] = read.delim(sprintf(
            "setting_%s/%s/N%d-%d/results_CS95_cross.txt",
            s, ld_methods[[ld]][method], k, N2_seq[m]))
          
          results_CS95[["shared"]] = read.delim(sprintf(
            "setting_%s/%s/N%d-%d/results_CS95_shared.txt",
            s, ld_methods[[ld]][method], k, N2_seq[m]))
          
          if(method %in% methods[1:7])
          {
            results_CS95[["pop1"]] = read.delim(sprintf(
              "setting_%s/%s/N%d-%d/results_CS95_1.txt",
              s, ld_methods[[ld]][method], k, N2_seq[m]))
            
            results_CS95[["pop2"]] = read.delim(sprintf(
              "setting_%s/%s/N%d-%d/results_CS95_2.txt",
              s, ld_methods[[ld]][method], k, N2_seq[m]))
          } else {
            results_CS95[["pop1"]] = read.delim(sprintf(
              "setting_%s/%s/N1-200000/results_CS95.txt",
              s, ld_methods[[ld]][method]))
            
            results_CS95[["pop2"]] = read.delim(sprintf(
              "setting_%s/%s/N%d-%d/results_CS95.txt",
              s, ld_methods[[ld]][method], k, N2_seq[m]))
          }
          
          foreach(pid = names(tasks[[k-1]]), .combine = "rbind") %do%
          {
            foreach(g = names(size_range), .combine = "rbind") %do%
            {
              idx = which(results_CS95[[pid]]$size >= size_range[[g]][1] &
                            results_CS95[[pid]]$size <= size_range[[g]][2])
              n = length(idx)
              if(n > 0)
              {
                c = mean(results_CS95[[pid]]$coverage[idx] > 0)
                data.frame("LD" = ld,
                           "scenario" = s,
                           "pops" = sprintf("EUR+%s", pids[k]),
                           "N2" = N2_seq[m],
                           "task" = tasks[[k-1]][pid],
                           "method" = method,
                           "size" = g,
                           "n_CS95" = n,
                           "coverage" = c,
                           "coverage_sd" = sqrt(c * (1-c) / n))
              } else {
                data.frame("LD" = ld,
                           "scenario" = s,
                           "pops" = sprintf("EUR+%s", pids[k]),
                           "N2" = N2_seq[m],
                           "task" = tasks[[k-1]][pid],
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
}

data$LD %<>% factor(levels = lds)
data$pops %<>% factor(levels = c("EUR+AFR","EUR+EAS"))
data$method %<>% factor(levels = methods)
data$task %<>% factor(levels = c("Cross", "EUR", "AFR", "EAS", "Shared"))
data$size %<>% factor(levels = names(size_range))
data %<>% arrange(LD, scenario, pops, N2, task, method, size)

saveRDS(data, "coverage-size.RData")

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

colors = c("#521A13","#996035FF","#DA2222",
           "#FF9933FF","#FFDD00FF","#00AD00FF","#5D7A2BFF",
           "#009393FF","#69D2E7FF","#0066CCFF",
           "#1E2085FF","#8785B2FF","#953272FF")

data$task2 = as.character(data$task)
data$task2[data$task %in% c("AFR","EAS")] = "AFR / EAS"
data$task2 %<>% factor(levels = c("Cross","EUR","AFR / EAS","Shared"))

data$N2 = sprintf("N^{(2)}==%d*k", data$N2 / 1e3)
data$N2 %<>% factor(levels = sprintf("N^{(2)}==%d*k", N2_seq / 1e3))

pdf("coverage-size_UKBB.pdf", width = 9, height = 8.4, onefile = T, bg = "white")
for(s in 1:3)
{
  ylimits = filter(data, LD == lds[1] & scenario == s) %>%
    group_by(task2, N2) %>% summarise(
      coverage_l = floor(min(coverage - coverage_sd, na.rm = T) * 20) / 20,
      coverage_u = ceiling(max(coverage + coverage_sd, na.rm = T) * 20) / 20)
  
  print(ggplot(filter(data, LD == lds[1] & scenario == s),
               aes(x = size, y = coverage)) +
          geom_bar(aes(fill = method), stat = "identity", width = 0.8, 
                   position = position_dodge2()) +
          geom_errorbar(aes(ymin = coverage - coverage_sd, ymax = coverage + coverage_sd),
                        width = 0.8, position = position_dodge2(), 
                        linewidth = 0.25, color = "black") +
          facet_grid(task2 + N2 ~ pops, scales = "free", labeller = label_parsed) +
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
          scale_fill_manual(values = colors[c(1:11,13)]) +
          guides(fill = guide_legend(nrow = 2)) +
          geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey", 
                     linewidth = 0.25) +
          geom_vline(xintercept = 1.5, color = "grey", linewidth = 0.25) +
          geom_vline(xintercept = 2.5, color = "grey", linewidth = 0.25) +
          geom_vline(xintercept = 3.5, color = "grey", linewidth = 0.25) +
          geom_vline(xintercept = 4.5, color = "grey", linewidth = 0.25) +
          labs(x = "Size of 95% CSs", y = "Coverage of 95% CSs", fill = NULL, 
               title = scenarios[s]))
}
dev.off()

pdf("coverage-size_1kG.pdf", width = 9, height = 8.1, onefile = T, bg = "white")
for(s in 1:3)
{
  ylimits = filter(data, LD == lds[2] & scenario == s) %>%
    group_by(task2, N2) %>% summarise(
      coverage_l = floor(min(coverage - coverage_sd, na.rm = T) * 20) / 20,
      coverage_u = ceiling(max(coverage + coverage_sd, na.rm = T) * 20) / 20)
  
  print(ggplot(filter(data, LD == lds[2] & scenario == s),
               aes(x = size, y = coverage)) +
          geom_bar(aes(fill = method), stat = "identity", width = 0.8, 
                   position = position_dodge2()) +
          geom_errorbar(aes(ymin = coverage - coverage_sd, ymax = coverage + coverage_sd),
                        width = 0.8, position = position_dodge2(), 
                        linewidth = 0.25, color = "black") +
          facet_grid(task2 + N2 ~ pops, scales = "free", labeller = label_parsed) +
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
          scale_fill_manual(values = colors[c(1:3,8:10,12,13)]) +
          guides(fill = guide_legend(nrow = 1)) +
          geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey", 
                     linewidth = 0.25) +
          geom_vline(xintercept = 1.5, color = "grey", linewidth = 0.25) +
          geom_vline(xintercept = 2.5, color = "grey", linewidth = 0.25) +
          geom_vline(xintercept = 3.5, color = "grey", linewidth = 0.25) +
          geom_vline(xintercept = 4.5, color = "grey", linewidth = 0.25) +
          labs(x = "Size of 95% CSs", y = "Coverage of 95% CSs", fill = NULL, 
               title = scenarios[s]))
}
dev.off()
