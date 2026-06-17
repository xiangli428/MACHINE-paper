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

setwd("real_data")


data = foreach(method = c("MACHINE", "MESuSiE"), .combine = "rbind") %do%
{
  results = read.delim(sprintf("scz2022/results/%s/results.txt", method))
  min_p = read.delim("scz2022/data/min_p_sub.txt")
  if(method == "MACHINE")
  {
    M = sqrt(rowSums(min_p[,4:5]^2) / 2)
  } else {
    M = results$num_variants
  }
  
  data.frame("method" = method,
             "pheno" = "SCZ",
             "chromosome" = results$chromosome,
             "block" = results$block,
             "num_variants" = results$time,
             "num_ancestries" = 2,
             "time" = time) %>% 
    rbind(
      foreach(pheno = c("HDL", "LDL", "TG", "TC"), .combine = "rbind") %do%
        {
          min_p = read.delim(sprintf("glgc/data/GLGC/%s/min_p_sub.txt", pheno))
          
          results = read.delim(sprintf("glgc/results/%s/%s/GLGC/results.txt", method, pheno))
          
          if(method == "MACHINE")
          {
            P = rowSums(min_p[,11:13] != 0)
            M = sqrt(rowSums(min_p[,11:13]^2) / P)
          } else {
            min_p = min_p[results$L == 5,]
            results %<>% filter(L == 5)
            P = results$num_ancestries
            M = results$num_variants
          }
          
          data.frame("method" = method,
                     "pheno" = pheno,
                     "chromosome" = results$chromosome,
                     "block" = results$block,
                     "num_variants" = M,
                     "num_ancestries" = P,
                     "time" = results$time)
        })
}

data$method %<>% factor(levels = c("MACHINE", "MESuSiE"))
saveRDS(data, "computing_time.RData")

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
    legend.title = element_text(size = 8),
    legend.position = "bottom",
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA))
}

p1 = ggplot(filter(data, method == "MACHINE"), 
            aes(x = log10(num_variants), y = log10(time))) +
  geom_point(aes(color = factor(num_ancestries)), size = 0.5, alpha = 0.5) +
  theme_classic() + custom_theme() +
  scale_x_continuous(limits = c(1.3,4.2)) +
  scale_y_continuous(limits = c(-3,4.7)) +
  labs(x = TeX("$log_{10}$ (CPU time/s)"), y = TeX("$log_{10}$(Num variants)"),
       color = "Num ancestries", title = "MACHINE")

p2 = ggplot(filter(data, method == "MESuSiE"), 
            aes(x = log10(num_variants), y = log10(time))) +
  geom_point(aes(color = factor(num_ancestries)), size = 0.5, alpha = 0.5) +
  theme_classic() + custom_theme() +
  scale_x_continuous(limits = c(1.3,4.2)) +
  scale_y_continuous(limits = c(-3,4.7)) +
  labs(x = TeX("$log_{10}$ (CPU time/s)"), y = TeX("$log_{10}$(Num variants)"),
       color = "Num ancestries", title = "MESuSiE")

pdf("computing_time.pdf", width = 9, height = 5.1, useDingbats = F, bg = "white")
ggarrange(p1, p2, nrow = 1, ncol = 2, common.legend = T, legend = "bottom",
          align = "hv")
dev.off()


model_info = foreach(method = c("MACHINE", "MESuSiE"), .combine = "rbind") %do%
{
  as.data.frame(summary(lm(log10(time) ~ log10(num_variants) + factor(num_ancestries), 
                           data = data[data$method == method,]))$coefficients) %>%
    cbind(data.frame("method" = method, "variable" = rownames(.)))
}[,c(5,6,1:4)]

saveRDS(model_info, "computing_time_model_info.RData")
