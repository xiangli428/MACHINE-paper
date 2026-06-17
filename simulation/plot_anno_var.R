options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(tidyr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(scales)

setwd("simulation/results")

select_block = readRDS("../data/select_block.RData")

N1 = 2e5
N2_seq = c(2e4,2e5)

scenarios = c("1" = "5 shared causal variants", 
              "2" = "3 shared causal variants",
              "3" = "1 shared causal variant")

suffixes = c("N1-200000" = "N^{(EUR)}==200*k",
             "N2-20000" = "N^{(AFR)}==20*k",
             "N2-200000" = "N^{(AFR)}==200*k",
             "N3-20000" = "N^{(EAS)}==20*k",
             "N3-200000" = "N^{(EAS)}==200*k")
idx = c(1,2,2,2,2)

custom_theme = function()
{
  theme(
    aspect.ratio = 1,
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

registerDoParallel(20)


#####
var_all = foreach(s = 1:3, .combine = "rbind") %dopar%
{
  data = foreach(block = select_block, .combine = "rbind", .inorder = F) %dopar%
  {
    readRDS(sprintf("../data/setting_%s/%d.RData", s, block))
  }
  
  foreach(ld = c("UKBB","1kg"), .combine = "rbind") %do%
  {
    gLDSC_var = readRDS(sprintf(
      "../gLDSC/gLDSC_results/%s/setting_%d/gLDSC_var.RData", ld, s))
    for(k in 1:5)
    {
      gLDSC_var[,k] %<>% pmax(max(.) / 20)
    }
    
    polyfun_var = readRDS(sprintf(
      "../polyfun/output/%s/setting_%d/polyfun_var.RData", ld, s))
    
    foreach(k = 1:5, .combine = "rbind") %do%
    {
      data.frame("scenario" = s,
                 "suffix" = suffixes[k],
                 "causal" = data[,16+idx[k]],
                 "LD" = ld,
                 "g-LDSC" = gLDSC_var[,k],
                 "PolyFun" = polyfun_var[,k],
                 check.names = F)
    }
  }
}
var_all$suffix %<>% factor(levels = suffixes)
var_all$scenario = factor(scenarios[var_all$scenario], levels = scenarios)


pdf("var_violin.pdf", width = 9, height = 7, useDingbats = F, bg = "white")
for(ld in c("UKBB","1kg"))
{
  p1 = ggplot(filter(var_all, LD == ld), aes(x = causal, y = `g-LDSC`, fill = causal)) +
    geom_violin(scale = "width") +
    facet_grid(suffix ~ scenario, scales = "free", labeller = labeller(
      scenario = label_wrap_gen(width = 16), suffix = label_parsed)) +
    theme_classic() + custom_theme() + 
    scale_x_discrete(labels = NULL) +
    scale_fill_discrete(labels = c("FALSE" = "Non-causal", "TRUE" = "Causal")) +
    labs(x = NULL, y = TeX("$\\hat{\\sigma}$"), title = "g-LDSC")
  
  p2 = ggplot(filter(var_all, LD == ld), aes(x = causal, y = PolyFun, fill = causal)) +
    geom_violin(scale = "width") +
    facet_grid(suffix ~ scenario, scales = "free", labeller = labeller(
      scenario = label_wrap_gen(width = 16), suffix = label_parsed)) +
    theme_classic() + custom_theme() + 
    scale_x_discrete(labels = NULL) +
    scale_fill_discrete(labels = c("FALSE" = "Non-causal", "TRUE" = "Causal")) +
    labs(x = NULL, y = TeX("$\\hat{\\sigma}$"), title = "PolyFun")
  
  print(ggarrange(p1, p2, 
                  nrow = 1, ncol = 2, common.legend = T, legend = "bottom",
                  labels = c("a","b"),
                  font.label = list(size = 12, color = "black", face = "bold", 
                                    family = NULL), align = "v"))
}
dev.off()


fun = function(x)
{
  c(wilcox.test(x[x[,3],5], x[!x[,3],5], alternative = "greater", correct = F)$p.value,
    wilcox.test(x[x[,3],6], x[!x[,3],6], alternative = "greater", correct = F)$p.value)
}

wilcox_p = foreach(ld = c("UKBB","1kg"), .combine = "rbind") %do%
{
  foreach(s = 1:3, .combine = "rbind") %do%
  {
    foreach(k = 1:5, .combine = "rbind") %do%
    {
      pval = filter(var_all, scenario == scenarios[s] & suffix == suffixes[k] &
                   LD == ld) %>% fun
      
      data.frame("LD" = ld,
                 "scenario" = s,
                 "setting" = names(suffixes)[k],
                 "g-LDSC" = pval[1],
                 "PolyFun" = pval[2],
                 check.names = F)
    }
  }
}

wilcox_p %<>% gather(method, p, -LD, -scenario, -setting)
wilcox_p$LD = c("UKBB" = "In-sample LD", "1kg" = "1kG LD")[wilcox_p$LD]
wilcox_p$method %<>% factor(levels = c("g-LDSC","PolyFun"))
wilcox_p$setting %<>% factor(levels = names(suffixes))

saveRDS(wilcox_p, "var_wilcox_p.RData")

wilcox_p$setting = suffixes[wilcox_p$setting]
wilcox_p$setting %<>% factor(levels = suffixes)
wilcox_p$scenario = scenarios[wilcox_p$scenario]
wilcox_p$scenario %<>% factor(levels = scenarios)

pdf("var_wilcox_p.pdf", width = 9, height = 7.3, bg = "white")

p1 = ggplot(filter(wilcox_p, LD == "In-sample LD"), aes(x = method, y = -log10(p))) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.6,  
           position = position_dodge2()) +
  facet_grid(setting ~ scenario, scales = "free",
             labeller = labeller(scenario = label_wrap_gen(width = 16), 
                                 setting = label_parsed)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(labels = NULL) +
  labs(x = NULL, y = TeX("$-log_{10}(\\italic(P))$"), title = "In-sample LD")

p2 = ggplot(filter(wilcox_p, LD == "1kG LD"), aes(x = method, y = -log10(p))) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.6,  
           position = position_dodge2()) +
  facet_grid(setting ~ scenario, scales = "free",
             labeller = labeller(scenario = label_wrap_gen(width = 16), 
                                 setting = label_parsed)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(labels = NULL) +
  labs(x = NULL, y = TeX("$-log_{10}(\\italic(P))$"), title = "1kG LD")

print(ggarrange(p1, p2, 
                nrow = 1, ncol = 2, common.legend = T, legend = "bottom",
                labels = c("a","b"),
                font.label = list(size = 12, color = "black", face = "bold", 
                                  family = NULL), align = "v"))

dev.off()
