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

suffixes = c("N1-200000" = "N^{(1)}==200*k",
             "N2-20000" = "N^{(2)}==20*k",
             "N2-200000" = "N^{(2)}==200*k")
idx = c(1,2,2)

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

registerDoParallel(10)

#####
# UKBB

var_all = foreach(s = 1:3, .combine = "rbind") %dopar%
{
  data = foreach(block = select_block, .combine = "rbind", .inorder = F) %dopar%
  {
    readRDS(sprintf("../data/setting_%s/%d.RData", s, block))
  }
  
  gLDSC_var = readRDS(sprintf(
    "../gLDSC/gLDSC_results/UKBB/setting_%d/gLDSC_var.RData", s))[,c(1,2,5)]
  for(k in 1:3)
  {
    gLDSC_var[,k] %<>% pmax(max(.) / 20)
  }
  
  polyfun_var = readRDS(sprintf(
    "../polyfun_output/UKBB/setting_%d/polyfun_var.RData", s))[,c(1,2,5)]
  
  foreach(k = 1:3, .combine = "rbind") %do%
  {
    data.frame("scenario" = s,
               "suffix" = suffixes[k],
               "causal" = data[,8+idx[k]],
               "beta" = data[,10+idx[k]],
               "g-LDSC" = gLDSC_var[,k],
               "PolyFun" = polyfun_var[,k],
               check.names = F)
  }
}

p1 = ggplot(var_all, aes(x = causal, y = `g-LDSC`, fill = causal)) +
  geom_violin(scale = "width") +
  facet_grid(suffix ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios, 
                                 suffix = label_parsed)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(labels = NULL) +
  scale_fill_discrete(labels = c("FALSE" = "Non-causal", "TRUE" = "Causal")) +
  labs(x = NULL, y = TeX("$\\hat{\\sigma}$"), title = "g-LDSC")

p2 = ggplot(var_all, aes(x = causal, y = PolyFun, fill = causal)) +
  geom_violin(scale = "width") +
  facet_grid(suffix ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios, 
                                 suffix = label_parsed)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(labels = NULL) +
  scale_fill_discrete(labels = c("FALSE" = "Non-causal", "TRUE" = "Causal")) +
  labs(x = NULL, y = TeX("$\\hat{\\sigma}$"), title = "PolyFun")

ggsave("var_violin_UKBB.pdf",
       ggarrange(p1,p2, 
                 nrow = 2, ncol = 1, common.legend = T, legend = "bottom",
                 heights = c(1,1), labels = c("a","b"),
                 font.label = list(size = 12, color = "black", face = "bold", 
                                   family = NULL), align = "v"),
       device = "pdf", width = 6.3, height = 12.2, units = "in", bg = "white")


fun = function(x)
{
  c(wilcox.test(x[x[,3],5], x[!x[,3],5], alternative = "greater")$p.value,
    wilcox.test(x[x[,3],6], x[!x[,3],6], alternative = "greater")$p.value)
}


wilcox_p = foreach(s = 1:3, .combine = "rbind") %do%
{
  foreach(k = 1:3, .combine = "rbind") %do%
  {
    c = filter(var_all, scenario == s & suffix == suffixes[k]) %>% fun
    
    data.frame("scenario" = s,
               "setting" = names(suffixes)[k],
               "g-LDSC" = c[1],
               "PolyFun" = c[2],
               check.names = F)
  }
}

wilcox_p %<>% gather(method, p, -scenario, -suffix)
wilcox_p$method %<>% factor(levels = c("g-LDSC","PolyFun"))
wilcox_p$setting %<>% factor(levels = names(suffixes))

write_delim(wilcox_p, "var_wilcox_p_UKBB.txt", delim = '\t')
saveRDS(wilcox_p, "var_wilcox_p_UKBB.RData")

wilcox_p = readRDS("var_wilcox_p_UKBB.RData")
wilcox_p$setting = suffixes[wilcox_p$setting]
wilcox_p$setting %<>% factor(levels = suffixes)

pdf("var_wilcox_p_UKBB.pdf", width = 6.3, height = 6.42, bg = "white")
ggplot(wilcox_p, aes(x = method, y = -log10(p))) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.6,  
           position = position_dodge2()) +
  facet_grid(setting ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios, 
                                 setting = label_parsed)) +
  theme_classic() + custom_theme() + 
  scale_x_discrete(labels = NULL) +
  labs(x = NULL, y = TeX("$-log_{10}(\\italic(P))$"))
dev.off()



#####
# 1kG

var_all = foreach(s = 1:3, .combine = "rbind") %dopar%
{
  data = foreach(block = select_block, .combine = "rbind", .inorder = F) %dopar%
  {
    readRDS(sprintf("../data/setting_%s/%d.RData", s, block))
  }
  
  gLDSC_var = readRDS(sprintf(
    "../gLDSC/gLDSC_results/1kg/setting_%d/gLDSC_var.RData", s))[,c(1,2,5)]
  for(k in 1:3)
  {
    gLDSC_var[,k] %<>% pmax(max(.) / 20)
  }
  
  polyfun_var = readRDS(sprintf(
    "../polyfun_output/1kg/setting_%d/polyfun_var.RData", s))[,c(1,2,5)]
  
  foreach(k = 1:3, .combine = "rbind") %do%
  {
    data.frame("scenario" = s,
               "suffix" = suffixes[k],
               "causal" = data[,8+idx[k]],
               "beta" = data[,10+idx[k]],
               "g-LDSC" = gLDSC_var[,k],
               "PolyFun" = polyfun_var[,k],
               check.names = F)
  }
}

p1 = ggplot(var_all, aes(x = causal, y = `g-LDSC`, fill = causal)) +
  geom_violin(scale = "width") +
  facet_grid(suffix ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios, 
                                 suffix = label_parsed)) +
  theme_classic() + custom_theme() +
  scale_x_discrete(labels = NULL) +
  scale_fill_discrete(labels = c("FALSE" = "Non-causal", "TRUE" = "Causal")) +
  labs(x = NULL, y = TeX("$\\hat{\\sigma}$"), title = "g-LDSC")

p2 = ggplot(var_all, aes(x = causal, y = PolyFun, fill = causal)) +
  geom_violin(scale = "width") +
  facet_grid(suffix ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios, 
                                 suffix = label_parsed)) +
  theme_classic() + custom_theme() +
  scale_x_discrete(labels = NULL) +
  scale_fill_discrete(labels = c("FALSE" = "Non-causal", "TRUE" = "Causal")) +
  labs(x = NULL, y = TeX("$\\hat{\\sigma}$"), title = "PolyFun")

ggsave("var_violin_1kG.pdf",
       ggarrange(p1,p2, 
                 nrow = 2, ncol = 1, common.legend = T, legend = "bottom",
                 heights = c(1,1), labels = c("a","b"),
                 font.label = list(size = 12, color = "black", face = "bold", 
                                   family = NULL), align = "v"),
       device = "pdf", width = 6.3, height = 12.2, units = "in", bg = "white")


fun = function(x)
{
  c(wilcox.test(x[x[,3],5], x[!x[,3],5], alternative = "greater")$p.value,
    wilcox.test(x[x[,3],6], x[!x[,3],6], alternative = "greater")$p.value)
}


wilcox_p = foreach(s = 1:3, .combine = "rbind") %do%
{
  foreach(k = 1:3, .combine = "rbind") %do%
  {
    c = filter(var_all, scenario == s & suffix == suffixes[k]) %>% fun
    
    data.frame("scenario" = s,
               "setting" = names(suffixes)[k],
               "g-LDSC" = c[1],
               "PolyFun" = c[2],
               check.names = F)
  }
}

wilcox_p %<>% gather(method, p, -scenario, -suffix)
wilcox_p$method %<>% factor(levels = c("g-LDSC","PolyFun"))
wilcox_p$setting %<>% factor(levels = names(suffixes))

write_delim(wilcox_p, "var_wilcox_p_1kG.txt", delim = '\t')
saveRDS(wilcox_p, "var_wilcox_p_1kG.RData")

wilcox_p = readRDS("var_wilcox_p_1kG.RData")
wilcox_p$setting = suffixes[wilcox_p$setting]
wilcox_p$setting %<>% factor(levels = suffixes)

pdf("var_wilcox_p_1kG.pdf", width = 6.3, height = 6.42, bg = "white")
ggplot(wilcox_p, aes(x = method, y = -log10(p))) +
  geom_bar(aes(fill = method), stat = "identity", width = 0.6,  
           position = position_dodge2()) +
  facet_grid(setting ~ scenario, scales = "free",
             labeller = labeller(scenario = scenarios, 
                                 setting = label_parsed)) +
  theme_classic() + custom_theme() +
  scale_x_discrete(labels = NULL) +
  labs(x = NULL, y = TeX("$-log_{10}(\\italic(P))$"))
dev.off()
