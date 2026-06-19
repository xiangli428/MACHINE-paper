options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(tidyr)
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

registerDoParallel(20)


output = foreach(block = select_block) %dopar%
{
  R = readRDS(sprintf("../data/LD/%d/LD_UKBB.RData", block))
  R_1kG = readRDS(sprintf("../data/LD/%d/LD_1kg.RData", block))
  
  r = foreach(k = 1:3, .combine = "cbind") %do%
  {
    R[[k]] %<>% as.matrix()
    R[[k]][upper.tri(R[[k]])]
  } %>% cbind(foreach(k = 1:3, .combine = "cbind") %do% {
    R_1kG[[k]] %<>% as.matrix()
    R_1kG[[k]][upper.tri(R_1kG[[k]])]})
  
  freq = foreach(k = 1:6, .combine = "cbind") %do%
  {
    table(factor(pmin(floor(abs(r[,k]) * 100) + 1, 100), levels = 1:100))
  }
  
  list(freq, t(r) %*% r)
}


freq = foreach(k = 1:200, .combine = "+") %dopar%
{
  output[[k]][[1]]
}

freq %<>% as.data.frame()
names(freq) = c("EUR_UKBB","AFR_UKBB","EAS_UKBB","EUR_1kG","AFR_1kG","EAS_1kG")
freq$r = seq(from = 0.005, to = 0.995, by = 0.01)
freq = freq[,c(7,1:6)]

freq %<>% gather(pid_LD, density, -r) %>% 
  separate(pid_LD, into = c("pid", "LD"), sep = '_') %>%
  mutate(density = density / 121698792 / 0.01)
freq$pid %<>% factor(levels = c("EUR","AFR","EAS"))
freq$LD %<>% factor(levels = c("UKBB","1kG"))

saveRDS(freq, "LD_distribution.RData")

rr = foreach(k = 1:200, .combine = "+") %dopar%
{
  output[[k]][[2]]
}
LD_similarity = data.frame(t(rr / sqrt(diag(rr))) / sqrt(diag(rr)))
LD_similarity$pid_LD1 = names(LD_similarity) = c(
  "EUR_UKBB","AFR_UKBB","EAS_UKBB","EUR_1kG","AFR_1kG","EAS_1kG")
LD_similarity %<>% gather(pid_LD2, similarity, -pid_LD1)
LD_similarity %<>% separate(pid_LD1, into = c("pid1", "LD1"), sep = '_') %>%
  separate(pid_LD2, into = c("pid2", "LD2"), sep = '_')
LD_similarity$pid1 %<>% factor(levels = c("EUR","AFR","EAS"))
LD_similarity$pid2 %<>% factor(levels = c("EAS","AFR","EUR"))
LD_similarity$LD1 %<>% factor(levels = c("UKBB","1kG"))
LD_similarity$LD2 %<>% factor(levels = c("UKBB","1kG"))

saveRDS(LD_similarity, "LD_similarity.RData")

LD_similarity$label = sprintf("%.2f", LD_similarity$similarity)


custom_theme = function()
{
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 7),  
    axis.title = element_text(size = 8),
    axis.ticks = element_line(linewidth = 0.25),
    axis.line = element_blank(),
    strip.text = element_text(size = 9),
    strip.background = element_rect(
      fill = "lightgray", color = "lightgray"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    plot.title = element_text(size = 10, hjust = 0.5),
    panel.spacing = unit(0.01, "in"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA))
}


p1 = ggplot(freq, aes(x = r, y = density, color = pid)) +
  geom_line(aes(linetype = LD), linewidth = 0.25) + theme_classic() + custom_theme() +
  scale_x_continuous(limits = c(0.15,1)) +
  scale_y_continuous(limits = c(0,0.8)) +
  labs(x = TeX("$|r|$"), y = "Density", title = NULL, color = NULL)

p2 = ggplot(LD_similarity, aes(x = pid1, y = pid2, fill = pmin(similarity,1))) +
  geom_tile() + theme_classic() + custom_theme() +
  geom_text(aes(label = label), color = "black", parse = T, size = 3) +
  facet_grid(LD2 ~ LD1, scales = "free") +
  scale_fill_gradient(limits = c(0.4,1), low = "white", high = "blue") +
  labs(x = NULL, y = NULL, title = NULL, fill = "Cosine similarity")


pdf("LD_compare.pdf", width = 9, height = 5.5, useDingbats = F, bg = "white")
p1 + p2 + plot_layout() &
  theme(legend.position = "bottom", 
        plot.tag = element_text(size = 12, face = "bold")) &
  plot_annotation(tag_level = "a")
dev.off()
