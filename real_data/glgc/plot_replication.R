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
library(RColorBrewer)
library(latex2exp)
library(scales)
library(Hmisc)

setwd("real_data/glgc")

dbs = c("UKBB", "GLGC")
phenos = c("HDL", "LDL", "TG", "TC")
pids = c("EUR", "AFR", "EAS")

methods = c("MACHINE + g-LDSC" = "MACHINE-gLDSC",
            "MACHINE" = "MACHINE",
            "MESuSiE" = "MESuSiE",
            "h2-D2 + g-LDSC" = "h2D2-gLDSC",
            "h2-D2" = "h2D2",
            "SuSiE" = "SuSiE",
            "RSparsePro" = "RSparsePro",
            "CARMA" = "CARMA")

range_UKBB = list("[0.5,0.6)" = c(0.5,0.6),
                  "[0.6,0.7)" = c(0.6,0.7),
                  "[0.7,0.8)" = c(0.7,0.8),
                  "[0.8,0.9)" = c(0.8,0.9),
                  "[0.9,1]" = c(0.9,1.1))
range_GLGC = list("[0,0.1)" = c(0,0.1),
                  "[0.1,0.5)" = c(0.1,0.5),
                  "[0.5,0.9)" = c(0.5,0.9),
                  "[0.9,1]" = c(0.9,1.1))

registerDoParallel(30)

data = foreach(pheno = phenos, .combine = "rbind") %do%
{
  foreach(method = names(methods), .combine = "rbind") %do%
  {
    if(method %in% c("MACHINE + g-LDSC", "MACHINE", "MESuSiE"))
    {
      results = read.delim(sprintf("results/%s/%s/%s/results.txt", 
                                   methods[method], pheno, "UKBB"))
      
      data = foreach(db = dbs) %do%
      {
        foreach(l = which(
          results$num_variants_0.5.EUR > 0 | results$num_variants_0.5.AFR > 0 |
            results$num_variants_0.5.EAS > 0), .combine = "rbind") %dopar%
        {
          i = results$chromosome[l]
          block = results$block[l]
          
          read.delim(sprintf("%s/%s/%s/chr%d/%d/variant_info.txt.gz", 
                             methods[method], pheno, db, i, block))
        }
      }
      
      foreach(pid = pids, .combine = "rbind") %do%
      {
        if(method == "MESuSiE")
        {
          CL_1 = na.omit(data[[1]][,sprintf("PIP.%s", pid)])
          CL_2 = na.omit(data[[2]][,sprintf("PIP.%s", pid)])
        } else {
          CL_1 = na.omit(data[[1]][,sprintf("CL.%s", pid)])
          CL_2 = na.omit(data[[2]][,sprintf("CL.%s", pid)])
        }
        
        foreach(g1 = names(range_UKBB), .combine = "rbind") %do%
        {
          foreach(g2 = names(range_GLGC), .combine = "rbind") %do%
          {
            data.frame(
              "pheno" = pheno, "pid" = pid, "method" = method,
              "group_UKBB" = g1, "group_GLGC" = g2,
              "nVar" = sum(CL_1 >= range_UKBB[[g1]][1] & 
                              CL_1 <= range_UKBB[[g1]][2] &
                              CL_2 >= range_GLGC[[g2]][1] & 
                              CL_2 <= range_GLGC[[g2]][2]))
          }
        }
      }
    } else {
      foreach(pid = pids, .combine = "rbind") %do%
      {
        results = read.delim(sprintf("results/%s/%s/%s/%s/results.txt", 
                                     methods[method], pheno, "UKBB", pid))
        
        data = foreach(db = dbs) %do%
        {
          foreach(l = which(results$num_variants_0.5 > 0), 
                  .combine = "rbind") %dopar%
          {
            i = results$chromosome[l]
            block = results$block[l]
            
            read.delim(sprintf("%s/%s/%s/%s/chr%d/%d/variant_info.txt.gz", 
                               methods[method], pheno, db, pid, i, block))
          }
        }
        
        if(method %in% c("h2-D2 + g-LDSC", "h2-D2"))
        {
          CL_1 = data[[1]]$CL
          CL_2 = data[[2]]$CL
        } else {
          CL_1 = data[[1]]$PIP
          CL_2 = data[[2]]$PIP
        }
        
        foreach(g1 = names(range_UKBB), .combine = "rbind") %do%
        {
          foreach(g2 = names(range_GLGC), .combine = "rbind") %do%
          {
            data.frame(
              "pheno" = pheno, "pid" = pid, "method" = method,
              "group_UKBB" = g1, "group_GLGC" = g2,
              "nVar" = sum(CL_1 >= range_UKBB[[g1]][1] & 
                              CL_1 <= range_UKBB[[g1]][2] &
                              CL_2 >= range_GLGC[[g2]][1] & 
                              CL_2 <= range_GLGC[[g2]][2]))
          }
        }
      }
    }
  }
}

data$pheno %<>% factor(levels = phenos)
data$pid %<>% factor(levels = pids)
data$method %<>% factor(levels = names(methods))
data$group_UKBB %<>% factor(levels = names(range_UKBB))
data$group_GLGC %<>% factor(levels = names(range_GLGC))

saveRDS(data, "replication.RData")

data = readRDS("replication.RData")

custom_theme = function()
{
  theme(
    aspect.ratio = 3/5,
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
    panel.spacing = unit(0.05, "in"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA))
}

bk = c("EUR" = 50, "AFR" = 5, "EAS" = 2)

pdf("replication.pdf", width = 9, height = 11.27, onefile = T, bg = "white")
for(pid in pids)
{
  print(ggplot(
    data[data$pid == pid,], aes(x = nVar, y = group_UKBB, fill = group_GLGC)) +
      geom_bar(stat = "identity") +
      facet_grid(method ~ pheno, scales = "free_x") +
      theme_classic() + custom_theme() +
      scale_x_continuous(breaks = seq(from = 0, to = 400, by = bk[pid])) +
      scale_fill_manual(values = brewer.pal(5, "Blues")[2:5]) +
      labs(x = "Number of variants", y = "UKBB CL or PIP bins", title = pid,
           fill = "GLGC CL or PIP bins"))
}
dev.off()
