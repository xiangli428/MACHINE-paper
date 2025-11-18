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
library(igraph)

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

data = foreach(db = dbs, .combine = "rbind") %do%
{
  foreach(pheno = phenos, .combine = "rbind") %do%
  {
    foreach(method = names(methods), .combine = "rbind") %do%
    {
      output = data.frame("db" = db, "pheno" = pheno, "method" = method,
                          "EUR_AFR_EAS" = 0, "EUR_AFR" = 0, "EUR_EAS" = 0,
                          "AFR_EAS" = 0, "EUR" = 0, "AFR" = 0, "EAS" = 0)
      
      CS_overlap = read.delim(sprintf("results/%s/%s/%s/CS_overlap.txt", 
                                      methods[method], pheno, db))
      if(nrow(CS_overlap) > 0)
      {
        chr_block = unique(CS_overlap[,3:4])
        CS_overlap[,5:7][is.na(CS_overlap[,5:7])] = ""

        CS_overlap_union = foreach(l = 1:nrow(chr_block), .combine = "rbind") %dopar%
        {
          i = chr_block[l,1]
          block = chr_block[l,2]

          if(method %in% c("MACHINE + g-LDSC", "MACHINE", "MESuSiE"))
          {
            CS_info = read.delim(sprintf("%s/%s/%s/chr%d/%d/CS_info.txt",
                                         methods[method], pheno, db, i, block))
          } else {
            CS_info = foreach(pid = pids, .combine = "rbind") %do%
            {
              if(file.exists(sprintf("%s/%s/%s/%s/chr%d/%d/CS_info.txt",
                                     methods[method], pheno, db, pid, i, block)))
              read.delim(sprintf("%s/%s/%s/%s/chr%d/%d/CS_info.txt",
                                 methods[method], pheno, db, pid, i, block))
            }
          }
          CS_info$CS.tags[is.na(CS_info$CS.tags)] = ""

          sub = filter(CS_overlap, chromosome == chr_block[l,1] &
                         block == chr_block[l,2])[,5:7]
          vs_names = unique(c(sub[,1],sub[,2],sub[,3])) %>% setdiff("")
          vs = 1:length(vs_names)
          names(vs) = vs_names
          g = make_undirected_graph(n = length(vs), edges = NULL)
          for(j in 1:nrow(sub))
          {
            if(all(sub[j,] != ""))
            {
              g %<>% add_edges(edges = c(vs[sub[j,1]],vs[sub[j,2]],vs[sub[j,1]],
                                         vs[sub[j,3]],vs[sub[j,2]],vs[sub[j,3]]))
            } else {
              Ks = which(sub[j,] != "")
              g %<>% add_edges(edges = c(vs[sub[j,Ks[1]]],vs[sub[j,Ks[2]]]))
            }
          }

          cc = components(g)
          foreach(n = 1:cc$no, .combine = "rbind") %do%
          {
            ids = vs_names[which(cc$membership == n)]
            tmp = data.frame("pheno" = pheno,
                             "db" = db,
                             "chromosome" = i,
                             "block" = block)
            for(pid in pids)
            {
              tmp[[sprintf("CS.%s",pid)]] = paste(grep(pid, ids, value = T),
                                                  collapse = ',')
            }

            set = foreach(s = match(ids, CS_info$CS), .combine = "union") %do%
            {
              c(strsplit(CS_info$CS.SNP[s], ',')[[1]],
                strsplit(CS_info$CS.tags[s], ',')[[1]])
            }

            tmp$union.size = length(set)
            tmp$union.SNP = paste(set, collapse = ',')
            tmp
          }
        }
        
        write_delim(CS_overlap_union, sprintf(
          "results/%s/%s/%s/CS_overlap_union.txt", methods[method], pheno, db),
          delim = '\t')
        
        output$EUR_AFR_EAS = sum((CS_overlap_union$CS.EUR != '') & 
                                   (CS_overlap_union$CS.AFR != '') &
                                   (CS_overlap_union$CS.EAS != ''))
        output$EUR_AFR = sum((CS_overlap_union$CS.EUR != '') & 
                               (CS_overlap_union$CS.AFR != '') &
                               (CS_overlap_union$CS.EAS == ''))
        output$EUR_EAS = sum((CS_overlap_union$CS.EUR != '') & 
                               (CS_overlap_union$CS.AFR == '') &
                               (CS_overlap_union$CS.EAS != ''))
        output$AFR_EAS = sum((CS_overlap_union$CS.EUR == '') & 
                               (CS_overlap_union$CS.AFR != '') &
                               (CS_overlap_union$CS.EAS != ''))
      }
      
      if(method %in% c("MACHINE + g-LDSC", "MACHINE", "MESuSiE"))
      {
        CS_info = read.delim(sprintf("results/%s/%s/%s/CS_info.txt", 
                                     methods[method], pheno, db))
        for(pid in pids)
        {
          output[,pid] = length(setdiff(grep(pid, CS_info$CS, value = T), 
                                        CS_overlap[,sprintf("CS.%s",pid)]))
        }
      } else {
        for(pid in pids)
        {
          CS_info = read.delim(sprintf("results/%s/%s/%s/%s/CS_info.txt", 
                                       methods[method], pheno, db, pid))
          output[,pid] = length(setdiff(CS_info$CS, 
                                        CS_overlap[sprintf("CS.%s",pid)]))
        }
      }
      
      output
    }
  }
}

data$db %<>% factor(levels = dbs)
data$pheno %<>% factor(levels = phenos)
data$method %<>% factor(levels = names(methods))

saveRDS(data, "CS_overlap_num.RData")

data = readRDS("CS_overlap_num.RData")

data %<>% gather(pops, nCSs, -db, -pheno, -method)
data$pops %<>% gsub("_", " & ", .)
data$pops[grep("&", data$pops, invert = T)] %<>% sprintf("%s-specific", .)
data$pops %<>% factor(levels = c("EUR & AFR & EAS", "EUR & AFR", "EUR & EAS",
                                 "AFR & EAS", "EUR-specific", "AFR-specific", 
                                 "EAS-specific"))

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
    panel.grid = element_blank(),
    panel.spacing = unit(0.05, "in"),
    panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA))
}

pdf("CS_overlap_num.pdf", width = 9, height = 9.93, onefile = T, bg = "white")
for(db in dbs)
{
  print(ggplot(data[data$db == db,], aes(x = method, y = nCSs)) +
          geom_bar(aes(fill = method), stat = "identity", width = 0.8,
                   position = position_dodge2()) +
          geom_text(aes(label = nCSs), size = 2, vjust = -0.3) +
          facet_grid(pops ~ pheno, scales = "free_y") +
          theme_classic() + custom_theme() +
          scale_y_continuous(expand = expansion(mult = c(0,0.098))) +
          scale_fill_manual(values = hue_pal()(13)[c(1,3,4,8,10:13)]) +
          guides(fill = guide_legend(nrow = 1)) +
          labs(x = NULL, y = "Number of CSs", fill = "", title = db))
}
dev.off()
