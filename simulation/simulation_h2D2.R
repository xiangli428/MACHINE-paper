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
library(scales)
library(h2D2)

setwd("simulation")

select_block_info = read.delim("data/select_block_info.txt")
select_block_info %<>% arrange(num_variants)
select_block_info[134:200,] %<>% arrange(-num_variants)
select_block = sort(select_block_info$block)

n_causal = data.frame("n" = c(5,3,1),
                      "n1" = c(0,1,2),
                      "n2" = c(0,1,2))

N1 = 2e5
N2_seq = c(2e4,2e5)
suffix = c("N1-200000", sprintf("N2-%d", N2_seq))
idx = c(1,2,2)

group_bound = data.frame("l" = c(0,0.1,0.5,0.9),
                         "u" = c(0.1,0.5,0.9,1.1))

registerDoParallel(67)

Dir = "h2D2"

for(s in 1:3)
{
  foreach(block = select_block_info$block, .combine = "c", .inorder = F) %dopar%
  {
    dir.create(sprintf("results/setting_%s/%s/%s",
                       s, Dir, block), recursive = T)

    data_full = readRDS(sprintf("data/setting_%s/%d.RData", s, block))
    R = readRDS(sprintf("data/LD/%d.RData", block))

    for(k in 1:3)
    {
      data = data_full[,c(1:12,12+k)]

      h2D2 = Createh2D2Object(data[,13],
                              R[[idx[k]]],
                              c(N1, N2_seq)[k],
                              data$rsid,
                              in_sample_LD = T,
                              a = 0.005,
                              b = NULL)

      tic = Sys.time()
      h2D2 = h2D2_MCMC(h2D2, mcmc_n = 5500, burn_in = 500, thin = 2, get_CS = F)
      while(max(h2D2@mcmc_samples$PSRF_beta) > 1.2)
      {
        h2D2 = h2D2_MCMC(h2D2, mcmc_n = 1000, burn_in = 1000, thin = 2, get_CS = F)
      }
      h2D2@CS = h2D2_CS(h2D2)
      toc = Sys.time()

      saveRDS(h2D2, file = sprintf(
        "results/setting_%s/%s/%s/h2D2_%s.RData",
        s, Dir, block, suffix[k]))

      data$CL = h2D2@CL

      data$CS = NA
      if(length(h2D2@CS$sets) > 0)
      {
        for(l in 1:length(h2D2@CS$sets))
        {
          data$CS[h2D2@CS$sets[[l]]] = l
        }
      }

      saveRDS(data, file = sprintf(
        "results/setting_%s/%s/%s/data_%s.RData",
        s, Dir, block, suffix[k]))
      
      result_CS95 = data.frame("block" = block,
                               "coverage" = 0,
                               "size" = 0,
                               "min.abs.corr" = 0)[-1,]
      if(length(h2D2@CS$sets) > 0)
      {
        for(l in 1:length(h2D2@CS$sets))
        {
          result_CS95 %<>% rbind(
            data.frame("block" = block,
                       "coverage" = sum(
                         data[h2D2@CS$sets[[l]], 8+idx[k]]),
                       "size" = length(h2D2@CS$sets[[l]]),
                       "min.abs.corr" = h2D2@CS$purity$min.abs.corr[l]))
        }
      }
      
      write_delim(result_CS95, sprintf(
        "results/setting_%s/%s/%s/result_CS95_%s.txt",
        s, Dir, block, suffix[k]), delim = '\t')
      
      result = data.frame(
        "block" = block,
        "max_PSRF" = max(h2D2@mcmc_samples$PSRF_beta),
        "mcmc_n" = h2D2@mcmc_samples$n_samples + h2D2@mcmc_samples$n_burnin,
        "time" = difftime(toc, tic, units = "secs"))
      
      if(k > 1)
      {
        m = k - 1
        
        data = list()
        data[[1]] = readRDS(sprintf(
          "results/setting_%s/%s/%s/data_%s.RData",
          s, Dir, block, suffix[1]))
        data[[2]] = readRDS(sprintf(
          "results/setting_%s/%s/%s/data_%s.RData",
          s, Dir, block, suffix[m+1]))
        
        result_CS95_cross = data.frame("block" = block,
                                       "coverage" = 0,
                                       "size" = 0,
                                       "min.abs.corr" = 0)[-1,]
        result_CS95_cross %<>% rbind(
          foreach(k = 1:2, .combine = "rbind") %do%
            {
              if(length(na.omit(data[[k]]$CS)) > 0)
              {
                foreach(l = unique(na.omit(data[[k]]$CS)),
                        .combine = "rbind") %do%
                  {
                    set = which(data[[k]]$CS == l)
                    if(length(set) == 1)
                    {
                      purity = 1
                    } else {
                      R_sub = R[[k]][set, set]
                      purity = min(abs(R_sub[upper.tri(R_sub)]))
                    }
                    data.frame("block" = block,
                               "coverage" = sum(data[[1]][set,9] | data[[1]][set,10]),
                               "size" = length(set),
                               "min.abs.corr" = purity)
                  }
              }
            }
        )
        
        write_delim(result_CS95_cross, sprintf(
          "results/setting_%s/%s/%s/result_CS95_N2-%d_cross.txt",
          s, Dir, block, N2_seq[m]), delim = '\t')
        
        result_CS95_shared = data.frame("block" = block,
                                        "coverage" = 0,
                                        "size" = 0,
                                        "min.abs.corr" = 0)[-1,]
        if(length(h2D2[[1]]@CS$sets) > 0 & length(h2D2[[2]]@CS$sets) > 0)
        {
          for(l1 in 1:length(h2D2[[1]]@CS$sets))
          {
            set1 = h2D2[[1]]@CS$sets[[l1]]
            for(l2 in 1:length(h2D2[[2]]@CS$sets))
            {
              set2 = h2D2[[2]]@CS$sets[[l2]]
              if(length(intersect(set1,set2)) > 0)
              {
                result_CS95_shared %<>% rbind(data.frame(
                  "block" = block,
                  "coverage" = sum(data[[1]][intersect(set1,set2),9] &
                                     data[[1]][intersect(set1,set2),10]),
                  "size" = length(union(set1,set2)),
                  "min.abs.corr" = min(h2D2[[1]]@CS$purity$min.abs.corr[l1],
                                       h2D2[[2]]@CS$purity$min.abs.corr[l2]))
                )
              }
            }
          }
        }
        
        write_delim(result_CS95_shared, sprintf(
          "results/setting_%s/%s/%s/result_CS95_N2-%d_shared.txt",
          s, Dir, block, N2_seq[m]), delim = '\t')
        
        result$CS95_causal = sum((!is.na(data[[1]]$CS) | !is.na(data[[2]]$CS)) &
                                   (data[[1]]$causal_1 | data[[1]]$causal_2))
        result$CS95_causal_1 = sum(!is.na(data[[1]]$CS) & data[[1]]$causal_1)
        result$CS95_causal_2 = sum(!is.na(data[[2]]$CS) & data[[1]]$causal_2)
        result$CS95_causal_0 = sum((!is.na(data[[1]]$CS) & !is.na(data[[2]]$CS)) &
                                     (data[[1]]$causal_1 & data[[1]]$causal_2))
      }

      write_delim(result, sprintf(
        "results/setting_%s/%s/%s/result_%s.txt",
        s, Dir, block, suffix[k]), delim = '\t')
    }
    
    NULL
  }
  
  for(k in 1:3)
  {
    results = foreach(block = select_block, .combine = "rbind") %dopar%
    {
      read.delim(sprintf(
        "results/setting_%s/%s/%s/result_%s.txt",
        s, Dir, block, suffix[k]))
    }
    write_delim(results, sprintf(
      "results/setting_%s/%s/results_%s.txt",
      s, Dir, suffix[k]), delim = '\t')

    results_CS95 = foreach(block = select_block, .combine = "rbind") %dopar%
    {
      read.delim(sprintf(
        "results/setting_%s/%s/%s/result_CS95_%s.txt",
        s, Dir, block, suffix[k]))
    }
    write_delim(results_CS95, sprintf(
      "results/setting_%s/%s/results_CS95_%s.txt",
      s, Dir, suffix[k]), delim = '\t')
    
    results_CS95_cross = foreach(block = select_block, .combine = "rbind") %dopar%
    {
      read.delim(sprintf(
        "results/setting_%s/%s/%s/result_CS95_N2-%d_cross.txt",
        s, Dir, block, N2_seq[m]))
    }
    write_delim(results_CS95_cross, sprintf(
      "results/setting_%s/%s/results_CS95_N2-%d_cross.txt",
      s, Dir, N2_seq[m]), delim = '\t')
    
    results_CS95_shared = foreach(block = select_block, .combine = "rbind") %dopar%
    {
      read.delim(sprintf(
        "results/setting_%s/%s/%s/result_CS95_N2-%d_shared.txt",
        s, Dir, block, N2_seq[m]))
    }
    write_delim(results_CS95_shared, sprintf(
      "results/setting_%s/%s/results_CS95_N2-%d_shared.txt",
      s, Dir, N2_seq[m]), delim = '\t')
  }
  
  # Calibration
  
  data_1 = foreach(block = select_block, .combine = "rbind") %dopar%
  {
    readRDS(sprintf(
      "results/setting_%s/%s/%s/data_%s.RData",
      s, Dir, block, suffix[1]))
  }
  
  for(m in 1:4)
  {
    data_2 = foreach(block = select_block, .combine = "rbind") %dopar%
    {
      readRDS(sprintf(
        "results/setting_%s/%s/%s/data_%s.RData",
        s, Dir, block, suffix[m+1]))
    }
    
    calibration = list()
    
    CL = pmax(data_1$CL, data_2$CL)
    calibration[["cross"]] = foreach(g = 1:4, .combine = "rbind") %dopar%
    {
      idx = which(CL >= group_bound$l[g] & CL < group_bound$u[g])
      n = length(idx)
      
      data.frame("group" = g,
                 "n" = n,
                 "Expected" = sum(CL[idx]) / n,
                 "Prop" = sum((data_1$causal_1 | data_2$causal_2)[idx]) / n)
    }
    
    CL = pmin(data_1$CL, data_2$CL)
    calibration[["shared"]] = foreach(g = 1:4, .combine = "rbind") %dopar%
    {
      idx = which(CL >= group_bound$l[g] & CL < group_bound$u[g])
      n = length(idx)
      
      data.frame("group" = g,
                 "n" = n,
                 "Expected" = sum(CL[idx]) / n,
                 "Prop" = sum((data_1$causal_1 & data_2$causal_2)[idx]) / n)
    }
    
    CL = data_1$CL
    calibration[["pop1"]] = foreach(g = 1:4, .combine = "rbind") %dopar%
    {
      idx = which(CL >= group_bound$l[g] & CL < group_bound$u[g])
      n = length(idx)
      
      data.frame("group" = g,
                 "n" = n,
                 "Expected" = sum(CL[idx]) / n,
                 "Prop" = sum(data_1$causal_1[idx]) / n)
    }
    
    CL = data_2$CL
    calibration[["pop2"]] = foreach(g = 1:4, .combine = "rbind") %dopar%
    {
      idx = which(CL >= group_bound$l[g] & CL < group_bound$u[g])
      n = length(idx)
      
      data.frame("group" = g,
                 "n" = n,
                 "Expected" = sum(CL[idx]) / n,
                 "Prop" = sum(data_2$causal_2[idx]) / n)
    }
    
    saveRDS(calibration, sprintf(
      "results/setting_%s/%s/calibration_N2-%d.RData",
      s, Dir, N2_seq[m]))
  }
}
