options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(tidyr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)
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
suffixes = c("N1-200000", sprintf("N2-%d", N2_seq), sprintf("N3-%d", N2_seq))
idx = c(1,2,2,3,3)

group_bound = data.frame("l" = c(0,0.1,0.5,0.9),
                         "u" = c(0.1,0.5,0.9,1.1))

registerDoParallel(67)

for(anno in c("gLDSC","polyfun"))
{
  Dir = sprintf("h2D2-%s", anno)
  
  for(s in 1:3)
  {
    if(anno == "gLDSC")
    {
      gLDSC_var = readRDS(sprintf(
        "gLDSC/gLDSC_results/UKBB/setting_%d/gLDSC_var.RData", s))
      for(k in 1:5)
      {
        gLDSC_var[,k] %<>% pmax(max(gLDSC_var[,k]) / 20)
      }
    } else {
      gLDSC_var = readRDS(sprintf(
        "polyfun/output/UKBB/setting_%s/polyfun_var.RData", s))
    }
    
    foreach(block = select_block_info$block, .combine = "c", .inorder = F) %dopar%
    {
      data_full = readRDS(sprintf("data/setting_%s/%d.RData", s, block))
      R = readRDS(sprintf("data/LD/%d/LD_UKBB.RData", block))
      gLDSC_var_sub = gLDSC_var[data_full$rsid,]
  
      for(k in 1:5)
      {
        dir.create(sprintf("results/setting_%s/%s/%s/%s",
                           s, Dir, suffixes[k], block), recursive = T)
        
        data = data_full[,c(1:6,14:20,20+k)]
        n = c(N1, N2_seq, N2_seq)[k]
        z = data[,14]
        z = sqrt(n - 1) * z / sqrt(n - z^2)
        a = gLDSC_var_sub[,k] / mean(gLDSC_var_sub[,k]) * 0.005
        
        h2D2 = Createh2D2Object(z,
                                R[[idx[k]]],
                                c(N1, N2_seq, N2_seq)[k],
                                data$rsid,
                                in_sample_LD = T,
                                a = a)
  
        tic = Sys.time()
        h2D2 = h2D2_MCMC(h2D2, mcmc_n = 5500, burn_in = 500, thin = 2, get_CS = F)
        while(max(h2D2@mcmc_samples$PSRF_beta) > 1.2)
        {
          h2D2 = h2D2_MCMC(h2D2, mcmc_n = 1000, burn_in = 1000, thin = 2, get_CS = F)
        }
        h2D2@CS = h2D2_CS(h2D2)
        toc = Sys.time()
  
        h2D2@mcmc_samples$beta = NULL
        h2D2@mcmc_samples$tau = NULL
        h2D2@mcmc_samples$psi = NULL
        saveRDS(h2D2, file = sprintf(
          "results/setting_%s/%s/%s/%s/h2D2.RData",
          s, Dir, suffixes[k], block))
  
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
          "results/setting_%s/%s/%s/%s/data.RData",
          s, Dir, suffixes[k], block))
        
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
                           data[h2D2@CS$sets[[l]], 9+c(1,2,2,2,2)[k]]),
                         "size" = length(h2D2@CS$sets[[l]]),
                         "min.abs.corr" = h2D2@CS$purity$min.abs.corr[l]))
          }
        }
        
        write_delim(result_CS95, sprintf(
          "results/setting_%s/%s/%s/%s/result_CS95.txt",
          s, Dir, suffixes[k], block), delim = '\t')
        
        result = data.frame(
          "block" = block,
          "max_PSRF" = max(h2D2@mcmc_samples$PSRF_beta),
          "mcmc_n" = h2D2@mcmc_samples$n_samples + h2D2@mcmc_samples$n_burnin,
          "time" = difftime(toc, tic, units = "secs"))
        
        if(k > 1)
        {
          data_list = list()
          data_list[[1]] = readRDS(sprintf(
            "results/setting_%s/%s/%s/%s/data.RData",
            s, Dir, suffixes[1], block))
          data_list[[2]] = data
          
          h2D2_list = list()
          h2D2_list[[1]] = readRDS(sprintf(
            "results/setting_%s/%s/%s/%s/h2D2.RData",
            s, Dir, suffixes[1], block))
          h2D2_list[[2]] = h2D2
          
          result_CS95_cross = data.frame("block" = block,
                                         "coverage" = 0,
                                         "size" = 0,
                                         "min.abs.corr" = 0)[-1,]
          result_CS95_cross %<>% rbind(
            foreach(m = 1:2, .combine = "rbind") %do%
            {
              if(length(na.omit(data_list[[m]]$CS)) > 0)
              {
                foreach(l = unique(na.omit(data_list[[m]]$CS)),
                        .combine = "rbind") %do%
                  {
                    set = which(data_list[[m]]$CS == l)
                    purity = h2D2_list[[m]]@CS$purity$min.abs.corr[l]
                    data.frame("block" = block,
                               "coverage" = sum(data[set,10] | data[set,11]),
                               "size" = length(set),
                               "min.abs.corr" = purity)
                  }
              }
            }
          )
          
          write_delim(result_CS95_cross, sprintf(
            "results/setting_%s/%s/%s/%s/result_CS95_cross.txt",
            s, Dir, suffixes[k], block), delim = '\t')
          
          result_CS95_shared = data.frame("block" = block,
                                          "coverage" = 0,
                                          "size" = 0,
                                          "min.abs.corr" = 0)[-1,]
          if(length(h2D2_list[[1]]@CS$sets) > 0 & length(h2D2_list[[2]]@CS$sets) > 0)
          {
            for(l1 in 1:length(h2D2_list[[1]]@CS$sets))
            {
              set1 = h2D2_list[[1]]@CS$sets[[l1]]
              for(l2 in 1:length(h2D2_list[[2]]@CS$sets))
              {
                set2 = h2D2_list[[2]]@CS$sets[[l2]]
                if(length(intersect(set1,set2)) > 0)
                {
                  result_CS95_shared %<>% rbind(data.frame(
                    "block" = block,
                    "coverage" = sum(data[intersect(set1,set2),10] &
                                       data[intersect(set1,set2),11]),
                    "size" = length(union(set1,set2)),
                    "min.abs.corr" = min(h2D2_list[[1]]@CS$purity$min.abs.corr[l1],
                                         h2D2_list[[2]]@CS$purity$min.abs.corr[l2]))
                  )
                }
              }
            }
          }
          
          write_delim(result_CS95_shared, sprintf(
            "results/setting_%s/%s/%s/%s/result_CS95_shared.txt",
            s, Dir, suffixes[k], block), delim = '\t')
          
          result$CS95_causal = sum((!is.na(data_list[[1]]$CS) | 
                                      !is.na(data_list[[2]]$CS)) &
                                     (data$causal_1 | data$causal_2))
          result$CS95_causal_1 = sum(!is.na(data_list[[1]]$CS) & data$causal_1)
          result$CS95_causal_2 = sum(!is.na(data_list[[2]]$CS) & data$causal_2)
          result$CS95_causal_0 = sum((!is.na(data_list[[1]]$CS) & 
                                        !is.na(data_list[[2]]$CS)) &
                                       (data$causal_1 & data$causal_2))
        }
  
        write_delim(result, sprintf(
          "results/setting_%s/%s/%s/%s/result.txt",
          s, Dir, suffixes[k], block), delim = '\t')
      }
      
      NULL
    }
    
    for(k in 1:5)
    {
      results = foreach(block = select_block, .combine = "rbind") %dopar%
      {
        read.delim(sprintf(
          "results/setting_%s/%s/%s/%s/result.txt",
          s, Dir, suffixes[k], block))
      }
      write_delim(results, sprintf(
        "results/setting_%s/%s/%s/results.txt",
        s, Dir, suffixes[k]), delim = '\t')
  
      results_CS95 = foreach(block = select_block, .combine = "rbind") %dopar%
      {
        read.delim(sprintf(
          "results/setting_%s/%s/%s/%s/result_CS95.txt",
          s, Dir, suffixes[k], block))
      }
      write_delim(results_CS95, sprintf(
        "results/setting_%s/%s/%s/results_CS95.txt",
        s, Dir, suffixes[k]), delim = '\t')
      
      if(k > 1)
      {
        results_CS95_cross = foreach(block = select_block, .combine = "rbind") %dopar%
        {
          read.delim(sprintf(
            "results/setting_%s/%s/%s/%s/result_CS95_cross.txt",
            s, Dir, suffixes[k], block))
        }
        write_delim(results_CS95_cross, sprintf(
          "results/setting_%s/%s/%s/results_CS95_cross.txt",
          s, Dir, suffixes[k]), delim = '\t')
        
        results_CS95_shared = foreach(block = select_block, .combine = "rbind") %dopar%
        {
          read.delim(sprintf(
            "results/setting_%s/%s/%s/%s/result_CS95_shared.txt",
            s, Dir, suffixes[k], block))
        }
        write_delim(results_CS95_shared, sprintf(
          "results/setting_%s/%s/%s/results_CS95_shared.txt",
          s, Dir, suffixes[k]), delim = '\t')
      }
    }
    
    data_1 = foreach(block = select_block, .combine = "rbind") %dopar%
    {
      readRDS(sprintf(
        "results/setting_%s/%s/%s/%s/data.RData",
        s, Dir, suffixes[1], block))
    }
    
    for(m in 1:4)
    {
      data_2 = foreach(block = select_block, .combine = "rbind") %dopar%
      {
        readRDS(sprintf(
          "results/setting_%s/%s/%s/%s/data.RData",
          s, Dir, suffixes[m+1], block))
      }
      
      # Calibration
      
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
        "results/setting_%s/%s/%s/calibration.RData",
        s, Dir, suffixes[m+1]))
    }
  }
}
