options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)
library(susieR)

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

Dir = "SuSiE"

for(s in 1:3)
{
  foreach(block = select_block_info$block, .combine = "c", .inorder = F) %dopar%
  {
    data_full = readRDS(sprintf("data/setting_%s/%d.RData", s, block))
    R = readRDS(sprintf("data/LD/%d/LD_UKBB.RData", block))
    for(k in 1:3)
    {
      R[[k]] %<>% as.matrix()
      diag(R[[k]]) = 1
    }

    for(k in 1:5)
    {
      dir.create(sprintf("results/setting_%s/%s/%s/%s",
                         s, Dir, suffixes[k], block), recursive = T)
      
      data = data_full[,c(1:6,14:20,20+k)]
      n = c(N1, N2_seq, N2_seq)[k]
      
      z = data[,14]
      z = sqrt(n - 1) * z / sqrt(n - z^2)

      tic = Sys.time()
      fit_susie = susie_rss(z = z,
                            R = R[[idx[k]]],
                            var_y = 1,
                            n = n,
                            L = 5,
                            estimate_residual_variance = T,
                            refine = T)
      toc = Sys.time()

      saveRDS(fit_susie, file = sprintf(
        "results/setting_%s/%s/%s/%s/fit_susie.RData",
        s, Dir, suffixes[k], block))

      data$PIP = fit_susie$pip

      data$CS = NA
      if(length(fit_susie$sets$cs) > 0)
      {
        for(l in 1:length(fit_susie$sets$cs))
        {
          data$CS[fit_susie$sets$cs[[l]]] = names(fit_susie$sets$cs)[l]
        }
      }

      saveRDS(data, file = sprintf(
        "results/setting_%s/%s/%s/%s/data.RData",
        s, Dir, suffixes[k], block))
      
      result_CS95 = data.frame("block" = block,
                               "coverage" = 0,
                               "size" = 0,
                               "min.abs.corr" = 0)[-1,]
      if(length(fit_susie$sets$cs) > 0)
      {
        for(l in 1:length(fit_susie$sets$cs))
        {
          result_CS95 %<>% rbind(
            data.frame("block" = block,
                       "coverage" = sum(
                         data[fit_susie$sets$cs[[l]], 9+c(1,2,2,2,2)[k]]),
                       "size" = length(fit_susie$sets$cs[[l]]),
                       "min.abs.corr" = fit_susie$sets$purity$min.abs.corr[l]))
        }
      }
      
      write_delim(result_CS95, sprintf(
        "results/setting_%s/%s/%s/%s/result_CS95.txt",
        s, Dir, suffixes[k], block), delim = '\t')
      
      result = data.frame(
        "block" = block,
        "time" = difftime(toc, tic, units = "secs"))
      
      if(k > 1)
      {
        data_list = list()
        data_list[[1]] = readRDS(sprintf(
          "results/setting_%s/%s/%s/%s/data.RData",
          s, Dir, suffixes[1], block))
        data_list[[2]] = data
        
        fit_susie_list = list()
        fit_susie_list[[1]] = readRDS(sprintf(
          "results/setting_%s/%s/%s/%s/fit_susie.RData",
          s, Dir, suffixes[1], block))
        fit_susie_list[[2]] = fit_susie
        
        result_CS95_cross = data.frame("block" = block,
                                       "coverage" = 0,
                                       "size" = 0,
                                       "min.abs.corr" = 0)[-1,]
        result_CS95_cross %<>% rbind(
          foreach(m = 1:2, .combine = "rbind") %do%
          {
            if(length(fit_susie_list[[m]]$sets$cs) > 0)
            {
              foreach(l = 1:length(fit_susie_list[[m]]$sets$cs), .combine = "rbind") %do%
              {
                set = fit_susie_list[[m]]$sets$cs[[l]]
                
                data.frame("block" = block,
                           "coverage" = sum(data[set,10] | data[set,11]),
                           "size" = length(set),
                           "min.abs.corr" = fit_susie_list[[m]]$sets$purity$min.abs.corr[l])
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
        if(length(fit_susie_list[[1]]$sets$cs) > 0 & length(fit_susie_list[[2]]$sets$cs) > 0)
        {
          for(l1 in 1:length(fit_susie_list[[1]]$sets$cs))
          {
            set1 = fit_susie_list[[1]]$sets$cs[[l1]]
            for(l2 in 1:length(fit_susie_list[[2]]$sets$cs))
            {
              set2 = fit_susie_list[[2]]$sets$cs[[l2]]
              if(length(intersect(set1,set2)) > 0)
              {
                result_CS95_shared %<>% rbind(data.frame(
                  "block" = block,
                  "coverage" = sum(data[intersect(set1,set2),10] &
                                     data[intersect(set1,set2),11]),
                  "size" = length(union(set1,set2)),
                  "min.abs.corr" = min(fit_susie_list[[1]]$sets$purity$min.abs.corr[l1],
                                       fit_susie_list[[2]]$sets$purity$min.abs.corr[l2]))
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

    PIP = pmax(data_1$PIP, data_2$PIP)
    calibration[["cross"]] = foreach(g = 1:4, .combine = "rbind") %dopar%
    {
      idx = which(PIP >= group_bound$l[g] & PIP < group_bound$u[g])
      n = length(idx)

      data.frame("group" = g,
                 "n" = n,
                 "Expected" = sum(PIP[idx]) / n,
                 "Prop" = sum((data_1$causal_1 | data_2$causal_2)[idx]) / n)
    }

    PIP = pmin(data_1$PIP, data_2$PIP)
    calibration[["shared"]] = foreach(g = 1:4, .combine = "rbind") %dopar%
    {
      idx = which(PIP >= group_bound$l[g] & PIP < group_bound$u[g])
      n = length(idx)

      data.frame("group" = g,
                 "n" = n,
                 "Expected" = sum(PIP[idx]) / n,
                 "Prop" = sum((data_1$causal_1 & data_2$causal_2)[idx]) / n)
    }

    PIP = data_1$PIP
    calibration[["pop1"]] = foreach(g = 1:4, .combine = "rbind") %dopar%
    {
      idx = which(PIP >= group_bound$l[g] & PIP < group_bound$u[g])
      n = length(idx)

      data.frame("group" = g,
                 "n" = n,
                 "Expected" = sum(PIP[idx]) / n,
                 "Prop" = sum(data_1$causal_1[idx]) / n)
    }

    PIP = data_2$PIP
    calibration[["pop2"]] = foreach(g = 1:4, .combine = "rbind") %dopar%
    {
      idx = which(PIP >= group_bound$l[g] & PIP < group_bound$u[g])
      n = length(idx)

      data.frame("group" = g,
                 "n" = n,
                 "Expected" = sum(PIP[idx]) / n,
                 "Prop" = sum(data_2$causal_2[idx]) / n)
    }

    saveRDS(calibration, sprintf(
      "results/setting_%s/%s/%s/calibration.RData",
      s, Dir, suffixes[m+1]))
  }
}