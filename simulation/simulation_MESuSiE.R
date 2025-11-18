options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(scales)
library(MESuSiE)

setwd("simulation")

select_block_info = read.delim("data/select_block_info.txt")
select_block_info %<>% arrange(num_variants)
select_block_info[101:200,] %<>% arrange(-num_variants)
select_block = sort(select_block_info$block)

pids = c("EUR", "EAS")

n_causal = data.frame("n" = c(5,3,1),
                      "n1" = c(0,1,2),
                      "n2" = c(0,1,2))

N1 = 2e5
N2_seq = c(2e4,2e5)

group_bound = data.frame("l" = c(0,0.1,0.5,0.9),
                         "u" = c(0.1,0.5,0.9,1.1))

registerDoParallel(10)

Dir = "MESuSiE"

for(s in 1:3)
{
  foreach(block = select_block_info$block, .combine = "c", .inorder = F) %dopar%
  {
    dir.create(sprintf("results/setting_%s/%s/%s",
                       s, Dir, block), recursive = T)
    
    data_full = readRDS(sprintf("data/setting_%s/%d.RData", s, block))
    R = readRDS(sprintf("data/LD/%d.RData", block))
    
    for(k in 1:2)
    {
      R[[k]] = as.matrix(R[[k]])
      diag(R[[k]]) = 1
    }
    
    for(m in 1:2)
    {
      data = data_full[,c(1:13,13+m)]
      Z = data[,13:14]
      N = c(N1,N2_seq[m])
      
      sumstat = foreach(k = 1:2) %do%
      {
        tmp = data.frame("SNP" = data$rsid)
        tmp$Beta = Z[,k] / sqrt(N[k])
        mutate(tmp, Se = 1 / sqrt(N[k]),
               Z = Beta / Se,
               N = N[k])
      }
      rm(tmp)
      names(R) = names(sumstat) = 1:2
      
      tic = Sys.time()
      MESuSiE_res = meSuSie_core(R, sumstat, L = 5)
      toc = Sys.time()
      
      saveRDS(MESuSiE_res, file = sprintf(
        "results/setting_%s/%s/%s/MESuSiE_res_N2-%d.RData",
        s, Dir, block, N2_seq[m]))
      
      data$PIP = MESuSiE_res$pip
      data$PIP_0 = MESuSiE_res$pip_config[,3]
      data$PIP_1 = rowSums(MESuSiE_res$pip_config[,c(1,3)])
      data$PIP_2 = rowSums(MESuSiE_res$pip_config[,c(2,3)])
      data$CS_2 = data$CS_1 = NA
      
      result_CS95 = foreach(k = 1:2) %do%
      {
        data.frame("block" = block,
                   "coverage" = 0,
                   "size" = 0,
                   "min.abs.corr" = 0)[-1,]
      }
      
      result_CS95_cross = data.frame("block" = block,
                                     "coverage" = 0,
                                     "size" = 0,
                                     "min.abs.corr" = 0)[-1,]
      
      result_CS95_shared = data.frame("block" = block,
                                      "coverage" = 0,
                                      "size" = 0,
                                      "min.abs.corr" = 0)[-1,]
      
      if(length(MESuSiE_res$cs$cs) > 0)
      {
        for(l in 1:length(MESuSiE_res$cs$cs))
        {
          set = MESuSiE_res$cs$cs[[l]]
          
          purities = foreach(k = 1:2, .combine = "c") %do%
            {
              if(length(set) == 1)
              {
                1
              } else {
                min(abs(R[[k]][set, set]))
              }
            }
          
          if(MESuSiE_res$cs$cs_category[l] == "1") {
            result_CS95_cross %<>% rbind(data.frame(
              "block" = block,
              "coverage" = sum(data[set,9] | data[set,10]),
              "size" = length(set),
              "min.abs.corr" = purities[1]))
            
            result_CS95[[1]] %<>% rbind(data.frame(
              "block" = block,
              "coverage" = sum(data[set,9]),
              "size" = length(set),
              "min.abs.corr" = purities[1]))
            
            data$CS_1[set] = l
          } else if(MESuSiE_res$cs$cs_category[l] == "2") {
            result_CS95_cross %<>% rbind(data.frame(
              "block" = block,
              "coverage" = sum(data[set,9] | data[set,10]),
              "size" = length(set),
              "min.abs.corr" = purities[2]))
            
            result_CS95[[2]] %<>% rbind(data.frame(
              "block" = block,
              "coverage" = sum(data[set,10]),
              "size" = length(set),
              "min.abs.corr" = purities[2]))
            
            data$CS_2[set] = l
          } else {
            result_CS95_cross %<>% rbind(
              foreach(k = 1:2, .combine = "rbind") %do%
                {
                  data.frame("block" = block,
                             "coverage" = sum(data[set,9] | data[set,10]),
                             "size" = length(set),
                             "min.abs.corr" = purities[k])
                }
            )
            
            for(k in 1:2)
            {
              result_CS95[[k]] %<>% rbind(data.frame(
                "block" = block,
                "coverage" = sum(data[set,8+k]),
                "size" = length(set),
                "min.abs.corr" = purities[k]))
            }
            
            result_CS95_shared %<>% rbind(data.frame(
              "block" = block,
              "coverage" = sum(data[set,9] & data[set,10]),
              "size" = length(set),
              "min.abs.corr" = purities[k]))
            
            data$CS_1[set] = data$CS_2[set] = l
          }
        }
      }
      
      saveRDS(data, file = sprintf(
        "results/setting_%s/%s/%s/data_N2-%d.RData",
        s, Dir, block, N2_seq[m]))
      
      for(k in 1:2)
      {
        write_delim(result_CS95[[k]], sprintf(
          "results/setting_%s/%s/%s/result_CS95_%d_N2-%d.txt",
          s, Dir, block, k, N2_seq[m]), delim = '\t')
      }
      
      write_delim(result_CS95_cross, sprintf(
        "results/setting_%s/%s/%s/result_CS95_N2-%d_cross.txt",
        s, Dir, block, N2_seq[m]), delim = '\t')
      
      write_delim(result_CS95_shared, sprintf(
        "results/setting_%s/%s/%s/result_CS95_N2-%d_shared.txt",
        s, Dir, block, N2_seq[m]), delim = '\t')
      
      result = data.frame(
        "block" = block,
        "time" = difftime(toc, tic, units = "secs"),
        "CS95_causal" = sum((!is.na(data$CS_1) | !is.na(data$CS_2)) &
                              (data$causal_1 | data$causal_2)),
        "CS95_causal_1" = sum(!is.na(data$CS_1) & data$causal_1),
        "CS95_causal_2" = sum(!is.na(data$CS_2) & data$causal_2),
        "CS95_causal_0" = sum((!is.na(data$CS_1) & !is.na(data$CS_2)) &
                                (data$causal_1 & data$causal_2)))
      
      write_delim(result, sprintf(
        "results/setting_%s/%s/%s/result_N2-%d.txt",
        s, Dir, block, N2_seq[m]), delim = '\t')
      
      NULL
    }
    
    NULL
  }
  
  for(m in 1:2)
  {
    results = foreach(block = select_block, .combine = "rbind") %dopar%
    {
      read.delim(sprintf(
        "results/setting_%s/%s/%s/result_N2-%d.txt",
        s, Dir, block, N2_seq[m]))
    }
    write_delim(results, sprintf(
      "results/setting_%s/%s/results_N2-%d.txt",
      s, Dir, N2_seq[m]), delim = '\t')
    
    results_CS95 = list()
    
    for(k in 1:2)
    {
      results_CS95[[k]] = foreach(block = select_block, .combine = "rbind") %dopar%
      {
        read.delim(sprintf(
          "results/setting_%s/%s/%s/result_CS95_%d_N2-%d.txt",
          s, Dir, block, k, N2_seq[m]))
      }
      write_delim(results_CS95[[k]], sprintf(
        "results/setting_%s/%s/results_CS95_%d_N2-%d.txt",
        s, Dir, k, N2_seq[m]), delim = '\t')
    }
    
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
    
    # Calibration
    
    data = foreach(block = select_block, .combine = "rbind") %dopar%
    {
      readRDS(sprintf(
        "results/setting_%s/%s/%s/data_N2-%d.RData",
        s, Dir, block, N2_seq[m]))
    }
    
    calibration = list()
    
    PIP = data$PIP
    calibration[["cross"]] = foreach(g = 1:4, .combine = "rbind") %dopar%
    {
      idx = which(PIP >= group_bound$l[g] & PIP < group_bound$u[g])
      n = length(idx)
      
      data.frame("group" = g,
                 "n" = n,
                 "Expected" = sum(PIP[idx]) / n,
                 "Prop" = sum((data$causal_1 | data$causal_2)[idx]) / n)
    }
    
    PIP = data$PIP_0
    calibration[["shared"]] = foreach(g = 1:4, .combine = "rbind") %dopar%
    {
      idx = which(PIP >= group_bound$l[g] & PIP < group_bound$u[g])
      n = length(idx)
      
      data.frame("group" = g,
                 "n" = n,
                 "Expected" = sum(PIP[idx]) / n,
                 "Prop" = sum((data$causal_1 & data$causal_2)[idx]) / n)
    }
    
    PIP = data$PIP_1
    calibration[["pop1"]] = foreach(g = 1:4, .combine = "rbind") %dopar%
    {
      idx = which(PIP >= group_bound$l[g] & PIP < group_bound$u[g])
      n = length(idx)
      
      data.frame("group" = g,
                 "n" = n,
                 "Expected" = sum(PIP[idx]) / n,
                 "Prop" = sum(data$causal_1[idx]) / n)
    }
    
    PIP = data$PIP_2
    calibration[["pop2"]] = foreach(g = 1:4, .combine = "rbind") %dopar%
    {
      idx = which(PIP >= group_bound$l[g] & PIP < group_bound$u[g])
      n = length(idx)
      
      data.frame("group" = g,
                 "n" = n,
                 "Expected" = sum(PIP[idx]) / n,
                 "Prop" = sum(data$causal_2[idx]) / n)
    }
    
    saveRDS(calibration, sprintf(
      "results/setting_%s/%s/calibration_N2-%d.RData",
      s, Dir, N2_seq[m]))
  }
}