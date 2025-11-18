options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)

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

Dir = "MultiSuSiE"

for(s in 1:3)
{
  foreach(block = select_block_info$block, .combine = "c", .inorder = F) %dopar%
  {
    dir.create(sprintf("results/setting_%s/%s/%s",
                       s, Dir, block), recursive = T)
    
    data_full = readRDS(sprintf("data/setting_%s/%d.RData", s, block))
    R = readRDS(sprintf("data/LD/%d.RData", block))
    
    data = data_full[,c(1,13:17)]
    
    for(m in 1:5)
    {
      n = c(N1, N2_seq)[m]
      data[,1+m] = sqrt(n - 1) * data[,1+m] / sqrt(n - data[,1+m]^2)
    }
    write_delim(data, sprintf("results/setting_%s/%s/%s/data.txt",
                              s, Dir, block), delim = '\t')
    
    for(k in 1:2)
    {
      R[[k]] = as.matrix(R[[k]])
      diag(R[[k]]) = 1
      write.table(R[[k]], sprintf(
        "results/setting_%s/%s/%s/LD_%d.txt", s, Dir, block, k),
        sep = '\t', row.names = F, col.names = F, quote = F)
    }
    
    system(sprintf(
      "/home/r14user2/anaconda3/envs/MultiSuSiE/bin/python simulation_MultiSuSiE.py %s",
      sprintf("results/setting_%s/%s/%s/", s, Dir, block)))
    
    runtime = read.delim(sprintf("results/setting_%s/%s/%s/runtime.txt",
                                 s, Dir, block), header = F)[,1]
    
    for(m in 1:2)
    {
      system(sprintf(
        "/home/r14user2/anaconda3/envs/MultiSuSiE/bin/python simulation_MultiSuSiE.py %s %s",
        sprintf("results/setting_%s/%s/%s/", s, Dir, block), m-1))
      
      data = read.delim(sprintf("results/setting_%s/%s/%s/data_N2-%d.txt", 
                                s, Dir, block, N2_seq[m]))
      data = cbind(data_full[,c(1:13,13+m)], data[,4:5])
      data$CS[data$CS == 0] = NA
      
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
      
      if(sum(!is.na(data$CS)) > 0)
      {
        for(l in sort(unique(na.omit(data$CS))))
        {
          set = which(data$CS == l)
          
          purities = foreach(k = 1:2, .combine = "c") %do%
          {
            if(length(set) == 1)
            {
              1
            } else {
              min(abs(R[[k]][set, set]))
            }
          }
          
          for(k in 1:2)
          {
            result_CS95_cross %<>% rbind(data.frame(
              "block" = block,
              "coverage" = sum(data[set,9] | data[set,10]),
              "size" = length(set),
              "min.abs.corr" = purities[k]))
            
            result_CS95[[k]] %<>% rbind(data.frame(
              "block" = block,
              "coverage" = sum(data[set,8+k]),
              "size" = length(set),
              "min.abs.corr" = purities[k]))
            
            result_CS95_shared %<>% rbind(data.frame(
              "block" = block,
              "coverage" = sum(data[set,9] & data[set,10]),
              "size" = length(set),
              "min.abs.corr" = max(purities)))
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
        "CS95_causal" = sum(!is.na(data$CS) & (data$causal_1 | data$causal_2)),
        "CS95_causal_1" = sum(!is.na(data$CS) & data$causal_1),
        "CS95_causal_2" = sum(!is.na(data$CS) & data$causal_2),
        "CS95_causal_0" = sum(!is.na(data$CS) & (data$causal_1 & data$causal_2)))
      
      write_delim(result, sprintf(
        "results/setting_%s/%s/%s/result_N2-%d.txt",
        s, Dir, block, N2_seq[m]), delim = '\t')
    }
    
    system(sprintf("rm results/setting_%s/%s/%s/data.txt",
                   s, Dir, block))
    for(k in 1:2)
    {
      system(sprintf("rm results/setting_%s/%s/%s/LD_%d.txt", 
                     s, Dir, block, k))
    }
    for(m in c(1,4))
    {
      system(sprintf("rm results/setting_%s/%s/%s/data_N2-%d.txt", 
                     s, Dir, block, N2_seq[m]))
    }
    system(sprintf("rm results/setting_%s/%s/%s/runtime.txt",
                   s, Dir, block))
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
    
    calibration[["shared"]] = foreach(g = 1:4, .combine = "rbind") %dopar%
    {
      idx = which(PIP >= group_bound$l[g] & PIP < group_bound$u[g])
      n = length(idx)
      
      data.frame("group" = g,
                 "n" = n,
                 "Expected" = sum(PIP[idx]) / n,
                 "Prop" = sum((data$causal_1 & data$causal_2)[idx]) / n)
    }
    
    calibration[["pop1"]] = foreach(g = 1:4, .combine = "rbind") %dopar%
    {
      idx = which(PIP >= group_bound$l[g] & PIP < group_bound$u[g])
      n = length(idx)
      
      data.frame("group" = g,
                 "n" = n,
                 "Expected" = sum(PIP[idx]) / n,
                 "Prop" = sum(data$causal_1[idx]) / n)
    }
    
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