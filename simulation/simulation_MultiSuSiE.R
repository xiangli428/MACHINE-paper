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

pids = c("EUR", "AFR", "EAS")

n_causal = data.frame("n" = c(5,3,1),
                      "n1" = c(0,1,2),
                      "n2" = c(0,1,2))

N1 = 2e5
N2_seq = c(2e4,2e5)

suffixes = c("N1-200000", sprintf("N2-%d", N2_seq), sprintf("N3-%d", N2_seq))

group_bound = data.frame("l" = c(0,0.1,0.5,0.9),
                         "u" = c(0.1,0.5,0.9,1.1))

registerDoParallel(10)

Dir = "MultiSuSiE"

for(s in 1:3)
{
  foreach(block = select_block_info$block, .combine = "c", .inorder = F) %dopar%
  {
    dir.create(sprintf("results/setting_%s/%s/%s", s, Dir, block), recursive = T)
    
    data_full = readRDS(sprintf("data/setting_%s/%d.RData", s, block))
    R = readRDS(sprintf("data/LD/%d/LD_UKBB.RData", block))
    
    for(k in 1:3)
    {
      R[[k]] = as.matrix(R[[k]])
      diag(R[[k]]) = 1
      write.table(R[[k]], sprintf(
        "results/setting_%s/%s/%s/LD_%s.txt", s, Dir, block, pids[k]),
        sep = '\t', row.names = F, col.names = F, quote = F)
    }
    
    data_sub = data_full[,c(5,21:25)]
    for(k in 1:5)
    {
      n = c(N1, N2_seq, N2_seq)[k]
      data_sub[,k+1] = sqrt(n - 1) * data_sub[,1+k] / sqrt(n - data_sub[,1+k]^2)
    }
    
    for(p in 2:3)
    {
      for(m in 1:2)
      {
        suffix = sprintf("N%d-%d", p, N2_seq[m])
        dir.create(sprintf("results/setting_%s/%s/%s/%s",
                           s, Dir, suffix, block), recursive = T)
        
        write_delim(data_sub[,c(1,2,2*p+m-2)], sprintf(
          "results/setting_%s/%s/%s/%s/data.txt", s, Dir, suffix, block),
          delim = '\t')
        
        tic = Sys.time()
        system(sprintf(
          "MultiSuSiE/bin/python simulation_MultiSuSiE.py %s %s %s %d %d",
          sprintf("results/setting_%s/%s/%s/%s/", s, Dir, suffix, block),
          sprintf("results/setting_%s/%s/%s/LD_%s.txt", s, Dir, block, pids[1]),
          sprintf("results/setting_%s/%s/%s/LD_%s.txt", s, Dir, block, pids[p]),
          2e5, N2_seq[m]))
        toc = Sys.time()
        
        data = read.delim(sprintf(
          "results/setting_%s/%s/%s/%s/data.txt", s, Dir, suffix, block))
        data = cbind(data_full[,c(1:6,14:20)], data[,2:5])
        data$CS[data$CS == 0] = NA
        
        saveRDS(data, file = sprintf(
          "results/setting_%s/%s/%s/%s/data.RData",
          s, Dir, suffix, block))
        
        result = data.frame(
          "time" = difftime(toc, tic, units = "secs"),
          "CS95_causal" = sum(!is.na(data$CS) & (data$causal_1 | data$causal_2)),
          "CS95_causal_1" = sum(!is.na(data$CS) & data$causal_1),
          "CS95_causal_2" = sum(!is.na(data$CS) & data$causal_2),
          "CS95_causal_0" = sum(!is.na(data$CS) & (data$causal_1 & data$causal_2)))
        
        write_delim(result, sprintf(
          "results/setting_%s/%s/%s/%s/result.txt",
          s, Dir, suffix, block), delim = '\t')
        
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
                min(abs(R[[c(1,p)[k]]][set, set]))
              }
            }
            
            for(k in 1:2)
            {
              result_CS95_cross %<>% rbind(data.frame(
                "block" = block,
                "coverage" = sum(data[set,10] | data[set,11]),
                "size" = length(set),
                "min.abs.corr" = purities[k]))
              
              result_CS95[[k]] %<>% rbind(data.frame(
                "block" = block,
                "coverage" = sum(data[set,9+k]),
                "size" = length(set),
                "min.abs.corr" = purities[k]))
              
              result_CS95_shared %<>% rbind(data.frame(
                "block" = block,
                "coverage" = sum(data[set,10] & data[set,11]),
                "size" = length(set),
                "min.abs.corr" = max(purities)))
            }
          }
        }
        
        for(k in 1:2)
        {
          write_delim(result_CS95[[k]], sprintf(
            "results/setting_%s/%s/%s/%s/result_CS95_%d.txt",
            s, Dir, suffix, block, k), delim = '\t')
        }
        
        write_delim(result_CS95_cross, sprintf(
          "results/setting_%s/%s/%s/%s/result_CS95_cross.txt",
          s, Dir, suffix, block), delim = '\t')
        
        write_delim(result_CS95_shared, sprintf(
          "results/setting_%s/%s/%s/%s/result_CS95_shared.txt",
          s, Dir, suffix, block), delim = '\t')
      }
    }
    
    system(sprintf("rm -r results/setting_%s/%s/%s", s, Dir, block))
  }
  
  for(suffix in c("N2-20000","N2-200000","N3-20000","N3-200000"))
  {
    results = foreach(block = select_block, .combine = "rbind") %dopar%
    {
      read.delim(sprintf(
        "results/setting_%s/%s/%s/%s/result.txt",
        s, Dir, suffix, block))
    }
    write_delim(results, sprintf(
      "results/setting_%s/%s/%s/results.txt",
      s, Dir, suffix), delim = '\t')

    results_CS95 = list()

    for(k in 1:2)
    {
      results_CS95[[k]] = foreach(block = select_block, .combine = "rbind") %dopar%
      {
        read.delim(sprintf(
          "results/setting_%s/%s/%s/%s/result_CS95_%d.txt",
          s, Dir, suffix, block, k))
      }
      write_delim(results_CS95[[k]], sprintf(
        "results/setting_%s/%s/%s/results_CS95_%d.txt",
        s, Dir, suffix, k), delim = '\t')
    }

    results_CS95_cross = foreach(block = select_block, .combine = "rbind") %dopar%
    {
      read.delim(sprintf(
        "results/setting_%s/%s/%s/%s/result_CS95_cross.txt",
        s, Dir, suffix, block))
    }
    write_delim(results_CS95_cross, sprintf(
      "results/setting_%s/%s/%s/results_CS95_cross.txt",
      s, Dir, suffix), delim = '\t')

    results_CS95_shared = foreach(block = select_block, .combine = "rbind") %dopar%
    {
      read.delim(sprintf(
        "results/setting_%s/%s/%s/%s/result_CS95_shared.txt",
        s, Dir, suffix, block))
    }
    write_delim(results_CS95_shared, sprintf(
      "results/setting_%s/%s/%s/results_CS95_shared.txt",
      s, Dir, suffix), delim = '\t')

    data = foreach(block = select_block, .combine = "rbind") %dopar%
    {
      readRDS(sprintf(
        "results/setting_%s/%s/%s/%s/data.RData",
        s, Dir, suffix, block))
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
      "results/setting_%s/%s/%s/calibration.RData",
      s, Dir, suffix))
  }
}
