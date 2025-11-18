options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(tidyr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)
library(PRROC)

setwd("~/Documents/GWAS/project_4/simulation")

select_block_info = read.delim("data/select_block_info.txt")
select_block_info %<>% arrange(num_variants)
select_block_info[101:200,] %<>% arrange(-num_variants)
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

registerDoParallel(20)

Dir = "RSparsePro_1kg"

for(s in 1:3)
{
  foreach(block = select_block_info$block, .combine = "c", .inorder = F) %dopar%
  {
    data_full = readRDS(sprintf("data/setting_%s/%d.RData", s, block))
    R = readRDS(sprintf("data/LD_1kg/%d.RData", block))
    M = nrow(data_full)

    for(k in 1:2)
    {
      R[[k]] %<>% as.matrix()
      diag(R[[k]]) = 1
    }

    for(k in 1:3)
    {
      dir.create(sprintf("results/setting_%s/%s/%s/%s",
                         s, Dir, block, suffix[k]), recursive = T)

      data = data_full[,c(1:12,12+k)]
      n = c(N1, N2_seq)[k]

      z = data[,13]

      write_delim(data.frame("RSID" = data$rsid, "Z" = z),
                  sprintf("results/setting_%s/%s/%s/%s/z.txt",
                          s, Dir, block, suffix[k]), delim = '\t')
      write.table(R[[idx[k]]], sprintf(
        "results/setting_%s/%s/%s/%s/ld.txt", s, Dir, block, suffix[k]),
        sep = '\t', row.names = F, col.names = F, quote = F)

      tic = Sys.time()
      system(paste("/home/r9user9/anaconda3/bin/python",
                   "~/Documents/GWAS/packages/RSparsePro_LD/src/rsparsepro_ld.py",
                   "--z", sprintf("results/setting_%s/%s/%s/%s/z.txt",
                                  s, Dir, block, suffix[k]),
                   "--ld", sprintf("results/setting_%s/%s/%s/%s/ld.txt",
                                   s, Dir, block, suffix[k]),
                   "--save", sprintf("results/setting_%s/%s/%s/%s/result",
                                     s, Dir, block, suffix[k]),
                   "--K 5 --minldthres 0.5"))
      toc = Sys.time()

      system(sprintf("rm results/setting_%s/%s/%s/%s/z.txt",
                     s, Dir, block, suffix[k]))
      system(sprintf("rm results/setting_%s/%s/%s/%s/ld.txt",
                     s, Dir, block, suffix[k]))
      
      vare = system(sprintf(
        "cat results/setting_%s/%s/%s/%s/result.rsparsepro.log %s",
        s, Dir, block, suffix[k], "|grep vare |tail -n 1 |cut -d ' ' -f 6"), intern = T)

      output = read.delim(sprintf(
        "results/setting_%s/%s/%s/%s/result.rsparsepro.txt",
        s, Dir, block, suffix[k]))

      data$PIP = output$PIP
      data$CS = NA

      result_CS95 = data.frame("block" = block,
                               "coverage" = 0,
                               "size" = 0,
                               "min.abs.corr" = 0)[-1,]
      if(sum(output$cs > 0) > 0)
      {
        for(l in 1:max(output$cs))
        {
          set = which(output$cs == l)
          if(length(set) == 1)
          {
            purity = 1
          } else {
            R_sub = R[[idx[k]]][set,set]
            purity = min(abs(R_sub[upper.tri(R_sub)]))
          }

          if(purity >= 0.5)
          {
            data$CS[set] = l
            result_CS95 %<>% rbind(data.frame(
              "block" = block,
              "coverage" = sum(data[set,8+idx[k]]),
              "size" = length(set),
              "min.abs.corr" = purity))
          }
        }
      }
      
      write_delim(result_CS95, sprintf(
        "results/setting_%s/%s/%s/result_CS95_%s.txt",
        s, Dir, block, suffix[k]), delim = '\t')

      saveRDS(data, file = sprintf(
        "results/setting_%s/%s/%s/data_%s.RData",
        s, Dir, block, suffix[k]))

      result = data.frame(
        "block" = block,
        "time" = difftime(toc, tic, units = "secs"),
        "vare" = vare)
      
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
        if(sum(!is.na(data[[1]]$CS)) > 0 & sum(!is.na(data[[2]]$CS)) > 0)
        {
          for(l1 in unique(na.omit(data[[1]]$CS)))
          {
            set1 = which(data[[1]]$CS == l1)
            for(l2 in unique(na.omit(data[[2]]$CS)))
            {
              set2 = which(data[[2]]$CS == l2)
              if(length(intersect(set1,set2)) > 0)
              {
                purity = min(min(abs(R[[1]][set1,set1])), min(abs(R[[2]][set2,set2])))
                
                result_CS95_shared %<>% rbind(data.frame(
                  "block" = block,
                  "coverage" = sum(data[[1]][intersect(set1,set2),9] &
                                     data[[1]][intersect(set1,set2),10]),
                  "size" = length(union(set1,set2)),
                  "min.abs.corr" = purity)
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
      if(file.exists(sprintf(
        "results/setting_%s/%s/%s/result_CS95_%s.txt",
        s, Dir, block, suffix[k])))
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
  
  for(m in 1:2)
  {
    data_2 = foreach(block = select_block, .combine = "rbind") %dopar%
    {
      readRDS(sprintf(
        "results/setting_%s/%s/%s/data_%s.RData",
        s, Dir, block, suffix[m+1]))
    }
    
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
      "results/setting_%s/%s/calibration_N2-%d.RData",
      s, Dir, N2_seq[m]))
  }
}
