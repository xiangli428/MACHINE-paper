options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)
library(XMAP)

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

group_bound = data.frame("l" = c(0,0.1,0.5,0.9),
                         "u" = c(0.1,0.5,0.9,1.1))

registerDoParallel(50)

Dir = "XMAP"

for(s in 1:3)
{
  foreach(block = select_block_info$block, .combine = "c", .inorder = F) %dopar%
  {
    data_full = readRDS(sprintf("data/setting_%s/%d.RData", s, block))
    R = readRDS(sprintf("data/LD/%d/LD_UKBB.RData", block))
    
    for(k in 1:3)
    {
      R[[k]] = as.matrix(R[[k]])
      diag(R[[k]]) = 1
    }
    
    for(p in 2:3)
    {
      for(m in 1:2)
      {
        suffix = sprintf("N%d-%d", p, N2_seq[m])
        dir.create(sprintf("results/setting_%s/%s/%s/%s",
                           s, Dir, suffix, block), recursive = T)
        
        data = data_full[,c(1:6,14:21,17+2*p+m)]
        Z = data[,14:15]
        N = c(N1,N2_seq[m])
        
        tic = Sys.time()
        fit_XMAP <- XMAP(R = simplify2array(R[c(1,p)]), z = Z,
                         n = c(N1,N2_seq[m]), K = 5,
                         Omega = diag(1e-100, 2, 2),
                         Sig_E = c(1,1))
        CS = foreach(k = c(1,p)) %do%
        {
          get_CS(fit_XMAP, Xcorr = R[[k]], coverage = 0.95, min_abs_corr = 0.5)
        }
        toc = Sys.time()
        
        saveRDS(fit_XMAP, file = sprintf(
          "results/setting_%s/%s/%s/%s/fit_XMAP.RData",
          s, Dir, suffix, block))
        saveRDS(CS, file = sprintf(
          "results/setting_%s/%s/%s/%s/CS.RData",
          s, Dir, suffix, block))
        
        data$PIP = get_pip(fit_XMAP$gamma)
        
        result_CS95 = foreach(k = 1:2) %do%
        {
          data.frame("block" = block,
                     "coverage" = 0,
                     "size" = 0,
                     "min.abs.corr" = 0)[-1,]
        }
        
        for(k in 1:2)
        {
          data[sprintf("CS_%s", k)] = NA
          if(length(CS[[k]]$cs) > 0)
          {
            for(l in 1:length(CS[[k]]$cs))
            {
              data[CS[[k]]$cs[[l]], sprintf("CS_%s", k)] =
                CS[[k]]$cs_index[l]
              result_CS95[[k]] %<>% rbind(data.frame(
                "block" = block,
                "coverage" = sum(data[CS[[k]]$cs[[l]],sprintf("causal_%s", k)]),
                "size" = length(CS[[k]]$cs[[l]]),
                "min.abs.corr" = CS[[k]]$purity$min.abs.corr[l]))
            }
          }
          
          write_delim(result_CS95[[k]], sprintf(
            "results/setting_%s/%s/%s/%s/result_CS95_%d.txt",
            s, Dir, suffix, block, k), delim = '\t')
        }
        
        saveRDS(data, file = sprintf(
          "results/setting_%s/%s/%s/%s/data.RData",
          s, Dir, suffix, block))
        
        result_CS95_cross = data.frame("block" = block,
                                       "coverage" = 0,
                                       "size" = 0,
                                       "min.abs.corr" = 0)[-1,]
        result_CS95_cross %<>% rbind(
          foreach(k = 1:2, .combine = "rbind") %do%
          {
            if(length(CS[[k]]$cs) > 0)
            {
              foreach(l = 1:length(CS[[k]]$cs), .combine = "rbind") %do%
              {
                data.frame("block" = block,
                           "coverage" = sum(data[CS[[k]]$cs[[l]],10] |
                                              data[CS[[k]]$cs[[l]],11]),
                           "size" = length(CS[[k]]$cs[[l]]),
                           "min.abs.corr" = CS[[k]]$purity$min.abs.corr[l])
              }
            }
          })
        write_delim(result_CS95_cross, sprintf(
          "results/setting_%s/%s/%s/%s/result_CS95_cross.txt",
          s, Dir, suffix, block), delim = '\t')
        
        result_CS95_shared = data.frame("block" = block,
                                        "coverage" = 0,
                                        "size" = 0,
                                        "min.abs.corr" = 0)[-1,]
        if(sum(!is.na(data$CS_1) & !is.na(data$CS_2)) > 0)
        {
          for(l in unique(data$CS_1[!is.na(data$CS_1) & !is.na(data$CS_2)]))
          {
            set = which(data$CS_1 == l)
            
            purity = min(CS[[1]]$purity$min.abs.corr[which(CS[[1]]$cs_index==l)],
                         CS[[2]]$purity$min.abs.corr[which(CS[[2]]$cs_index==l)])
            
            result_CS95_shared %<>% rbind(data.frame(
              "block" = block,
              "coverage" = sum(data[set,10] & data[set,11]),
              "size" = length(set),
              "min.abs.corr" = purity))
          }
        }
        
        write_delim(result_CS95_shared, sprintf(
          "results/setting_%s/%s/%s/%s/result_CS95_shared.txt",
          s, Dir, suffix, block), delim = '\t')
        
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
          "results/setting_%s/%s/%s/%s/result.txt",
          s, Dir, suffix, block), delim = '\t')
      }
    }
    
    NULL
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