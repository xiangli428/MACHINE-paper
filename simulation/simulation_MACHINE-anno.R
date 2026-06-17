options(stringsAsFactors = F, check.names = F)

library(readr)
library(magrittr)
library(dplyr)
library(Matrix)
library(foreach)
library(doParallel)
library(MACHINE)

setwd("simulation")

select_block_info = read.delim("data/select_block_info.txt")
select_block_info %<>% arrange(num_variants)
select_block_info[134:200,] %<>% arrange(-num_variants)
select_block = sort(select_block_info$block)

pids = c("EUR", "AFR", "EAS")

n_causal = data.frame("n" = c(5,3,1),
                      "n1" = c(0,1,2),
                      "n2" = c(0,1,2))

N1 = 2e5
N2_seq = c(2e4,2e5)

group_bound = data.frame("l" = c(0,0.1,0.5,0.9),
                         "u" = c(0.1,0.5,0.9,1.1))

registerDoParallel(67)

for(anno in c("gLDSC","polyfun"))
{
  Dir = sprintf("MACHINE-%s", anno)
  
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
      a1 = gLDSC_var_sub[,1] / mean(gLDSC_var_sub[,1]) * 0.005
      
      for(p in 2:3)
      {
        for(m in 1:2)
        {
          suffix = sprintf("N%d-%d", p, N2_seq[m])
          dir.create(sprintf("results/setting_%s/%s/%s/%s",
                             s, Dir, suffix, block), recursive = T)
          
          data = data_full[,c(1:6,14:21,17+2*p+m)]
          Z = data[,14:15]
          n = c(N1,N2_seq[m])
          for(k in 1:2)
          {
            Z[,k] = sqrt(n[k] - 1) * Z[,k] / sqrt(n[k] - Z[,k]^2)
          }
          
          a2 = gLDSC_var_sub[,2*p+m-3] / mean(gLDSC_var_sub[,2*p+m-3]) * 0.005
          a = (a1+a2) / 2
          c = cbind(a1,a2)
          c = c / rowSums(c) * 0.4
          
          MACHINE = CreateMACHINEObject(Z, 
                                        R[c(1,p)], 
                                        n,
                                        SNP_ID = data$rsid,
                                        POP_ID = pids[c(1,p)],
                                        in_sample_LD = T,
                                        a = a,
                                        c = c)
          
          tic = Sys.time()
          MACHINE = MACHINE_MCMC(MACHINE, mcmc_n = 5500, burn_in = 500, thin = 2,
                                 stepsize = 1, get_CS = F)
          while(max(MACHINE@mcmc_samples$PSRF_beta) > 1.2)
          {
            MACHINE = MACHINE_MCMC(MACHINE, mcmc_n = 1000, burn_in = 1000, thin = 2,
                                   stepsize = 1, get_CS = F)
          }
          MACHINE@CS = MACHINE_CS(MACHINE, MACHINE@coverage, MACHINE@purity)
          toc = Sys.time()
          
          data$CL_1 = MACHINE@CL[,1]
          data$CL_2 = MACHINE@CL[,2]
          
          data$CS_2 = data$CS_1 = NA
          
          for(k in 1:2)
          {
            if(length(MACHINE@CS[[k]]$sets) > 0)
            {
              for(l in 1:length(MACHINE@CS[[k]]$sets))
              {
                data[MACHINE@CS[[k]]$sets[[l]],sprintf("CS_%d", k)] = l
              }
            }
          }
          
          saveRDS(data, file = sprintf(
            "results/setting_%s/%s/%s/%s/data.RData",
            s, Dir, suffix, block))
          
          result_CS95 = foreach(k = 1:2) %do%
          {
            data.frame("block" = block,
                       "coverage" = 0,
                       "size" = 0,
                       "min.abs.corr" = 0)[-1,]
          }
          
          for(k in 1:2)
          {
            if(length(MACHINE@CS[[k]]$sets) > 0)
            {
              result_CS95[[k]] = foreach(l = 1:length(MACHINE@CS[[k]]$sets), 
                                         .combine = "rbind") %do%
                {
                  data.frame("block" = block,
                             "coverage" = sum(data[MACHINE@CS[[k]]$sets[[l]],9+k]),
                             "size" = length(MACHINE@CS[[k]]$sets[[l]]),
                             "min.abs.corr" = MACHINE@CS[[k]]$purity$min.abs.corr[l])
                }
            }
            
            write_delim(result_CS95[[k]], sprintf(
              "results/setting_%s/%s/%s/%s/result_CS95_%d.txt",
              s, Dir, suffix, block, k), delim = '\t')
          }
          
          result_CS95_cross = data.frame("block" = block,
                                         "coverage" = 0,
                                         "size" = 0,
                                         "min.abs.corr" = 0)[-1,]
          result_CS95_cross %<>% rbind(
            foreach(k = 1:2, .combine = "rbind") %do%
            {
              if(length(MACHINE@CS[[k]]$sets) > 0)
              {
                foreach(l = 1:length(MACHINE@CS[[k]]$sets), 
                        .combine = "rbind") %do%
                {
                  data.frame("block" = block,
                             "coverage" = sum(data$causal_1[MACHINE@CS[[k]]$sets[[l]]] |
                                                data$causal_2[MACHINE@CS[[k]]$sets[[l]]]),
                             "size" = length(MACHINE@CS[[k]]$sets[[l]]),
                             "min.abs.corr" = MACHINE@CS[[k]]$purity$min.abs.corr[l])
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
          if(length(MACHINE@CS[[1]]$sets) > 0 & length(MACHINE@CS[[2]]$sets) > 0)
          {
            for(l1 in 1:length(MACHINE@CS[[1]]$sets))
            {
              set1 = MACHINE@CS[[1]]$sets[[l1]]
              for(l2 in 1:length(MACHINE@CS[[2]]$sets))
              {
                set2 = MACHINE@CS[[2]]$sets[[l2]]
                if(length(intersect(set1,set2)) > 0)
                {
                  result_CS95_shared %<>% rbind(data.frame(
                    "block" = block,
                    "coverage" = sum(data$causal_1[intersect(set1,set2)] &
                                       data$causal_2[intersect(set1,set2)]),
                    "size" = length(union(set1,set2)),
                    "min.abs.corr" = min(MACHINE@CS[[1]]$purity$min.abs.corr[l1],
                                         MACHINE@CS[[2]]$purity$min.abs.corr[l2]))
                  )
                }
              }
            }
          }
          
          write_delim(result_CS95_shared, sprintf(
            "results/setting_%s/%s/%s/%s/result_CS95_shared.txt",
            s, Dir, suffix, block), delim = '\t')
          
          result = data.frame(
            "block" = block,
            "max_PSRF" = max(MACHINE@mcmc_samples$PSRF_beta),
            "mcmc_n" = MACHINE@mcmc_samples$n_samples + MACHINE@mcmc_samples$n_burnin,
            "time" = difftime(toc, tic, units = "secs"),
            "CS95_causal" = sum((!is.na(data$CS_1) | !is.na(data$CS_2)) &
                                  (data$causal_1 | data$causal_2)),
            "CS95_causal_1" = sum(!is.na(data$CS_1) & data$causal_1),
            "CS95_causal_2" = sum(!is.na(data$CS_2) & data$causal_2),
            "CS95_causal_0" = sum(!is.na(data$CS_1) & !is.na(data$CS_2) &
                                    (data$causal_1 & data$causal_2)))
          
          write_delim(result, sprintf(
            "results/setting_%s/%s/%s/%s/result.txt",
            s, Dir, suffix, block), delim = '\t')
          
          MACHINE@mcmc_samples$beta = NULL
          MACHINE@mcmc_samples$tau = NULL
          MACHINE@mcmc_samples$psi = NULL
          saveRDS(MACHINE, file = sprintf(
            "results/setting_%s/%s/%s/%s/MACHINE.RData",
            s, Dir, suffix, block))
          
          NULL
        }
      }
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
  
      CL = pmax(data$CL_1, data$CL_2)
      calibration[["cross"]] = foreach(g = 1:4, .combine = "rbind") %dopar%
      {
        idx = which(CL >= group_bound$l[g] & CL < group_bound$u[g])
        n = length(idx)
  
        data.frame("group" = g,
                   "n" = n,
                   "Expected" = sum(CL[idx]) / n,
                   "Prop" = sum((data$causal_1 | data$causal_2)[idx]) / n)
      }
  
      CL = pmin(data$CL_1, data$CL_2)
      calibration[["shared"]] = foreach(g = 1:4, .combine = "rbind") %dopar%
      {
        idx = which(CL >= group_bound$l[g] & CL < group_bound$u[g])
        n = length(idx)
  
        data.frame("group" = g,
                   "n" = n,
                   "Expected" = sum(CL[idx]) / n,
                   "Prop" = sum((data$causal_1 & data$causal_2)[idx]) / n)
      }
  
      CL = data$CL_1
      calibration[["pop1"]] = foreach(g = 1:4, .combine = "rbind") %dopar%
      {
        idx = which(CL >= group_bound$l[g] & CL < group_bound$u[g])
        n = length(idx)
  
        data.frame("group" = g,
                   "n" = n,
                   "Expected" = sum(CL[idx]) / n,
                   "Prop" = sum(data$causal_1[idx]) / n)
      }
  
      CL = data$CL_2
      calibration[["pop2"]] = foreach(g = 1:4, .combine = "rbind") %dopar%
      {
        idx = which(CL >= group_bound$l[g] & CL < group_bound$u[g])
        n = length(idx)
  
        data.frame("group" = g,
                   "n" = n,
                   "Expected" = sum(CL[idx]) / n,
                   "Prop" = sum(data$causal_2[idx]) / n)
      }
  
      saveRDS(calibration, sprintf(
        "results/setting_%s/%s/%s/calibration.RData",
        s, Dir, suffix))
    }
  }
}
