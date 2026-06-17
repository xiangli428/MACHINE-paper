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

registerDoParallel(20)

Dir = "SuSiEx"
wkd = "simulation/results"

lpexp = function(x,y)
{
  if(x > y)
  {
    x + log1p(exp(y-x))
  } else {
    y + log1p(exp(x-y))
  }
}

for(s in 1:3)
{
  for(block in select_block)
  {
    dir.create(sprintf("results/setting_%s/%s/%s", s, Dir, block), recursive = T)
    
    data_full = readRDS(sprintf("data/setting_%s/%d.RData", s, block))
    R = readRDS(sprintf("data/LD/%d/LD_UKBB.RData", block))
    
    for(k in 1:3)
    {
      R[[k]] = as.matrix(R[[k]])
      diag(R[[k]]) = 1
      writeBin(as.vector(R[[k]]), file(sprintf(
        "results/setting_%s/%s/%s/%s.ld.bin", s, Dir, block, pids[k]), "wb"),
        size = 4)
    }
    
    for(k in 1:5)
    {
      data = data_full[,c(1,5,2:4)]
      n = c(N1, N2_seq, N2_seq)[k]
      data$betaHat = data_full[,20+k] / sqrt(n)
      data %<>% mutate(sigmaHat = sqrt(1 - betaHat^2) / sqrt(n-1),
                       pval = 2 * pnorm(-abs(betaHat / sigmaHat)))
      write_delim(data, sprintf("results/setting_%s/%s/%s/data_%s.txt",
                                s, Dir, block, suffixes[k]), delim = '\t')
    }
    
    for(p in 2:3)
    {
      for(m in 1:2)
      {
        suffix = sprintf("N%d-%d", p, N2_seq[m])
        dir.create(sprintf("results/setting_%s/%s/%s/%s",
                           s, Dir, suffix, block), recursive = T)
        
        tic = Sys.time()
        system(paste(
          "~/Documents/GWAS/packages/SuSiEx/bin/SuSiEx ",
          "--sst_file=",
          sprintf("%s/setting_%s/%s/%s/data_%s.txt,", wkd,s,Dir,block,suffixes[1]),
          sprintf("%s/setting_%s/%s/%s/data_%s.txt ", wkd,s,Dir,block,suffixes[2*p+m-3]),
          sprintf("--n_gwas=%d,%d ",N1,N2_seq[m]),
          "--ref_file=",
          sprintf("%s/%s/genotype_bed/chr1,",
                  "/home/r14user2/Documents/GWAS/data/UKBB/1",
                  "imputed_genotype_unique_info-0.8_maf-0.001_hwe-1e-6"),
          sprintf("%s/%s/%s/genotype_bed/chr1 ",
                  "/home/r14user2/Documents/GWAS/data/UKBB", p+2,
                  "imputed_genotype_unique_info-0.8_maf-0.001_hwe-1e-6"),
          "--ld_file=",
          sprintf("%s/setting_%s/%s/%s/%s,", wkd,s,Dir,block,pids[1]),
          sprintf("%s/setting_%s/%s/%s/%s ", wkd,s,Dir,block,pids[p]),
          "--out_dir=",
          sprintf("%s/setting_%s/%s/%s/%s ", wkd,s,Dir,suffix,block),
          "--out_name=output --chr=1 ",
          sprintf("--bp=%d,%d ", min(data_full$position)-1,max(data_full$position)+1),
          "--chr_col=1,1 --snp_col=2,2 --bp_col=3,3 --a1_col=4,4 --a2_col=5,5 ",
          "--eff_col=6,6 --se_col=7,7 --pval_col=8,8 ",
          "--plink=",
          "/home/r14user2/Documents/GWAS/packages/SuSiEx/utilities/plink ",
          "--keep-ambig=True --maf=0.0001 --threads=1  --pval_thresh=1", sep = ''))
        toc = Sys.time()
        
        data = data_full[,c(1:6,14:21,17+2*p+m)]
        
        output_SNP = read.delim(sprintf(
          "results/setting_%s/%s/%s/%s/output.snp",
          s, Dir, suffix, block), check.names = F)
        output_CS = read.delim(sprintf(
          "results/setting_%s/%s/%s/%s/output.cs",
          s, Dir, suffix, block), check.names = F)
        output_summary = read.delim(sprintf(
          "results/setting_%s/%s/%s/%s/output.summary",
          s, Dir, suffix, block), check.names = F, skip = 1)
        
        L = nrow(output_summary)
        
        rownames(data) = data$rsid
        data$PIP_0 = data$PIP_2 = data$PIP_1 = data$PIP = 0
        
        data$CS_2 = data$CS_1 = data$CS = NA
        
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
        
        if(L > 0)
        {
          probs = foreach(l = 1:L, .combine = "rbind") %do%
          {
            lbf_1 = foreach(x = sort(output_SNP[,3*l+1]), .combine = "lpexp") %do%
            {
              x
            }
            lbf_2 = foreach(x = sort(output_SNP[,3*l+2]), .combine = "lpexp") %do%
            {
              x
            }
            lbf = foreach(x = sort(output_SNP[,3*l+1]+output_SNP[,3*l+2]), 
                            .combine = "lpexp") %do%
            {
              x
            }
            
            lsum = foreach(x = c(lbf_1, lbf_2, lbf), .combine = "lpexp") %do%
            {
              x
            }
            
            c(exp(lbf_1 - lsum), exp(lbf_2 - lsum), exp(lbf - lsum))
          }
          
          if(L == 1) probs = matrix(probs, nrow = L, ncol = 3)
          
          data$PIP = 1 - foreach(l = 1:L, .combine = "*") %do%
          {
            1 - output_SNP[,3*l]
          }
          
          data$PIP_1 = 1 - foreach(l = 1:L, .combine = "*") %do%
          {
            1 - output_SNP[,3*l] * (probs[l,1] + probs[l,3])
          }
          
          data$PIP_2 = 1 - foreach(l = 1:L, .combine = "*") %do%
          {
            1 - output_SNP[,3*l] * (probs[l,2] + probs[l,3])
          }
          
          data$PIP_0 = 1 - foreach(l = 1:L, .combine = "*") %do%
          {
            1 - output_SNP[,3*l] * probs[l,3]
          }
          
          for(l in 1:L)
          {
            set = which(is.element(data$rsid, output_CS$SNP[output_CS$CS_ID == l]))
            data$CS[set] = l
            
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
              if(probs[l,k] + probs[l,3] >= 0.8)
              {
                data[set, sprintf("CS_%s", k)] = l
                result_CS95[[k]] %<>% rbind(data.frame(
                  "block" = block,
                  "coverage" = sum(data[set,sprintf("causal_%s", k)]),
                  "size" = length(set),
                  "min.abs.corr" = purities[k]))
                
                result_CS95_cross %<>% rbind(data.frame(
                  "block" = block,
                  "coverage" = sum(data[set,10] | data[set,11]),
                  "size" = length(set),
                  "min.abs.corr" = purities[k]))
              }
            }
            
            if(probs[l,1] + probs[l,3] >= 0.8 & probs[l,2] + probs[l,3] >= 0.8)
            {
              result_CS95_shared %<>% rbind(data.frame(
                "block" = block,
                "coverage" = sum(data[set,10] & data[set,11]),
                "size" = length(set),
                "min.abs.corr" = max(purities)))
            } else if(probs[l,1] + probs[l,3] < 0.8 & probs[l,2] + probs[l,3] < 0.8) {
              result_CS95_cross %<>% rbind(data.frame(
                "block" = block,
                "coverage" = sum(data[set,10] | data[set,11]),
                "size" = length(set),
                "min.abs.corr" = max(purities)))
            }
          }
        }
        
        saveRDS(data, file = sprintf(
          "results/setting_%s/%s/%s/%s/data.RData",
          s, Dir, suffix, block))
        
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
        
        result = data.frame(
          "block" = block,
          "time" = difftime(toc, tic, units = "secs"),
          "CS95_causal" = sum(!is.na(data$CS) & (data$causal_1 | data$causal_2)),
          "CS95_causal_1" = sum(!is.na(data$CS_1) & data$causal_1),
          "CS95_causal_2" = sum(!is.na(data$CS_2) & data$causal_2),
          "CS95_causal_0" = sum((!is.na(data$CS_1) & !is.na(data$CS_2)) &
                                (data$causal_1 & data$causal_2)))
        
        write_delim(result, sprintf(
          "results/setting_%s/%s/%s/%s/result.txt",
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
      "results/setting_%s/%s/%s/calibration.RData",
      s, Dir, suffix))
  }
}