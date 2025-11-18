source('gLDSC/make_left_right_seperate3.R')

setwd("real_data/scz2022")

pids = c("EUR", "EAS")
N = c("EUR" = 130644, "EAS" = 30761)

for(pid in pids)
{
  dir.create(sprintf("gLDSC_results_new/%s", pid), recursive = T)
  for(i in 1:22)
  {
    panel = sprintf("~/Documents/GWAS/project_4/annotations/LDSM/%s/chr%s",
                    pid, i)
    
    data = read.delim(gzfile(sprintf("data/%s/chr%s.txt.gz", pid, i)))
    gwas = data.frame("SNP" = data$rsid,
                      "A1" = data$first_allele,
                      "A2" = data$alternative_alleles,
                      "N" = N[pid],
                      "Z" = data$Z)
    
    out_path = sprintf("gLDSC_results/%s/chr%s", pid, i)
    
    #for one SNP
    result=gldsc2(panel = panel, gwas = gwas, jackknife = F)
    saveRDS(result, paste0(out_path, '.Rdata'))
  }
}