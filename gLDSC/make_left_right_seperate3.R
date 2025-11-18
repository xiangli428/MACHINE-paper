#' Title
#'
#' @param panel Path of LDSM.Rdata file.
#' @param gwas Path of GWAS summary input file.
#' @param jackknife Whether use Jackknife to estimate standard error.
#' @param intercept Whether consider confounding bias in the analysis.
#' @param numCores Number of cores used in analysis.
#'
#' @return
#' @export
#'
#' @examples
source('~/gLDSC_scripts/gldsc/R/middle_fun_big_sep.R')
#panel='~/Documents/GWAS/data/UKBB/5/imputed_genotype_unique_info-0.8_maf-0.001_hwe-1e-6/LD/chr22/LDSM_chr22_EAS.Rdata'
#gwas_pre_loc='/home/r9user9/Documents/GWAS/project_4/real_data/scz2022/data/asian/chr22.txt.gz'
#out='/home/r9user9/Documents/GWAS/project_4/gLDSC_result/chr22_pioneer/EAS/enrichment_res.Rdata'

#pre-process GWAS input
#file1<-gzfile(gwas_pre_loc,'rt')
#gwas_pre=read.table(file1,header=T)
#gwas=gwas_pre[,c('rsid','first_allele','alternative_alleles')]
#names(gwas)=c('SNP','A1','A2')
#gwas[,'N']=gwas_pre[,'NCON']+gwas_pre[,'NCAS']
#gwas[,'Z']=gwas_pre[,'BETA']/gwas_pre[,'SE']
#close(file1)

gldsc2<-function(outpath='',panel,gwas,jackknife=T,intercept=T,numCores=50){
  start_time=Sys.time()
  message('Analysis begin')
  #panel1=paste0(panel,'.Rdata')
  ref.pannel<-readRDS(paste0(panel,'.Rdata'))
  #ref.pannel<-readRDS('/home/r14user2/Documents/GWAS/project_4/annotations/LDSM/asian/chr1.Rdata')
  message('load reference panel...')
  snplist.pan<-unlist(ref.pannel[['snp.list']])
  message(paste0('read ',length(snplist.pan),' SNPs from reference panel'))
  registerDoParallel(numCores)
  x.sum<-Reduce(sum,lapply(ref.pannel[['x.gls']], sum))
  
  
  #gls
  #gwas.df<-read.table(gwas,header = T,fill = TRUE)
  gwas.df<-gwas
  message('read GWAS...')
  gwas.df<-na.omit(gwas.df)
  message(paste0('read ',nrow(gwas.df),' SNPs from GWAS'))
  snp.used=intersect(snplist.pan,gwas.df$SNP)
  message(paste0('After merging with reference panel ',length(snp.used),' SNPs remain'))
  #
  merged.loc<-match(snp.used,snplist.pan)
  
  gwas.df<-gwas.df[match(snp.used, gwas.df$SNP),]
  y.sum<-sum(gwas.df[,'Z']**2-1,na.rm = T)
  raw.Ntau=y.sum/x.sum
  
  gwasN<-median(gwas.df$N,na.rm = T)
  gls.left<-gls.right<-0
  raw.result.all<-list()
  for (chr in 'used chr') {
    #chr.exact<-substr(ref.pannel[['chrs']][chr],nchar('ldblk_1kg_')+1,nchar(ref.pannel[['chrs']][chr])-nchar('.hdf5'))
    chr.exact=chr
    #chr.exact<-ref.pannel[['chrs']][chr]
    #block.names.all<-names(ref.pannel[['LDSM']])[grep(x = names(ref.pannel[['LDSM']]), pattern = paste0(chr.exact,'_'))]
    panel2=list.files(panel)
    
    total.block<-length(panel2)
    
    if(total.block>0){
	    #output_loc=''
	    raw.result<-foreach (i=1:total.block) %dopar% gls.left.right(i=i,y=gwas.df,LDpath=panel,
                                                                   raw.Ntau = raw.Ntau,intercept = intercept
      )
      for (batch in 1:length(raw.result)) {
        gls.left<-gls.left+raw.result[[batch]][[1]]
        gls.right<-gls.right+raw.result[[batch]][[2]]
      }
      raw.result.all[[chr.exact]]<-raw.result
    }
    message(paste0(chr,' calculation complited'))
  }
  
  
  #first time result
  result.one<-list()
  #result.one<-gls.estimator(left = gls.left,right = gls.right,
  #                          A=ref.pannel[['Atotal']],M.anno = ref.pannel[['M.anno']],anno.total = ref.pannel[['anno.total']],N = gwasN,intercept = intercept)
  #result.one[['left_part']]=gls.left
  #result.one[['right_part']]=gls.right
  result.one[['Atotal']]=ref.pannel[['Atotal']]
  message('point estimate completed')
  #result
  if(jackknife==F){
    message(Sys.time()-start_time)
    return(raw.result)
    #return(result.one)
  }else if(jackknife==T){
    message('start Jackknife estimation')
    #jackknife starts here
    jack.tau<-jack.h2<-jack.e<-jack.intersecp<-jack.stat<-NA
    jack.time=0
    #chr in 1:length(ref.pannel[['chrs']])
    for (chr in 15) {
      #chr.exact<-substr(ref.pannel[['chrs']][chr],nchar('ldblk_1kg_')+1,nchar(ref.pannel[['chrs']][chr])-nchar('.hdf5'))
      chr.exact=ref.pannel[['chrs']][chr]
      raw.result.temp=raw.result.all[[chr.exact]]
      total.block2<-length(raw.result.temp)
      jack.result<-foreach (i=1:total.block2) %dopar% oneout(i,raw.result=raw.result.temp,left.total = gls.left,right.total = gls.right,
                                                             A=ref.pannel[['Atotal']],M.anno = ref.pannel[['M.anno']],
                                                             anno.total = ref.pannel[['anno.total']],N = gwasN,intercept = intercept)
      for (i in 1:total.block2) {
        jack.tau<-rbind(jack.tau,as.vector(jack.result[[i]][[1]]))
        jack.h2<-rbind(jack.h2,as.vector(jack.result[[i]][[2]]))
        jack.e<-rbind(jack.e,as.vector(jack.result[[i]][[3]]))
        #jack.intersecp<-c(jack.intersecp,jack.result[[i]][[4]])
        #jack.stat<-rbind(jack.stat,as.vector(jack.result[[i]][[5]]))
        jack.time=jack.time+1
      }
    }
    jack.tau.var<-apply(jack.tau[-1,], 2, var)*(jack.time-1)^2/jack.time
    jack.h2.var<-apply(jack.h2[-1,], 2, var)*(jack.time-1)^2/jack.time
    jack.e.var<-apply(jack.e[-1,], 2, var)*(jack.time-1)^2/jack.time
    #jack.intersecp.var<-var(jack.intersecp[-1])*(jack.time-1)^2/jack.time
    #jack.stat.var<-apply(jack.stat[-1,], 2, var)*(jack.time-1)^2/jack.time
    
    
    result.one[['tau.var.jack']]=jack.tau.var
    result.one[['h2.var.jack']]=jack.h2.var
    result.one[['e.var.jack']]=jack.e.var
    #result.one[['inter.var.jack']]=jack.intersecp.var
    #result.one[['P']]=pnorm(result.one[['e.stat']]/sqrt(jack.stat.var),lower.tail=F)
    message(Sys.time()-start_time)
    #return(result.one)
    return(raw.result)
  }else{print('Jackknife should be TRUE of FALSE')}
}


