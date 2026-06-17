library(magrittr)

source('~/gLDSC_scripts/gldsc/R/middle_fun_big.R')

mkLDSM <- function(LD.path, anno.path, maf.path=NULL, annotation='all',
                   snp=NULL, MAF=NULL, numCores=30, chr_used=1, out_path){
  #register cores
  registerDoParallel(numCores)
  
  #Check reference files#####
  LD.files <- list.files(LD.path)
  if(any(grepl(x = LD.files, pattern = "chr*"))){
    chrs <- LD.files[grep(x = LD.files, pattern = "chr*")]
  }else{
    error.message <- "Ref pannel (LDfiles) input error"
    stop(error.message)
  }
  
  anno.files<- list.files(anno.path)
  if(any(grepl(x = anno.files, pattern = "*.annot"))){
    annos <- anno.files[grep(x = anno.files, pattern = "*.annot")]
  }else{
    error.message <- "Ref pannel (ANNOfile) input error"
    stop(error.message)
  }
  
  #Preparation#####
  if(is.null(snp)==F){
    gwas.df<-read.table(snp,header = T)
    gwas.df<-na.omit(gwas.df)
    gwas.snp<-as.vector(gwas.df[,'SNP'])
    print(paste('load',length(gwas.snp),'SNPs from custom snplist',sep = ' '))
  }
  
  #Starts#####
  snp.list<-LDmatrix<-LDSM<-x.gls<-list()
  Atotal<-NA
  start_time <- Sys.time()
  #chr in 1:chr.length
  for (chr in chr_used) {
    #read block info
    chr.exact<-LD.files[chr]
    #block.info<-read.table(paste0(LD.path,'/',LD.files[chr],'/reblock.txt'),header = T)
    block.info<-list.files(paste0(LD.path,'/',LD.files[chr]))
    #chr.exact<-substr(chrs[chr],nchar('ldblk_1kg_')+1,nchar(chrs[chr])-nchar('.hdf5'))
    #read ref by chr
    LDinfo<-paste0(LD.path,'/',LD.files[chr])
    #######
    ##remove .txt files, only usable for this simulation
    block.info = 1:200
    #######
    #LDinfo<-h5ls(file = paste(LD.path,chrs[chr], sep = '/'),recursive=F)
    
    #file1<-gzfile(paste(anno.path,annos[chr], sep = '/'),'rt')
    file1<-paste(anno.path,annos[chr], sep = '/')
    #file3<-gzfile(paste(maf.path,mafs[chr], sep = '/'),'rt')
    
    anno2<-read.table(file1, header = T,sep = '')
    #maf<-read.table(file3,header = T,sep = '')
    #file4<-gzfile(paste0(LD.path,'/',LD.files[chr],'/variant_reblock.txt.gz'),'rt')
    #maf<-read.table(file4,header = T,sep = '')
    
    #close(file1)
    #close(file2)
    #close(file3)
    #close(file4)
    #select annotation
    if(annotation[1]=='all'){used.a<-names(anno2)[-2:-1]
    }else{used.a<-annotation}
    
    #select chosen annotation
    #anno<-anno2[,c(names(anno2)[1:4],used.a)]
    #maf.snp<-as.vector(maf[maf$FRQ<(1-MAF)&maf$FRQ>MAF,'SNP'])
    #ldsc.snp<-ldsc$SNP
    #grab rsid from file 4
    #maf.snp<-as.vector(maf[,'rsid'])
    
    #generate by block
    # out_path2=paste0(out_path,'_chr',chr.exact)
    if(is.null(snp)){gwas.snp<-maf.snp<-1}
    raw.result=foreach (batch=block.info) %dopar% 
      refpannel.par(batch, LDinfo = LDinfo, LD.path = LD.path,
                    chrs = chrs[chr], out = out_path,
                    anno = anno2,gwas.snp = gwas.snp, used.a = used.a)
    #record result
    for (batch in 1:length(raw.result)) {
      if(raw.result[[batch]][[1]]==1){
        #batch.name<-LDinfo[batch,'name']
        batch.name<-batch
        snp.list[[paste(chr.exact,batch.name,sep = '_')]]<-raw.result[[batch]][[2]]
        #    LDmatrix[[paste(chr.exact,batch.name,sep = '_')]]<-raw.result[[batch]][[3]]
        #    LDSM[[paste(chr.exact,batch.name,sep = '_')]]<-raw.result[[batch]][[4]]
        #    x.gls[[paste(chr.exact,batch.name,sep = '_')]]<-raw.result[[batch]][[5]]
        #    
        x.gls[[paste(chr.exact,batch.name,sep = '_')]]<-raw.result[[batch]][[3]]
        Anno.temp<-anno2[match(raw.result[[batch]][[2]],anno2[,'SNP']),used.a]
        Anno.temp<-as.matrix(Anno.temp)
        Atotal<-rbind(Atotal,Anno.temp)
      }
    }
    print(paste0('LD Score Matrix of CHR ',chr,' has been complited'))
  }
  
  Atotal<-Atotal[-1,]
  snplist<-unlist(snp.list)
  M<-length(snplist)
  anno.total<-apply(Atotal, 2, sum)
  M.anno<-anno.total/M
  print(paste(M,'SNPs left for analysis',sep = ' '))
  print(Sys.time()-start_time)
  
  #output
  result<-list()
  #result[['LDinfo']]<-LDinfo
  #result[['LDmatrix']]<-LDmatrix
  #result[['LDSM']]<-LDSM
  result[['snp.list']]<-snp.list
  result[['x.gls']]<-x.gls
  result[['anno.total']]<-anno.total
  result[['M.anno']]<-M.anno
  result[['Atotal']]<-Atotal
  result[['chrs']]<-chrs
  return(result)
  #saveRDS(result,paste0(out,'.Rdata'))
}

race_eng=c("EUR", "AFR", "EAS")
race_num = c(1,4,5)

for(race in 1:3)
{
  LD.path = sprintf(
    "~/Documents/GWAS/project_4/simulation_new/data/LD_1kg_batch/%s/",
    race_eng[race])
  anno.path = paste0(
    '~/Documents/GWAS/project_4/annotations/baseline_bed_intersect/maf_g_input_53_',
    race_num[race])
  
  snp=NULL
  
  annotation = read.table("~/Documents/GWAS/project_4/annotations/annotation_list.txt")[,1]
  annotation %<>% c("MAF_adj")
  
  out = sprintf('~/Documents/GWAS/project_4/simulation_new/gLDSC/LDSM_1kg/%s/chr1', 
                race_eng[race])
  dir.create(out, recursive = T)
  out_sep = sprintf("%s/block", out)
  
  print('Start...')
  result = mkLDSM(LD.path, anno.path, annotation = annotation, chr_used = 1,
                  out_path = out_sep)
  saveRDS(result, paste0(out, '.Rdata'))
  
  print(paste0('Race ',race,' CHR ', 1, ' done.'))
}