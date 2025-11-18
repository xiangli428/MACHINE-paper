#' Title
#'
#' @param LD.path The LD matrix files to use in calculating LD score matrix. Under this folder all LD matrixs file should be in .hdf5 format.
#' @param anno.path The location of the annotation. The input format here remain the same as .annot in ldsc.
#' @param maf.path The location of the MAF of SNPs. The input format here remain the same as .M_5_50 in ldsc.
#' @param annotation Specify a vector of annotation name to be used in the analysis. 'all' means use all of them.
#' @param snp Specify a list of SNPs 
#' @param MAF MAF for filtering SNPs in the analysis. Defalt set as 0.05. 
#' @param numCores Number of cores used in the analysis.
#'
#' @return A Rdata file contains all LD score matrix information used in gldsc analysis.
#' @export
#'
#' @examples

#LD files:
#LDpath/blocks/LD.mtx.gz; variant_reblock.txt; variant.txt
#LDpath/reblock.txt; variant_reblock.txt.gz
source('~/gLDSC_scripts/gldsc/R/middle_fun_big.R')
#LD.path='~/Documents/GWAS/data/UKBB/5/imputed_genotype_unique_info-0.8_maf-0.001_hwe-1e-6/LD'
#anno.path='/home/r9user9/Documents/GWAS/project_4/annotations/baseline_bed_intersect/g_input_52_5'
#snp=NULL
#annotation='all'
#out='~/Documents/GWAS/data/UKBB/5/imputed_genotype_unique_info-0.8_maf-0.001_hwe-1e-6/LD/chr22/LDSM_chr22_EAS'
mkLDSM <- function(LD.path, anno.path, maf.path=NULL, annotation='all',
                   snp=NULL, MAF=NULL, numCores=30, chr_used=15, out_path){
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
    #ldscs <- anno.files[grep(x = anno.files, pattern = "*.ldscore.gz")]
  }else{
    error.message <- "Ref pannel (ANNOfile) input error"
    stop(error.message)
  }
  
  #maf.files<- list.files(maf.path)
  #if(any(grepl(x = maf.files, pattern = "1000G.mac5eur.*"))){
  #  mafs <- maf.files[grep(x = maf.files, pattern = "1000G.mac5eur.*")]
  #}else{
  #  error.message <- "Ref pannel (MAFfiles) input error"
  #  stop(error.message)
  #}
  
  #if(length(chrs)==length(annos)&length(chrs)==length(mafs)){
  #  chr.length<-length(chrs)
  #}else{
  #  error.message <- "Ref pannel (number of chr) input error"
  #  stop(error.message)
  #}
  
  
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
    block.info<-block.info[1:(length(block.info)-2)]
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

race_eng=c("european", "afram", "asian")
race_num = c(1,4,5)
chr_order=c(1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,3,4,5,6,7,8,9)

#anno1=as.vector(read.table('/home/r9user9/Documents/GWAS/project_4/simulation/data/annotation_list_european.txt')[,1])
#anno2=as.vector(read.table('/home/r9user9/Documents/GWAS/project_4/simulation/data/annotation_list_asian.txt')[,1])

for (race in c(3)){
  for (chrom in c(19:22)){
    LD.path=paste0('~/Documents/GWAS/data/UKBB/',race_num[race],'/imputed_genotype_unique_info-0.8_maf-0.001_hwe-1e-6/LD/')
    anno.path=paste0('~/Documents/GWAS/project_4/annotations/baseline_bed_intersect/maf_g_input_53_',race_num[race])
    snp=NULL
    
    annotation='all'
    # out=paste0('/home/r14user2/Documents/GWAS/project_4/annotations/LDSM/',
    #            race_eng[race], "_chr", chr_order[chrom])
    out = sprintf('/home/r14user2/Documents/GWAS/project_4/annotations/LDSM/%s/chr%d', 
                  race_eng[race], chr_order[chrom])
    dir.create(out, recursive = T)
    out_sep = sprintf("%s/block", out)
    # out_sep=paste0('/home/r14user2/Documents/GWAS/project_4/annotations/LDSM/sepLD_',race_eng[race],
    #            "_chr", chr_order[chrom])
    #rm(result)
    print('Start...')
    result=mkLDSM(LD.path,anno.path,chr_used=chrom,out_path=out_sep)
    saveRDS(result,paste0(out,'.Rdata'))
    
    print(paste0('Race ',race,' CHR ', chrom, ' done.'))
    
  }
}





#result=mkLDSM(LD.path,anno.path)
#saveRDS(result,paste0(out,'.Rdata'))
