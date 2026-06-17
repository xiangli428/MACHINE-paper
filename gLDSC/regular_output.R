source('~/gLDSC_scripts/gldsc/R/middle_fun_big.R')
for(i in 1:3){
race=c('EUR','AFR','EAS')
race_num=c(1,4,5)
result_path=paste0('~/Documents/GWAS/project_4/real_data/scz2022/gLDSC_results/',race[i],'/chr')
anno_path=paste0('~/Documents/GWAS/project_4/annotations/baseline_bed_intersect/maf_g_input_53_',race_num[i])
out=paste0('~/Documents/GWAS/project_4/real_data/scz2022/gLDSC_results/')

#sum up all the lefts and rights
gls.left<-gls.right<-0
Atotal=NULL
for (chr in 1:22){
temp_data=readRDS(paste0(result_path,'/chr',chr,'.Rdata'))
temp_anno=read.table(paste0(anno_path,'/Amatrix.',chr,'.annot'),header=T)
Atotal=rbind(Atotal,temp_anno[,c(-1,-2)])

for(batch in 1:length(temp_data)){
        gls.left<-gls.left+temp_data[[batch]][[1]]
        gls.right<-gls.right+temp_data[[batch]][[2]]
}
}
Atotal=as.matrix(Atotal)
M=nrow(Atotal)
anno.total<-apply(Atotal, 2, sum)
M.anno<-anno.total/M

#0624: dim of left right doesn't match with annotation file
#result
result.one<-gls.estimator(left = gls.left,right = gls.right,
                            A=Atotal, M.anno = M.anno, anno.total = anno.total, N = 100000,intercept = T)

saveRDS(result.one,paste0(out,race[i],'_enrichment_result.Rdata'))
print(paste0(i,' Done'))
}
