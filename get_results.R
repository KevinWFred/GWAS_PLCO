#!/usr/bin/env Rscript

setwd("/data/BB_Bioinformatics/Kevin/GWAS_PLCO/code")

get_pvalues=function(opt="cancerdeath")
{
  res=NULL
  for (i in 1:708)
  {
    if (opt=="cancerdeath")
    {
      tmp=read.table(paste0("../result/splited/","processed.traw__",i,".txt"),header=T)
    }else
    {
      tmp=read.table(paste0("../result/splited/","processed.traw__",i,"_alldeath.txt"),header=T)
    }
    
    res=rbind(res,tmp)
  }
  #res=res[!is.na(res$P),]
  return(res)
}

res=get_pvalues()
res_alldeath=get_pvalues(opt="alldeath")

qqplot=function(pvalue=NULL,main="",xlim=NULL,ylim=NULL,title="")
  
{
  
  pvalue=pvalue[!is.na(pvalue)]
  n=length(pvalue)
  par(mar=c(5,5,2,1))
  plot(-log((1:n)/n,base=10),-log(pvalue[order(pvalue)],base=10),xlab="Expected p-value (log base 10)",
       
       ylab="Observed p-value (log base 10)",main=main,xlim=xlim,ylim=ylim,cex.lab=1.2,cex.axis=1.2)
  title(main=title,cex=1.2)
  abline(0,1,lty=2)
  chisq <- qchisq(1-pvalue,1)
  lambda = median(chisq)/qchisq(0.5,1)
  print(lambda)
  
}
png("../result/prostate_qqplot.png",res=100)
qqplot(res$P,title="Prostate cancer death")
dev.off()
res[which(res$P<5e-8),]
png("../result/alldeath_qqplot.png",res=100)
qqplot(res_alldeath$P,title="All death")
dev.off()
res_alldeath[which(res_alldeath$P<5e-8),]

tmp=unlist(strsplit(res$snp,":"))
res$CHR=tmp[seq(1,length(tmp),4)]
res$CHR=as.integer(gsub("chr","",res$CHR))
res$BP=as.integer(tmp[seq(2,length(tmp),4)])
res$A1=tmp[seq(3,length(tmp),4)]
res$A2=tmp[seq(4,length(tmp),4)]
res_alldeath$CHR=res$CHR
res_alldeath$BP=res$BP
#manhattan plot
library(qqman)
png("../result/manhattan.png",res=100, width=1200)
manhattan(res, chr="CHR", bp="BP", snp="snp", p="P", suggestiveline = -log10(5e-5),col=c("navyblue", "springgreen3"),cex.axis=1.2,cex.lab=1.2)
dev.off()
png("../result/alldeath_manhattan.png",res=100, width=1200)
manhattan(res_alldeath, chr="CHR", bp="BP", snp="snp", p="P", suggestiveline = -log10(5e-5),col=c("navyblue", "springgreen3"),cex.axis=1.2,cex.lab=1.2)
dev.off()

#get rsid
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)

## this step is optional!
## here we just simplify the names of the objects, making the code neater
genome <- BSgenome.Hsapiens.UCSC.hg38
all_snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
seqlevelsStyle(genome) <- "NCBI"

## construct a GPos object containing all the positions we're interested in
positions <- GPos(seqnames = res$CHR, pos = res$BP)

## query the genome with out positions
my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)

## this gives us a GPos object
my_snps = as.data.frame(my_snps)

idx = match(paste0(res$CHR,":",res$BP),paste0(my_snps$seqnames,":",my_snps$pos))
res$rsid=my_snps$RefSNP_id[idx]
save(res,file="../result/prostate_death_gwas_batch.RData")
load("/data/DCEGLeiSongData/Kevin/HBVloadGwas/result/refgeneshg38codinggenes.RData") #allgenes
topsnps=res[res$P<5e-5,]
topsnps$gene=NA
topsnps$closestgene=NA
topsnps$gene_dist=NA
gr_topsnps=GRanges(seqnames = topsnps$CHR,ranges = IRanges(topsnps$BP,width = 1))
gr_allgenes=GRanges(seqnames = allgenes$chr,ranges = IRanges(start=allgenes$start,end=allgenes$end))
for (i in 1:nrow(topsnps))
{
  tmp=distance(gr_allgenes,gr_topsnps[i])
  idx=which.min(tmp)
  if (tmp[idx]>0)
  {
    if (allgenes$start[idx]>topsnps$BP[i])
    {
      topsnps$gene[i]=paste0(allgenes$gene[idx-1],";",allgenes$gene[idx])
    }
    if (allgenes$end[idx]<topsnps$BP[i])
    {
      topsnps$gene[i]=paste0(allgenes$gene[idx],";",allgenes$gene[idx+1])
    }
  }else
  {
    topsnps$gene[i]=allgenes$gene[idx]
  }
  topsnps$closestgene[i]=allgenes$gene[idx]
  topsnps$gene_dist[i]=min(tmp,na.rm = T)
}
topsnps=topsnps[order(topsnps$P),]
head(topsnps[,c("rsid","CHR","BP","A1","A2","Beta","SE","P","gene")])
