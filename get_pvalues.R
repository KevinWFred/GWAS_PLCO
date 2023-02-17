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
res[which(res$P<5e-8)]
qqplot(res_alldeath$P,title="")


plot(-log10(res_alldeath$P),-log10(res$P),xlab="All death -log10(P)",ylab="Prostate cancer death -log10(P)",cex.lab=1.2,cex.axis=1.2)
plot(res_alldeath$Beta,res$Beta,xlab="All death effect size",ylab="Prostate cancer death effect size", cex.lab=1.2,cex.axis=1.2)
