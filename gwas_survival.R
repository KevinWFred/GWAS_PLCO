#!/usr/bin/env Rscript
#survival analysis
library(data.table)
library("survival")
library("survminer")
setwd("/data/BB_Bioinformatics/Kevin/GWAS_PLCO/code")
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) stop("Two inputs are required: file names of genotype and output")
SNP_TRAW_FILE = args[1]
OUT_FILE = args[2]
opt = args[3]
gwas_surv<- function(snp.traw.file = "/data/BB_Bioinformatics/Kevin/GWAS_PLCO/result/splited/processed.traw__1.gz", 
                       out.file = "/data/BB_Bioinformatics/Kevin/GWAS_PLCO/result/splited/pradeath__1.txt",
                       opt="cancerdeath")
{
  print(snp.traw.file)
  print(out.file)
  print(opt)
  #phenotype data
  load("../result/pheno.RData")
  geno=fread(snp.traw.file)
  #geno=fread("../result/prostatemaf01.traw.gz",nrows = 10)
  rownames(geno)=geno$SNP
  geno=geno[,7:ncol(geno)]  
  tmp=unlist(strsplit(colnames(geno),"_"))
  colnames(geno)=tmp[seq(1,length(tmp),2)]
  if (any(colnames(geno)!=pheno$plco_id)) stop("sample orders of geno and pheno not consistent!")
  out=data.frame(snp=rownames(geno),P=rep(NA,nrow(geno)),Beta=NA,SE=NA)
  
  for (i in 1:nrow(geno))
  {
    dat=cbind(pheno,snp=unlist(geno[i,]))
    if (opt=="cancerdeath")
    {
      m1 = tryCatch(
        expr = {
          coxph(Surv(time, status) ~ snp+psa_grp+agelevel+primary_trtp+PC1+PC2+PC3+PC4+PC5, data = dat)
        },
        error = function(e){ 
          return(NULL)
        }
      )
    }else
    {
      m1 = tryCatch(
        expr = {
          coxph(Surv(time, status1) ~ snp+psa_grp+agelevel+primary_trtp+PC1+PC2+PC3+PC4+PC5, data = dat)
        },
        error = function(e){ 
          return(NULL)
        }
      )
    }
    if (!is.null(m1))
    {
      ctable = tryCatch(
        expr= {coef(summary(m1))
        },
        error = function(e){ 
          return(NULL)
        }
      )
      if (!is.null(ctable))
      {
        
        out[i, "SE"] = ctable["snp", "se(coef)"]
        out[i, "P"] = ctable["snp", "Pr(>|z|)"]
        out[i, "Beta"] = ctable["snp", "coef"]
      }
    }
  }
  write.table(out,file=out.file,row.names = F,sep="\t",quote=F)
  print("done")
}
gwas_surv(snp.traw.file = args[1], 
                       out.file = args[2],
                       opt=args[3])