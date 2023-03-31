#!/usr/bin/env Rscript
#the code is used to read phenotypes

setwd("/data/BB_Bioinformatics/Kevin/GWAS_PLCO/code")
library(gtsummary)

#genotype data from imputation
famGSA1=read.table("/data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/GSA/batch1/European/merged_prostate.fam") #505
famGSA2=read.table("/data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/GSA/batch2/European/merged_prostate.fam") #523
famGSA3=read.table("/data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/GSA/batch3/European/merged_prostate.fam") #505
famGSA4=read.table("/data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/GSA/batch4/European/merged_prostate.fam") #491
famomni25=read.table("/data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/Omni25/European/merged_prostate.fam") #4496
famomnix=read.table("/data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/OmniX/European/merged_prostate.fam") #26
famonco=read.table("/data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/Oncoarray/European/merged_prostate.fam") #942
allfam=rbind(famGSA1,famGSA2,famGSA3,famGSA4,famomni25,famonco)
allfam$platform=c(rep("GSA1",nrow(famGSA1)),rep("GSA2",nrow(famGSA2)),rep("GSA3",nrow(famGSA3)),rep("GSA4",nrow(famGSA4)),
                  rep("Omni25",nrow(famomni25)),rep("Oncoarray",nrow(famonco)))

oldpheno=read.table("../data/PLCO_survival_PCA.txt",header=T)
library(sas7bdat)
library(dplyr)
alldata=read.sas7bdat("../data/plco_845_pros_aug21_101421.sas7bdat")
table(alldata$j_ph_pros_trial)
table(alldata$j_ph_any_trial)
table(alldata$j_pros_cancer)
table(alldata$j_pros_cancer_first)
table(alldata$mortality_exitstat,alldata$pros_death_stat)
table(alldata$pros_death_stat,alldata$is_dead)
pros_samples=alldata$plco_id[which(alldata$j_pros_cancer_first==1 & alldata$j_ph_any_trial==0 & alldata$mortality_exitdays > alldata$j_pros_cancer_diagdays)]
pros_samples=intersect(pros_samples,allfam$V1)

peidat=read.csv("../data/pros_data_7parts.csv")
plcofam=read.table("/data/DCEGLeiSongData/Kevin/GWAS_PLCO/result/plcomaf01.fam",sep=" ")
tmp=unlist(strsplit(plcofam$V2,"_"))
plcofam$V2=tmp[seq(1,length(tmp),2)]
sum(peidat$plco_id %in% plcofam$V2)
write.table(plcofam,file="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/result/plcomaf01.fam",sep=" ",col.names=F,row.names=F,quote=F)
tmp=data.frame(FID=0,IID=peidat$plco_id)
write.table(tmp,file="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/result/pei.sample",sep=" ",col.names=F,row.names=F,quote=F)
tmp=data.frame(sample=paste0(peidat$plco_id,"_",peidat$plco_id))
write.table(tmp,file="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/result/pei_plink2.sample",sep=" ",col.names=F,row.names=F,quote=F)

sum(alldata$plco_id %in% peidat$plco_id)
tmp=intersect(alldata$plco_id,peidat$plco_id)

pheno=alldata[,c("plco_id", "pros_death_stat", "j_ph_any_trial","j_pros_cancer_first","agelevel","center", "primary_trtp", "dth_days",
                 "mortality_exitdays", "j_pros_cancer_diagdays","j_dx_psa","j_dx_psa_gap","race7")]
pheno=pheno %>%
  mutate(status=ifelse(pros_death_stat %in% c(0,2),1,2))
pheno$time=pheno$mortality_exitdays-pheno$j_pros_cancer_diagdays

pheno=pheno[match(pros_samples,pheno$plco_id),]
notincluded=pheno[!pheno$plco_id %in% oldpheno$plco_id,] #10, they all have time=0


#j_dx_psa: PSA Closest To Diagnosis (ng/mL)
pheno=alldata[,c("plco_id", "pros_death_stat", "j_ph_any_trial","j_pros_cancer_first","agelevel","center", "primary_trtp", "dth_days",
                 "mortality_exitstat","mortality_exitdays", "j_pros_cancer_diagdays","j_dx_psa","j_dx_psa_gap","race7")]
pheno=pheno %>%
  mutate(status=ifelse(pros_death_stat %in% c(0,2),1,2)) #1 died of prostate
pheno$time=pheno$mortality_exitdays-pheno$j_pros_cancer_diagdays
pheno=pheno[pheno$plco_id %in% pros_samples,]
pheno$batch=allfam$platform[match(pheno$plco_id,allfam$V1)]
pheno$batch=as.factor(pheno$batch)
pheno$psa_grp=cut(pheno$j_dx_psa,breaks=c(0,4,10,20,100000),labels=c("1: 0-3.9","2: 4-9.9","3: 10-19.9","4: 20+"),include.lowest = T,right=F)
table(pheno$pros_death_stat,pheno$mortality_exitstat)
pheno$agelevel=factor(pheno$agelevel,labels=c("0_59","60_64","65_69","GE70"))
pheno$center=factor(pheno$center,labels=c("UColorado","GeorgetownU","PacificHealth","HenryFordHealth","UMinnesota","WashingtonU",
                                         "UPittsburgh","UUtah","MarshfieldClinic","UAlabama"))
pheno$primary_trtp=factor(pheno$primary_trtp,labels = c("NoCancer","Prostatectomy","Radiation","RadiationandHormone","Hormone","OtherAblative","NoKnownTreatmentWithCurativeIntent","notcollected"))
pheno$primary_trtp[as.character(pheno$primary_trtp)=="notcollected"]=NA
pheno$primary_trtp=droplevels(pheno$primary_trtp)
#      1    2    3    4
# 0    0 3610  506   37
# 1  427    0    0    0
# 2 2139    0    0    0
table(pheno$status,pheno$center)
table(pheno$center)
table(pheno$primary_trtp,useNA="ifany")
table(pheno$status,pheno$primary_trtp,useNA="ifany")
table(pheno$psa_grp,useNA="ifany")
table(pheno$status,pheno$psa_grp,useNA="ifany")
table(pheno$batch)
table(pheno$status,pheno$batch)
idx=match(pheno$plco_id,alldata$plco_id)
table(alldata$mortality_exitstat[idx],alldata$pros_death_stat[idx])
table(pheno$status,pheno$pros_death_stat)
#count all death events
pheno$status1=1
pheno$status1[which(pheno$pros_death_stat %in% c(1,2))]=2
table(pheno$status1)
save(pheno,file="../result/pheno.RData")
write.table(data.frame(FID=pheno$plco_id,IID=pheno$plco_id),file="../result/Prostate_plinksample.txt",row.names = F,col.names = F,sep="\t",quote=F)
#check covariates
library("survival")
library("survminer")
fit0=coxph(Surv(time, status) ~ 1, data = pheno)
fit=coxph(Surv(time, status) ~ agelevel, data = pheno)
fit1=coxph(Surv(time=j_pros_cancer_diagdays,time2=mortality_exitdays, event=status) ~agelevel,data=pheno)
anova(fit0,fit,test = "Chisq")

fit=coxph(Surv(time, status) ~ center, data = pheno)
fit=coxph(Surv(time, status) ~ primary_trtp, data = pheno)
fit=coxph(Surv(time, status) ~ psa_grp, data = pheno)
fit=coxph(Surv(time, status) ~ batch, data = pheno)%>% 
  tbl_regression(exp = TRUE) 
idx=!is.na(pheno$primary_trtp)
fita=coxph(Surv(time, status) ~ agelevel, data = pheno[idx,])
fitb=coxph(Surv(time, status) ~ agelevel+primary_trtp, data = pheno[idx,])
anova(fita,fitb,test = "Chisq")

fit=coxph(Surv(time, status) ~ psa_grp+agelevel+primary_trtp+batch+center+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data =pheno)
summary(fit)

#PCA plots
pcadat=read.table("../result/prostate.eigenvec")
rownames(pcadat)=pcadat$V1
idx=match(pheno$plco_id,pcadat$V1)
pcadat1=pcadat[idx,3:ncol(pcadat)]
colnames(pcadat1)=paste0("PC",1:20)
pheno=cbind(pheno,pcadat1)
save(pheno,file="../result/pheno.RData")
tmp=cor(as.matrix(oldpheno[,12:21]),as.matrix(pcadat1[,1:10]))
# idx=match(nooldphenosamples,pcadat$V1)
# pcadat$col="green"
# pcadat$col[idx]="red"
# library(scales)
# 
# png(filename = "../result/PC12.png", width = 8, height = 8, units = "in",res=300)
# plot(pcadat$V3,pcadat$V4,xlab="PC1",ylab="PC2",col=alpha(pcadat$col,0.9),cex.axis=1.3,cex.lab=1.3)
# points(pcadat$V3[idx],pcadat$V4[idx],col="red")
# legend("bottomleft",legend = c("included","not included"),col=c("green","red"),pch=1,cex=1.2)
# dev.off()
# png(filename = "../result/PC13.png", width = 8, height = 8, units = "in",res=300)
# plot(pcadat$V3,pcadat$V5,xlab="PC1",ylab="PC3",col=alpha(pcadat$col,0.9),cex.axis=1.3,cex.lab=1.3)
# points(pcadat$V3[idx],pcadat$V5[idx],col="red")
# legend("bottomleft",legend = c("included","not included"),col=c("green","red"),pch=1,cex=1.2)
# dev.off()
# png(filename = "../result/PCAeigenvalue.png", width = 8, height = 8, units = "in",res=300)
# tmp=read.table("../result/prostate.eigenval")
# plot(tmp$V1,ylab="Eigen value",xlab="PCs",cex.axis=1.3,cex.lab=1.3)
# dev.off()

pcadat=read.table("/data/DCEGLeiSongData/Kevin/GWAS_PLCO/result/peiplco.eigenvec")
colnames(pcadat)[3:22]=paste0("PC",1:20)
png(filename="../result/peiplco_pca.png",res=100)
plot(pcadat[, c("PC1", "PC2", "PC3", "PC4", "PC5")])
dev.off()
pcabim=read.table("/data/DCEGLeiSongData/Kevin/GWAS_PLCO/result/peiplco.bim")
library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)

## this step is optional!
## here we just simplify the names of the objects, making the code neater
genome <- BSgenome.Hsapiens.UCSC.hg38
all_snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38

## By default the genome we're using follows the UCSC convention for
## naming chromosome e.g. "chr8".  This step changes that to match our
## SNP data which uses NCBI naming e.g. "8"
seqlevelsStyle(genome) <- "NCBI"

## construct a GPos object containing all the positions we're interested in
positions <- GPos(seqnames = pcabim$V1, pos =pcabim$V4)

## query the genome with out positions
my_snps <- snpsByOverlaps(all_snps, positions, genome = genome)

## this gives us a GPos object
my_snps = as.data.frame(my_snps)
tmp1=paste0(pcabim$V1,"_",pcabim$V4)
tmp2=paste0(my_snps$seqnames,"_",my_snps$pos)
tmp3=intersect(tmp1,tmp2)
idx1 = match(tmp3,tmp1)
idx2=match(tmp3,tmp2)

pcabim$V2[idx1]=my_snps$RefSNP_id[idx2]
write.table(pcabim,file="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/result/peiplco.bim",sep=" ",row.names=F,col.names=F,quote=F)


# #process LeiPLCO data
# tmp=read.table("/data/DCEGLeiSongData/Kevin/GWAS_PLCO/result/plcomaf01.pvar")
# write.table(tmp$V3,file="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/result/plcoma01_bivariate.snp",row.names = F,col.names = F,quote=F)
