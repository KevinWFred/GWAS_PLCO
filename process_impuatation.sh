#!/usr/bin/env bash
#the code is used to process imputed PLCO EUROPEAN data
ml plink/1.9

#QC for imputed data
filtering(){
  local chr="$1"
  echo $chr
  local infolder="$2"
  local outfolder="$3"
  prefix="$outfolder"chr${chr}_filtered
  #vcf file
  vcffile="$infolder"chr${chr}-filtered.dose.vcf.gz
  #info file
  infofile="$infolder"chr${chr}-filtered.info.gz
  infofile1="$infolder"chr${chr}-filtered.info
  if [[ ! -f $infofile1 ]];then
    gunzip -c $infofile > $infofile1
  fi
  #plink --vcf $vcffile  --maf 0.05 --biallelic-only strict --snps-only --hwe 0.00001 --geno 0.05 --make-bed --out tmp_s1_$chr --memory 6400
  #plink --vcf $vcffile  --maf 0.05 --hwe 0.000001 --geno 0.05 --make-bed --out tmp_s1_$chr --memory 6400
  plink --vcf $vcffile --biallelic-only --snps-only --hwe 0.000001 -geno 0.1 --make-bed --out tmp_s1_$chr --memory 6400 
  plink --bfile tmp_s1_$chr --qual-scores $infofile1 7 1 1 --qual-threshold 0.3 --make-bed --out "$prefix" --memory 6400
  
  #If binary merging fails because at least one variant would have more than two alleles, a list of offending variant(s) will be written to plink.missnp.multi-allelic
  #plink --bfile tmp_s1_$chr  --bmerge tmp_s1_$chr --merge-mode 6 --out tmp_$chr   --memory 6400
  #if [[ -f tmp_$chr.missnp ]];then
  #  plink --bfile tmp_s1_$chr  --exclude tmp_$chr.missnp --make-bed --out tmp_s1_$chr --memory 6400
  #fi
  #When multiple variants share the same bp coordinate and allele codes, it is likely that they are not actually distinct, and the duplicates should be merged or removed. (In fact, some tools, such as BEAGLE 4, require such duplicate variants to be merged/removed during preprocessing.) --list-duplicate-vars identifies them, and writes a report to plink.dupvar.
  #plink --bfile tmp_s1_$chr  --list-duplicate-vars --out tmp_$chr --memory 6400
  #plink --bfile tmp_s1_$chr  --exclude tmp_$chr.dupvar --make-bed --out tmp_s2_$chr --memory 6400
  #plink --bfile tmp_s2_$chr  --maf 0.05 --biallelic-only strict --snps-only --hwe 0.00001 --geno 0.05 --qual-scores $infofile1 7 1 1 --qual-threshold 0.3 --make-bed --out $prefix --memory 6400
  #filter samples
  #plink --bfile $prefix --make-bed --keep /data/BB_Bioinformatics/Kevin/GWAS_PLCO/result_plco_prcasample_pink.txt --out "$prefix" --memory 6400
  #create dosage
  #plink --bfile $prefix --recode A-transpose --out "$prefix" --memory 6400

  rm $infofile1
  rm tmp_s*_$chr.bed
  rm tmp_s*_$chr.bim
  rm tmp_s*_$chr.fam
  rm tmp_s*_$chr.nosex 
}

mergechr(){
  local outfolder="$1"
  rm ${outfolder}mergelist.txt
  for chr in {1..22}
  do
    echo ${outfolder}chr${chr}_filtered  >> ${outfolder}mergelist.txt
  done
  #plink --merge-list ${outfolder}mergelist.txt --make-bed --out ${outfolder}merged --memory 30000
  #plink --bfile ${outfolder}merged --maf 0.05 --make-bed --out ${outfolder}merged_maf05 --memory 30000
  plink --bfile ${outfolder}merged --keep /data/BB_Bioinformatics/Kevin/GWAS_PLCO/result/Prostate_plinksample.txt --make-bed --out ${outfolder}merged_prostate --memory 30000
  mv ../*.log .
}
infolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Post_Imputation_QCed/Oncoarray/European/"
outfolder="/data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/Oncoarray/European/"
mkdir -p ${outfolder}log
cd ${outfolder}log
for chr in {1..22}
do
  #filtering $chr $infolder $outfolder &
done
#wait
mergechr $outfolder


infolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Post_Imputation_QCed/OmniX/European/"
outfolder="/data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/OmniX/European/"
mkdir -p ${outfolder}log
cd ${outfolder}log
for chr in {1..22}
do
  #filtering $chr $infolder $outfolder &
done
#wait
mergechr $outfolder

infolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Post_Imputation_QCed/Omni25/European/"
outfolder="/data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/Omni25/European/"
mkdir -p ${outfolder}log
cd ${outfolder}log
for chr in {1..22}
do
  #filtering $chr $infolder $outfolder &
done
#wait
mergechr $outfolder

infolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Post_Imputation_QCed/GSA/batch3/European/"
outfolder="/data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/GSA/batch3/European/"
mkdir -p ${outfolder}log
cd ${outfolder}log
for chr in {1..22}
do
  #filtering $chr $infolder $outfolder &
done
#wait
mergechr $outfolder

infolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Post_Imputation_QCed/GSA/batch2/European/"
outfolder="/data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/GSA/batch2/European/"
mkdir -p ${outfolder}log
cd ${outfolder}log
for chr in {1..22}
do
  #filtering $chr $infolder $outfolder &
done
#wait
mergechr $outfolder

infolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Post_Imputation_QCed/GSA/batch1/European/"
outfolder="/data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/GSA/batch1/European/"
mkdir -p ${outfolder}log
cd ${outfolder}log
for chr in {1..22}
do
  #filtering $chr $infolder $outfolder &
done
#wait
mergechr $outfolder

infolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Post_Imputation_QCed/GSA/batch4/European/"
outfolder="/data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/GSA/batch4/European/"
mkdir -p ${outfolder}log
cd ${outfolder}log
for chr in {1..22}
do
  #filtering $chr $infolder $outfolder &
done
mergechr $outfolde

#get prostate data
echo /data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/Oncoarray/European/merged_prostate > ../result/prostate_merglist.txt
echo /data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/OmniX/European/merged_prostate >> ../result/prostate_merglist.txt
echo /data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/Omni25/European/merged_prostate >> ../result/prostate_merglist.txt
echo /data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/GSA/batch1/European/merged_prostate >> ../result/prostate_merglist.txt
echo /data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/GSA/batch2/European/merged_prostate >> ../result/prostate_merglist.txt
echo /data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/GSA/batch3/European/merged_prostate >> ../result/prostate_merglist.txt
echo /data/BB_Bioinformatics/Kevin/GWAS_PLCO/data/GSA/batch4/European/merged_prostate >> ../result/prostate_merglist.txt

plink --merge-list ../result/prostate_merglist.txt --geno 0.1 --hwe 0.000001 --make-bed --out ../result/prostate --memory 30000
plink --bfile ../result/prostate --maf 0.01 --geno 0.05 --keep ../result/Prostate_plinksel6997sample.txt --make-bed --out ../result/prostatemaf01
#generate PCs
# First, we need to perform prunning
plink \
    --bfile ../result/prostatemaf01 \
    --indep-pairwise 200 50 0.25 \
    --out ../result/prostatemaf01
# Then we calculate the first 20 PCs
plink \
    --bfile ../result/prostatemaf01 \
    --extract ../result/prostatemaf01.prune.in \
    --pca 20 \
    --out ../result/prostate
    