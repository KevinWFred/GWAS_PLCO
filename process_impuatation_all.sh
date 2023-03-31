#!/usr/bin/env bash
#the code is used to process imputed PLCO EUROPEAN data
#ml plink/1.9
plink=/usr/local/apps/plink/1.9/plink
plink2=/usr/local/apps/plink/2.3-alpha/plink2
ml samtools
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
  #if [[ ! -f "$prefix".pgen ]];then
    #read dosage
    $plink2 --vcf $vcffile dosage=DS --max-alleles 2 --snps-only --hwe 0.000001 --maf 0.01 -geno 0.1 --make-pgen --out tmp_s1_$chr --memory 32000 --threads 2
    #$plink2 --vcf $vcffile  --max-alleles 2 --snps-only --hwe 0.000001 --maf 0.001 -geno 0.1 --make-pgen --out tmp_s1_$chr --memory 32000 --threads 2
  
    $plink2 --pfile tmp_s1_$chr --qual-scores $infofile1 7 1 1 --qual-threshold 0.3 --make-pgen --out "$prefix" --memory 32000 --threads 2
  #fi
 
  rm $infofile1
  rm tmp_s*_$chr.*
}

mergechr(){
  local outfolder="$1"
  rm ${outfolder}mergelist.txt
  for chr in {1..22}
  do
    echo ${outfolder}chr${chr}_filtered  >> ${outfolder}mergelist.txt
  done
  
  $plink2 --pmerge-list ${outfolder}mergelist.txt --make-pgen --out ${outfolder}merged_plco --memory 300000
  #plink --bfile ${outfolder}merged --maf 0.05 --make-bed --out ${outfolder}merged_maf05 --memory 30000
  #plink --bfile ${outfolder}merged_plco  --make-bed --out ${outfolder}merged_plco --memory 30000
  #mv ../*.log .
}

echo "work on Oncoarray----------"
infolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Post_Imputation_QCed/Oncoarray/European/"
outfolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Oncoarray/European/"
mkdir -p ${outfolder}log
cd ${outfolder}log
for chr in {1..22}
do
  filtering $chr $infolder $outfolder &
done
wait
mergechr $outfolder

echo "work on OmniX--------------"
infolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Post_Imputation_QCed/OmniX/European/"
outfolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/OmniX/European/"
mkdir -p ${outfolder}log
cd ${outfolder}log
for chr in {1..22}
do
  filtering $chr $infolder $outfolder &
done
wait
mergechr $outfolder

echo "work on Omni25--------------"
infolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Post_Imputation_QCed/Omni25/European/"
outfolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Omni25/European/"
mkdir -p ${outfolder}log
cd ${outfolder}log
for chr in {1..22}
do
  filtering $chr $infolder $outfolder &
done
wait
mergechr $outfolder

echo "work on GSA batch3------------"
infolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Post_Imputation_QCed/GSA/batch3/European/"
outfolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch3/European/"
mkdir -p ${outfolder}log
cd ${outfolder}log
for chr in {1..22}
do
  filtering $chr $infolder $outfolder &
done
wait
mergechr $outfolder

echo "work on GSA batch2------------"
infolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Post_Imputation_QCed/GSA/batch2/European/"
outfolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch2/European/"
mkdir -p ${outfolder}log
cd ${outfolder}log
for chr in {1..22}
do
  filtering $chr $infolder $outfolder &
done
wait
mergechr $outfolder

echo "work on GSA batch1------------"
infolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Post_Imputation_QCed/GSA/batch1/European/"
outfolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch1/European/"
mkdir -p ${outfolder}log
cd ${outfolder}log
for chr in {1..22}
do
  filtering $chr $infolder $outfolder &
done
wait
mergechr $outfolder

echo "work on GSA batch4------------"
infolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Post_Imputation_QCed/GSA/batch4/European/"
outfolder="/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch4/European/"
mkdir -p ${outfolder}log
cd ${outfolder}log
for chr in {1..22}
do
  filtering $chr $infolder $outfolder &
done
wait
mergechr $outfolder

cd /data/DCEGLeiSongData/Kevin/GWAS_PLCO/result
#$plink2 --pfile /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Oncoarray/European/merged_plco --make-bed --out /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Oncoarray/European/merged_plco
#only extract pei samples for bcfs to merge
$plink2 --pfile /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Oncoarray/European/merged_plco --keep /data/DCEGLeiSongData/Kevin/GWAS_PLCO/result/pei_plink2.sample --export bcf vcf-dosage=DS --out /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Oncoarray/European/merged_plco
tabix /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Oncoarray/European/merged_plco.bcf
#$plink2 --pfile /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/OmniX/European/merged_plco --make-bed --out /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/OmniX/European/merged_plco
$plink2 --pfile /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/OmniX/European/merged_plco --keep /data/DCEGLeiSongData/Kevin/GWAS_PLCO/result/pei_plink2.sample --export bcf vcf-dosage=DS --out /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/OmniX/European/merged_plco
#$plink2 --pfile /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/OmniX/European/merged_plco --keep /data/DCEGLeiSongData/Kevin/GWAS_PLCO/result/pei_plink2.sample --recode A-transpose --out /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/OmniX/European/merged_plco
tabix /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/OmniX/European/merged_plco.bcf
#$plink2 --pfile /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Omni25/European/merged_plco --make-bed --out /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Omni25/European/merged_plco
$plink2 --pfile /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Omni25/European/merged_plco  --keep /data/DCEGLeiSongData/Kevin/GWAS_PLCO/result/pei_plink2.sample --export bcf vcf-dosage=DS --out /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Omni25/European/merged_plco
tabix /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Omni25/European/merged_plco.bcf
#$plink2 --pfile /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch1/European/merged_plco --make-bed --out /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch1/European/merged_plco
$plink2 --pfile /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch1/European/merged_plco  --keep /data/DCEGLeiSongData/Kevin/GWAS_PLCO/result/pei_plink2.sample --export bcf vcf-dosage=DS --out /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch1/European/merged_plco
tabix /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch1/European/merged_plco.bcf
#$plink2 --pfile /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch2/European/merged_plco --make-bed --out /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch2/European/merged_plco
$plink2 --pfile /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch2/European/merged_plco  --keep /data/DCEGLeiSongData/Kevin/GWAS_PLCO/result/pei_plink2.sample --export bcf vcf-dosage=DS --out /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch2/European/merged_plco
tabix /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch2/European/merged_plco.bcf
#$plink2 --pfile /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch3/European/merged_plco --make-bed --out /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch3/European/merged_plco
$plink2 --pfile /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch3/European/merged_plco  --keep /data/DCEGLeiSongData/Kevin/GWAS_PLCO/result/pei_plink2.sample --export bcf vcf-dosage=DS --out /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch3/European/merged_plco
tabix /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch3/European/merged_plco.bcf
#$plink2 --pfile /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch4/European/merged_plco --make-bed --out /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch4/European/merged_plco
$plink2 --pfile /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch4/European/merged_plco  --keep /data/DCEGLeiSongData/Kevin/GWAS_PLCO/result/pei_plink2.sample --export bcf vcf-dosage=DS --out /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch4/European/merged_plco
tabix /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch4/European/merged_plco.bcf

#get all data
echo /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Oncoarray/European/merged_plco > ../result/plco_merglist.txt
echo /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/OmniX/European/merged_plco >> ../result/plco_merglist.txt
echo /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Omni25/European/merged_plco >> ../result/plco_merglist.txt
echo /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch1/European/merged_plco >> ../result/plco_merglist.txt
echo /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch2/European/merged_plco >> ../result/plco_merglist.txt
echo /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch3/European/merged_plco >> ../result/plco_merglist.txt
echo /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch4/European/merged_plco >> ../result/plco_merglist.txt

$plink --merge-list ../result/plco_merglist.txt --geno 0.1 --hwe 0.000001 --make-bed --out ../result/plco
#this doesn't work
#$plink2 --pmerge-list ../result/plco_merglist.txt --geno 0.1 --hwe 0.000001 --make-pgen --out ../result/plco
ml samtools
#merge plink2 bcfs
bcftools merge /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Oncoarray/European/merged_plco.bcf /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/OmniX/European/merged_plco.bcf \
/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/Omni25/European/merged_plco.bcf /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch1/European/merged_plco.bcf \
/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch2/European/merged_plco.bcf /data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch3/European/merged_plco.bcf \
/data/DCEGLeiSongData/Kevin/GWAS_PLCO/data/GSA/batch4/European/merged_plco.bcf -o ../result/plco_plink2.bcf
tabix ../result/plco_plink2.bcf
#to deal with GparseWriteByteCt multiallelic dosage size request
bcftools view --types snps -m 2 -M 2  ../result/plco_plink2.bcf -Ou -o ../result/plco_plink2_bi.bcf
$plink2 --bcf ../result/plco_plink2_bi.bcf dosage=DS --maf 0.01 --geno 0.05 --make-pgen --out ../result/plcomaf01


$plink --bfile ../result/plco --maf 0.01 --geno 0.05  --make-bed --out ../result/plcomaf01
#plink --bfile ../result/plcomaf01 --recode A-transpose -out ../result/plcomaf01
#gzip ../result/plcomaf01.traw
#generate PCs
# First, we need to perform prunning
$plink \
    --bfile ../result/plcomaf01 \
    --indep-pairwise 200 50 0.25 \
    --out ../result/plcomaf01
# Then we calculate the first 20 PCs
$plink \
    --bfile ../result/plcomaf01 \
    --extract ../result/plcomaf01.prune.in \
    --pca 20 \
    --threads 3\
    --out ../result/plco

$plink --bfile ../result/plcomaf01 --keep /data/DCEGLeiSongData/Kevin/GWAS_PLCO/result/pei.sample --maf 0.01 --make-bed --out ../result/peiplco    
$plink --bfile ../result/peiplco --recode A-transpose -out ../result/peiplco
gzip -f ../result/peiplco.traw
$plink --bfile ../result/peiplco --freq --out ../result/peiplco    

   
$plink2 --pfile ../result/plcomaf01  --recode A-transpose -out ../result/peiplco_plink2
gzip -f ../result/peiplco_plink2.traw
$plink2 --pfile ../result/plcomaf01 --freq --out ../result/peiplco_plink2  

$plink \
    --bfile ../result/peiplco \
    --indep-pairwise 200 50 0.25 \
    --out ../result/peiplco
# Then we calculate the first 20 PCs
$plink \
    --bfile ../result/peiplco \
    --extract ../result/peiplco.prune.in \
    --pca 20 \
    --out ../result/peiplco