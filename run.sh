#!/usr/bin/env bash
#to split genomewide genotype into multiple jobs (each job has 10000 SNPs) and submit

basedir=/data/BB_Bioinformatics/Kevin/GWAS_PLCO/
cd ${basedir}code

#split
# if [ -d ${basedir}result/splited ]
# then 
# rm -r ${basedir}result/splited
# fi

if [ ! -d ${basedir}result/splited ]
then
mkdir -p ${basedir}result/splited
cd ${basedir}result/splited
zcat ${basedir}result/prostatemaf01.traw.gz | parallel --header : --pipe -N10000 'cat | gzip > processed.traw__{#}.gz'
fi

#create swarm file
if [ -f ${basedir}code/run.swarm ]
then
rm  ${basedir}code/run.swarm
fi

rscript=" ${basedir}code/gwas_survival.R"

splited_file_dir="${basedir}result/splited"
out_dir="${basedir}result/splited"

opt="cancerdeath"
ls --color=never $splited_file_dir/*.gz | sort -V | while read ifile
do
ofile=`basename $ifile | sed 's/.gz/.txt/g' -`
printf "Rscript $rscript $ifile ${out_dir}/$ofile $opt \n" >> ${basedir}code/run.swarm
#echo $ifile
#echo $ofile
done

#rm ${basedir}result/processed.traw__*.txt
#submit job
if [ -d ${basedir}logs ]
then
rm -r ${basedir}logs
fi

mkdir -p ${basedir}code/logs
cd ${basedir}code/logs

swarm -f ${basedir}code/run.swarm -g 64 --module R/4.2.0 --time=10:00:00 --gres=lscratch:64

#for all death
opt="alldeath"
rm ${basedir}code/run_alldeath.swarm
ls --color=never $splited_file_dir/*.gz | sort -V | while read ifile
do
ofile=`basename $ifile | sed 's/.gz/_alldeath.txt/g' -`
printf "Rscript $rscript $ifile ${out_dir}/$ofile $opt \n" >> ${basedir}code/run_alldeath.swarm
#echo $ifile
#echo $ofile
done
swarm -f ${basedir}code/run_alldeath.swarm -g 64 --module R/4.2.0 --time=10:00:00 --gres=lscratch:64