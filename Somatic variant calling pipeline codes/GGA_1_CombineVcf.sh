#!/bin/sh
#PBS -V
#PBS -j oe
#PBS -l nodes=1:ppn=1,mem=10gb,walltime=120:00:00
#PBS -q batch
export normal=$1
config_path=$2
source ${config_path}/config.sh

set -e

rm -rf  ${tmp_path}/Combine_vcf_${normal}.list

cat ${config_path}/tumor_normal.list | grep $normal | \
awk '{print vcf_path"/"$1"_"$2"_QC.vcf"}' vcf_path=${vcf_path} \
> ${tmp_path}/Combine_vcf_${normal}.list
cd ${vcf_path}

cat ${config_path}/tumor_normal.list | grep ${normal}


for samples in ` cat ${config_path}/tumor_normal.list | grep ${normal} | tr '\t' ',' ` 
do
echo $samples
Tumor=` echo $samples | awk -F, '{print$1}' `
Normal=` echo $samples | awk -F, '{print$2}' `

file=${vcf_path}/${Tumor}_${Normal}_QC.vcf

tumor_pos=` cat $file | grep CH | awk '{for(i=1;i<=NF;i++){if($i==tumor){print i}}}' tumor=$Tumor `
awk -F"\t" '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$tumor_pos}'  tumor_pos=$tumor_pos $file | \
awk '{if($1!~"X" && $1!~"Y" )print}' \
> ${tmp_path}/${Tumor}_${Normal}_combine.vcf

bgzip -f ${tmp_path}/${Tumor}_${Normal}_combine.vcf
tabix -f ${tmp_path}/${Tumor}_${Normal}_combine.vcf.gz
done

cat ${config_path}/tumor_normal.list | grep $normal | \
awk '{print tmp_path"/"$1"_"$2"_combine.vcf.gz"}' tmp_path=${tmp_path} \
> ${tmp_path}/Combine_vcf_${normal}.list
bcftools merge -l ${tmp_path}/Combine_vcf_${normal}.list --force-samples \
 -Oz -o ${vcf_path}/Combine_${normal}_forGGA.vcf.gz
tabix -f ${vcf_path}/Combine_${normal}_forGGA.vcf.gz

zcat ${vcf_path}/Combine_${normal}_forGGA.vcf.gz | grep -v "#" | awk '{OFS="\t"}{print $1,$2,$2}' >  ${tmp_path}/Combine_${normal}_GGA_site.bed

