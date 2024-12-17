#!/bin/sh
#PBS -V
#PBS -j oe
#PBS -l nodes=1:ppn=1,mem=10gb,walltime=120:00:00
#PBS -q batch
SAMPLE1=$1
SAMPLE2=$2

config_path=$3
source ${config_path}/config.sh

echo ${SAMPLE1}_${SAMPLE2}

${gatk4} --java-options "-Xmx10g" Mutect2 \
-R ${seq_ref}/GRCh37.fa \
-I ${bam_path}/${SAMPLE1}-recal.bam \
-tumor ${SAMPLE1} \
--alleles ${vcf_path}/Combine_${SAMPLE2}_forGGA.vcf.gz \
--native-pair-hmm-threads 10 \
-L ${vcf_path}/Combine_${SAMPLE2}_forGGA.vcf.gz \
-O ${vcf_path}/${SAMPLE1}_${SAMPLE2}_GGA.vcf.gz \
-genotype-filtered-alleles true \
-bamout ${vcf_path}/${SAMPLE1}_${SAMPLE2}_GGA.bam \
2>  ${log_path}/${SAMPLE1}_${SAMPLE2}_GGA.log



