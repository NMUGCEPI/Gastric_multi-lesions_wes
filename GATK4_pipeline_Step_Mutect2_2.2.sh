#!/bin/sh
#PBS -V
#PBS -j oe
#PBS -l nodes=1:ppn=15,mem=20gb,walltime=120:00:00
#PBS -q batch

SAMPLE1=$1
SAMPLE2=$2
config_path=$3
source ${config_path}/config.sh

echo ${sample}

${gatk4} --java-options "-Xmx10g" CalculateContamination \
-I ${vcf_path}/${SAMPLE1}_getpileupsummaries.table \
-tumor-segmentation ${vcf_path}/${SAMPLE1}_segments.table \
-matched ${vcf_path}/${SAMPLE2}_getpileupsummaries.table \
-O ${vcf_path}/${SAMPLE1}_contamination.table \
2>  ${log_path}/${SAMPLE1}_contamination.log

