#!/bin/sh
#PBS -V
#PBS -j oe
#PBS -l nodes=1:ppn=1,mem=1gb,walltime=120:00:00
#PBS -q batch

set -e

SAMPLE1=$1
SAMPLE2=$2
config_path=$3
source ${config_path}/config.sh

echo ${SAMPLE1}_${SAMPLE2}

${gatk4} --java-options "-Xmx1g" FilterMutectCalls \
-V ${vcf_path}/${SAMPLE1}_${SAMPLE2}_unfiltered.vcf.gz \
-R ${seq_ref}/GRCh37.fa \
--contamination-table ${vcf_path}/${SAMPLE1}_contamination.table \
--tumor-segmentation ${vcf_path}/${SAMPLE1}_segments.table \
--ob-priors ${vcf_path}/${SAMPLE1}_${SAMPLE2}_read-orientation-model.tar.gz \
-O ${vcf_path}/${SAMPLE1}_${SAMPLE2}_filtered.vcf \
2>  ${log_path}/${SAMPLE1}_FilterMutectCalls.log