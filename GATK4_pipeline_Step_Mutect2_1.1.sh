#!/bin/sh
#PBS -V
#PBS -j oe
#PBS -l nodes=1:ppn=1,mem=10gb,walltime=120:00:00
#PBS -q batch

sample=$1
config_path=$2
source ${config_path}/config.sh

set -e

#1. Create a sites-only PoN with CreateSomaticPanelOfNormals
${gatk4} --java-options "-Xmx20g" Mutect2 \
-R ${seq_ref}/GRCh37.fa \
-I ${bam_path}/${sample}-recal.bam \
--max-mnp-distance 0 \
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
-L ${target_bed} \
--native-pair-hmm-threads 16 \
-O ${vcf_path}/${sample}_pon.vcf.gz \
2>  ${log_path}/${sample}_pon.log

