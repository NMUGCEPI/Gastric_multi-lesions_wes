#!/bin/sh
#PBS -V
#PBS -j oe
#PBS -l nodes=1:ppn=1,mem=10gb,walltime=120:00:00
#PBS -q batch


SAMPLE1=$1
SAMPLE2=$2
config_path=$3

source ${config_path}/config.sh

set -e

echo ${SAMPLE1}_${SAMPLE2}

${gatk4} --java-options "-Xmx10g" Mutect2 \
-R ${seq_ref}/GRCh37.fa \
-I ${bam_path}/${SAMPLE1}-recal.bam \
-I ${bam_path}/${SAMPLE2}-recal.bam \
-tumor ${SAMPLE1} \
-normal ${SAMPLE2} \
--panel-of-normals ${vcf_path}/PON.vcf.gz \
--germline-resource ${variants_ref}/af-only-gnomad.raw.sites.b37.vcf.gz \
--af-of-alleles-not-in-resource 0.0000025 \
-L ${target_bed} \
--native-pair-hmm-threads 16 \
--f1r2-tar-gz ${vcf_path}/${SAMPLE1}_${SAMPLE2}_f1r2.tar.gz \
-O ${vcf_path}/${SAMPLE1}_${SAMPLE2}_unfiltered.vcf.gz \
-bamout ${vcf_path}/${SAMPLE1}_${SAMPLE2}_mutect2.bam \
2>  ${log_path}/${SAMPLE1}_${SAMPLE2}_mutect2.log



${gatk4} LearnReadOrientationModel \
-I ${vcf_path}/${SAMPLE1}_${SAMPLE2}_f1r2.tar.gz \
-O ${vcf_path}/${SAMPLE1}_${SAMPLE2}_read-orientation-model.tar.gz \
2>  ${log_path}/${SAMPLE1}_${SAMPLE2}_LearnReadOrientationModel.log
