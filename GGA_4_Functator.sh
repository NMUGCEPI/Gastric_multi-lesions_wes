#!/bin/sh
#PBS -V
#PBS -j oe
#PBS -l nodes=1:ppn=10,mem=10gb,walltime=120:00:00
#PBS -q batch

set -e

SAMPLE1=$1
SAMPLE2=$2
config_path=$3
source ${config_path}/config.sh


echo ${vcf_path}/${SAMPLE1}_${SAMPLE2}_GGA_Filter.vcf.gz
zcat ${vcf_path}/${SAMPLE1}_${SAMPLE2}_GGA_Filter.vcf.gz > ${vcf_path}/${SAMPLE1}_${SAMPLE2}_GGA_Filter.vcf

${gatk4} --java-options "-Xmx10g" Funcotator \
--variant ${vcf_path}/${SAMPLE1}_${SAMPLE2}_GGA_Filter.vcf \
--reference ${seq_ref}/GRCh37.fa \
--ref-version hg19 \
--data-sources-path ${funcotator_ref} \
--output ${maf_path}/${SAMPLE1}_${SAMPLE2}_GGA_Filter_funcotated.maf \
--output-file-format MAF \
--remove-filtered-variants true \
--transcript-selection-mode CANONICAL \
--force-b37-to-hg19-reference-contig-conversion \
2>  ${log_path}/${SAMPLE1}_${SAMPLE2}_GGA_Filter_funcotated.log


## https://gatk.broadinstitute.org/hc/en-us/community/posts/360060979451-Funcotator-b37-and-hg19-contig-compatibility-issue
## https://gatkforums.broadinstitute.org/gatk/discussion/11193/funcotator-information-and-tutorial

## grch37 and hg19 question
## WARN  FuncotatorEngine - WARNING: You are using B37 as a reference.  
## Funcotator will convert your variants to GRCh37, and this will be fine in the vast majority of cases.  
## There MAY be some errors (e.g. in the Y chromosome, but possibly in other places as well) due to changes between the two references.
