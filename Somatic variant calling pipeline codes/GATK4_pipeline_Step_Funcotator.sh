#!/bin/sh
#PBS -V
#PBS -j oe
#PBS -l nodes=1:ppn=10,mem=10gb,walltime=120:00:00
#PBS -q batch
#PBS -d /public/home/xxf/20200513_colr_GATKNew/qsub_scripts

source /public/home/xxf/20200513_colr_GATKNew/config/config.sh
set -e

SAMPLE1=$1
SAMPLE2=$2

echo ${SAMPLE1}_${SAMPLE2}

${gatk4} --java-options "-Xmx10g" Funcotator \
--variant ${vcf_path}/${SAMPLE1}_${SAMPLE2}_filtered.vcf \
--reference ${seq_ref}/GRCh37.fa \
--ref-version hg19 \
--data-sources-path ${funcotator_ref} \
--output ${vcf_path}/${SAMPLE1}_${SAMPLE2}_funcotated.maf \
--output-file-format MAF \
--remove-filtered-variants true \
--transcript-selection-mode CANONICAL \
--force-b37-to-hg19-reference-contig-conversion \
2>  ${log_path}/${SAMPLE1}_${SAMPLE2}_Funcotator.log


## https://gatk.broadinstitute.org/hc/en-us/community/posts/360060979451-Funcotator-b37-and-hg19-contig-compatibility-issue
## https://gatkforums.broadinstitute.org/gatk/discussion/11193/funcotator-information-and-tutorial

## grch37 and hg19 question
## WARN  FuncotatorEngine - WARNING: You are using B37 as a reference.  
## Funcotator will convert your variants to GRCh37, and this will be fine in the vast majority of cases.  
## There MAY be some errors (e.g. in the Y chromosome, but possibly in other places as well) due to changes between the two references.
