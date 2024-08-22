#!/bin/sh
#PBS -V
#PBS -j oe
#PBS -l nodes=1:node8:ppn=16,mem=10gb,walltime=120:00:00
#PBS -q batch
#PBS -d /public/home/xxf/20200513_colr_GATKNew/Qsub_log

config_path=$1
source ${config_path}/config.sh

set -e

Normal_vcfs=`ls ${vcf_path} | grep pon | awk -F'_' '{print $1}' | sort -u | awk '{print "-V",vcf_path"/"$0"_pon.vcf.gz"}' vcf_path=${vcf_path} `

${gatk4} --java-options "-Xmx10g" GenomicsDBImport \
--genomicsdb-workspace-path ${vcf_path}/pon_db \
-R ${seq_ref}/GRCh37.fa \
-L ${target_bed} \
--merge-input-intervals true \
${Normal_vcfs} \
--reader-threads 16 \
2>  ${log_path}/PoN_GenomicsDBImport.log


${gatk4} --java-options "-Xmx10g" CreateSomaticPanelOfNormals \
--germline-resource ${variants_ref}/af-only-gnomad.raw.sites.b37.vcf.gz \
-R ${seq_ref}/GRCh37.fa \
-V gendb://${vcf_path}/pon_db \
-O ${vcf_path}/PON.vcf.gz \
2>  ${log_path}/PoN_CreateSomaticPanelOfNormals.log



# https://gatk.broadinstitute.org/hc/en-us/articles/360042477052-GenomicsDBImport
# --merge-input-intervals false
# Boolean flag to import all data in between intervals. 
# Improves performance using large lists of intervals, as in exome sequencing, especially if GVCF data only exists for specified intervals.