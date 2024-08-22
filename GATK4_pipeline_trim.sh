#!/bin/sh
#PBS -V
#PBS -j oe
#PBS -l nodes=1:ppn=1,mem=10gb,walltime=120:00:00
#PBS -q batch

sample=$1
config_path=$2
source ${config_path}/config.sh

java -Xmx32g -jar ${trimmomatic} \
PE -phred33 -threads ${trim_threads} \
${raw_fastq_path}/${sample}.R1.clean.fastq.gz \
${raw_fastq_path}/${sample}.R2.clean.fastq.gz \
${trim_fastq_path}/${sample}_R1_trim.fq.gz ${trim_fastq_path}/${sample}_R1_trim_up.fq.gz \
${trim_fastq_path}/${sample}_R2_trim.fq.gz ${trim_fastq_path}/${sample}_R2_trim_up.fq.gz \
ILLUMINACLIP:${trim_adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 \
2> ${log_path}/${sample}-trim.log 



