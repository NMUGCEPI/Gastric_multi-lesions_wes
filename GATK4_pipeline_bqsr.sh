#!/bin/sh
#PBS -V
#PBS -j oe
#PBS -l nodes=1:ppn=1,mem=10gb,walltime=120:00:00
#PBS -q batch

sample=$1
config_path=$2
source ${config_path}/config.sh

fasta = "${params.genome_37}/GRCh37.fa"
fai = "${params.genome_37}/GRCh37.fa.fai"
dbSNP_vcf="${params.variants_ref_37}/dbsnp_138.b37.vcf.gz"
known_indels_sites_VCFs_1="${params.variants_ref_37}/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
known_indels_sites_VCFs_2="${params.variants_ref_37}/1000G_phase1.indels.b37.vcf.gz"
heterochromatic="${params.genome_37}/GRCh37_otherChr_region.bed"

## baseRecalibrator
${params.gatk4} --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
-XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
-Xloggc:${log_path}/${sample}_gc_log.BaseRecalibrator.log -Xms5g \
-Djava.io.tmpdir=${tmp_path}/ " \
BaseRecalibrator \
-R ${ref_fasta} \
-I ${sample}.aligned.duplicate_marked.sorted.bam \
--use-original-qualities \
-O ${sample}.recal_data.table \
--known-sites ${params.dbSNP_vcf} \
--known-sites ${params.known_indels_sites_VCFs_1} \
--known-sites ${params.known_indels_sites_VCFs_2} \
-L ${target_region} \
2> ${log_path}/${sample}.BaseRecalibrator.log

## applyBQSR
${params.gatk4} --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
-XX:+PrintGCDetails -Xloggc:${log_path}/${sample}_gc_log.ApplyBQSR.log \
-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level=${compression_level} \
-Djava.io.tmpdir=${PWD}/tmp/ \
-Xms5g" \
ApplyBQSR \
-R ${ref_fasta} \
-I ${sample}.aligned.duplicate_marked.sorted.bam \
-O ${sample}-recal.bam \
--create-output-bam-index true \
-L ${target_region} \
-bqsr ${sample}.recal_data.table \
--add-output-sam-program-record \
--create-output-bam-md5 \
--use-original-qualities \
--emit-original-quals true \
2> ${log_path}/${sample}.ApplyBQSR.log

${params.samtools} index ${sample}-recal.bam 
cp ${sample}-recal.bam.bai ${sample}-recal.bai
md5sum ${sample}-recal.bam > ${sample}-recal.bam.md5
