#!/bin/sh
#PBS -V
#PBS -j oe
#PBS -l nodes=1:ppn=1,mem=10gb,walltime=120:00:00
#PBS -q batch

sample=$1
config_path=$2
source ${config_path}/config.sh
compression_level = 5

## samToFastqAndBwaMem
${params.bwa} mem -K 100000000 -v 3 -t ${task.cpus} -R "@RG\\tID:${sample}\\tPL:illumina\\tPU:${sample}_LCB_WGS\\tSM:${sample}" \
${ref_fasta} ${params.inputdir}/${sample}_R1_trim.fq.gz ${params.inputdir}/${sample}_R2_trim.fq.gz \
| ${params.samtools} view -S -b - \
1> ${sample}.unmerged.bam \
2> ${log_path}/${sample}.unmerged.log

## bamTofastqAndBwaMem
java -Djava.io.tmpdir=${tmp_path} -Dsamjdk.compression_level=${compression_level} -Xms5G \
-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
-XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
-Xloggc:${log_path}/${sample}_gc_log.RevertSam.log \
-jar ${params.picard} \
RevertSam \
INPUT=${params.inputdir}/${sample}.raw.bam \
OUTPUT=${sample}.unmapped_bam.bam \
2> ${log_path}/${sample}-RevertSam.log
java -Djava.io.tmpdir=${tmp_path} -Dsamjdk.compression_level=${compression_level} -Xms5G \
-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
-XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
-Xloggc:${log_path}/${sample}_gc_log.SamToFastq.log \
-jar ${params.picard}  \
SamToFastq \
INPUT=${sample}.unmapped_bam.bam \
FASTQ=${sample}.ALIGNED.fq \
INTERLEAVE=true \
NON_PF=true \
VALIDATION_STRINGENCY=SILENT \
2> ${log_path}/${sample}_SamToFastq.log
${params.bwa} mem -K 100000000 -p -v 3 -t ${task.cpus} -R "@RG\\tID:${sample}\\tPL:illumina\\tPU:${sample}_LCB_WGS\\tSM:${sample}" \
-Y ${ref_fasta} ${sample}.ALIGNED.fq | \
${params.samtools} view -S -b - \
1> ${sample}.aligned.unsorted.bam \
2> ${log_path}/${sample}.aligned.unsorted.log

## mergeBamAlignment
${params.gatk4} --java-options " -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
-XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
-Xloggc:${log_path}/${sample}_gc_log.MergeBamAlignment.log \
-Dsamjdk.compression_level=${compression_level} -Xms5G -Xmx10G -Djava.io.tmpdir=${tmp_path}" \
MergeBamAlignment \
--VALIDATION_STRINGENCY SILENT \
--EXPECTED_ORIENTATIONS FR \
--ATTRIBUTES_TO_RETAIN X0 \
--ALIGNED_BAM ${sample}.unmerged.bam \
--UNMAPPED_BAM ${sample}.unmapped_bam.bam \
--OUTPUT ${sample}.aligned.unsorted.bam \
--REFERENCE_SEQUENCE ${ref_fasta} \
--PAIRED_RUN true \
--SORT_ORDER "unsorted" \
--IS_BISULFITE_SEQUENCE false \
--ALIGNED_READS_ONLY false \
--CLIP_ADAPTERS false \
--MAX_RECORDS_IN_RAM 2000000 \
--ADD_MATE_CIGAR true \
--MAX_INSERTIONS_OR_DELETIONS -1 \
--PRIMARY_ALIGNMENT_STRATEGY MostDistant \
--PROGRAM_RECORD_ID "bwamem" \
--PROGRAM_GROUP_VERSION "0.7.17-r1188" \
--PROGRAM_GROUP_COMMAND_LINE "bwa mem -K 100000000 -p -v 3 -t ${task.cpus} -R "@RG\\tID:${sample}\\tPL:illumina\\tPU:${sample}_LCB_WGS\\tSM:${sample}" -Y ${ref_fasta}" \
--PROGRAM_GROUP_NAME "bwamem" \
--UNMAPPED_READ_STRATEGY COPY_TO_TAG \
--ALIGNER_PROPER_PAIR_FLAGS true \
--UNMAP_CONTAMINANT_READS true \
2> ${log_path}/${sample}_MergeBamAlignment.log

## markDuplicates
${params.gatk4} --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
-XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
-Xloggc:${log_path}/${sample}_gc_log.MarkDuplicates.log \
-Dsamjdk.compression_level=${compression_level} -Xms5G -Djava.io.tmpdir=${tmp_path}" \
MarkDuplicates \
--INPUT ${sample}.aligned.unsorted.bam \
--OUTPUT ${sample}.aligned.unsorted.duplicates_marked.bam \
--METRICS_FILE ${bam_path}/${sample}.duplicate_metrics \
--VALIDATION_STRINGENCY SILENT \
--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
--ASSUME_SORT_ORDER "queryname" \
--CREATE_MD5_FILE true \
2> ${log_path}/${sample}_MarkDuplicates.log

## sortAndFixTags
${params.gatk4} --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
-XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
-Xloggc:${log_path}/${sample}_gc_log.SortSam.log \
-Dsamjdk.compression_level=${compression_level} -Xms5G -Xmx10G -Djava.io.tmpdir=${tmp_path}" \
SortSam \
--INPUT ${sample}.aligned.unsorted.duplicates_marked.bam \
--OUTPUT ${sample}.aligned.unsorted.duplicates_marked.tmp.bam \
--SORT_ORDER "coordinate" \
--CREATE_INDEX false \
--CREATE_MD5_FILE false 

${params.gatk4} --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
-XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
-Xloggc:${log_path}/${sample}_gc_log.SetNmMdAndUqTags.log \
-Dsamjdk.compression_level=${compression_level} -Xms5G -Xmx10G -Djava.io.tmpdir=${tmp_path}" \
SetNmMdAndUqTags \
--INPUT ${sample}.aligned.unsorted.duplicates_marked.tmp.bam \
--OUTPUT ${sample}.aligned.duplicate_marked.sorted.bam \
--CREATE_INDEX true \
--CREATE_MD5_FILE true \
--REFERENCE_SEQUENCE ${ref_fasta} \
2> ${log_path}/${sample}_sortAndFixTags.log
