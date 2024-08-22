# gastric_mutipleregion_wes

## Main code for paper "Evolving genomic characteristics of intestinal metaplasia and gastric cancer"

### Trim adapter and low-quality sequence (Trimmomatic (v0.36))

GATK4_pipeline_trim.sh

### Map read sequences to the human reference genome (BWA-MEM (v0.7.15), GATK (v4.1.7))

#### Map reads and mark duplicates 

GATK4_pipeline_bwamap.sh

#### Local realignment and base quality score recalibration

GATK4_pipeline_bqsr.sh

### Call somatic variants (GATK (v4.1.7))

#### Create PON

GATK4_pipeline_Step_Mutect2_1.1.sh 

GATK4_pipeline_Step_Mutect2_1.2.sh

#### Call raw variants

GATK4_pipeline_Step_Mutect2_1.3.sh

#### Calculate contamination

GATK4_pipeline_Step_Mutect2_2.1.sh

GATK4_pipeline_Step_Mutect2_2.2.sh

#### Filter somatic variants

GATK4_pipeline_Step_Mutect2_2.3.sh

GATK4_pipeline_Step_Mutect2_2.4.sh

#### Force-calls mode

GGA_1_CombineVcf.sh

GGA_2_Recall.sh

GGA_3_Filter.sh

GGA_4_Functator.sh

#### Annotation somatic SNVs and indels

GGA_4_Functator.sh
