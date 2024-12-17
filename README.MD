# Gastric_MultipleRegion_WES Project

## Overview
This repository contains the codebase used in the study  
**"Evolving Genomic Characteristics of Intestinal Metaplasia and Gastric Cancer"**,  
covering both data preprocessing pipelines and figure generation scripts.

The project focuses on exploring the genomic characteristics of intestinal metaplasia (IM), intestinal gastric cancer (IGC), and diffuse gastric cancer (DGC). The code is divided into two main sections:

1. **Data Processing Pipelines**: Scripts for preprocessing raw sequencing data, mapping reads, and identifying somatic variants.  
2. **Figure Generation**: Scripts used to produce the figures presented in the study.

---

## **1. Data Processing Pipelines**

### **Purpose**  
These scripts are used to process whole-exome sequencing (WES) data, from trimming raw reads to annotating somatic mutations.

### **Steps and Scripts**
#### **Step 1: Trim Adapters and Low-Quality Reads**  
- **Script**: `GATK4_pipeline_trim.sh`  
  - Tool: Trimmomatic (v0.36)  
  - Function: Removes adapter sequences and filters low-quality reads.

#### **Step 2: Map Reads to Reference Genome**  
- **Scripts**:  
  - `GATK4_pipeline_bwamap.sh`: Map reads and mark duplicates using BWA-MEM (v0.7.15).  
  - `GATK4_pipeline_bqsr.sh`: Perform local realignment and base quality score recalibration (BQSR) using GATK (v4.1.7).

#### **Step 3: Somatic Variant Calling**  
- **Scripts**:  
  - `GATK4_pipeline_Step_Mutect2_*.sh`: Perform somatic variant calling using GATK Mutect2, including:  
    - Creating a Panel of Normals (PON).  
    - Calling raw variants.  
    - Contamination estimation.  
    - Filtering somatic variants.  
  - `GGA_*.sh`: Refine and annotate variants through force-calling and Funcotator.

---

## **2. Figure Generation**

### **Purpose**  
These scripts are used to analyze processed data and generate publication-quality figures for the study.

### **Figures and Scripts**
#### **Figure 1**  
Focus: Mutation burden, significantly mutated genes, and clonal evolution patterns in IM, IGC, and DGC.  
- **Scripts**:  
  - `Figure1_gcburden.R`  
  - `Figure1_imburden.R`  
  - `Figure1_waterfall.R`  
  - `Figure2_evolutionarymode.R`

#### **Figure 2**  
Focus: Comparative mutation burden and CNVs in IGC and DGC.  
- **Scripts**:  
  - `Figure2_mutationburden.R`  
  - `Figure2_cnvburden.R`

#### **Figure 3**  
Focus: Mutation burden, tumor cell fractions, and CNVs in IM.  
- **Scripts**:  
  - `Figure3_mutationburden.R`  
  - `Figure3_cnvburden.R`

#### **Figure 4**  
Focus: Mutation waterfall plots, mutation rate comparisons, and clonal evolution patterns in IM, IGC, and DGC.  
- **Scripts**:  
  - `Figure4_waterfall.R`  
  - `Figure4_mutrate.R`  
  - `Figure4_modes.R`

#### **Figure 5**  
Focus: Mutation waterfall plot.  
- **Script**:  
  - `Figure5.R`

#### **Figure 6**  
Focus: Comparison of clonal evolution timing.  
- **Script**:  
  - `Figure6.R`

#### **Figure 7**  
Focus: MUC6 mutation distribution, epithelial proportions, cell-cell communication, and differential expression.  
- **Scripts**:  
  - `Figure7_mutation.R`  
  - `Figure7_compareratio.R`  
  - `Figure7_scissor.R`  
  - `Figure7_cellphonedb.R`  
  - `Figure7_diffexp.R`  
  - `Figure7_pitratio.R`  
  - `Figure7_bulk.R`

---

## **Usage Instructions**

### **Data Processing**
1. Clone the repository.
2. Install required tools (e.g., Trimmomatic, BWA-MEM, GATK).
3. Execute scripts in the specified order under the "Data Processing Pipelines" section.

### **Figure Generation**
1. Ensure processed data is available.
2. Install R and required R packages.
3. Run the relevant R scripts to generate each figure.

---

## **Contact**
For questions or issues regarding the code, please contact [Your Name/Email/Institution].

---
