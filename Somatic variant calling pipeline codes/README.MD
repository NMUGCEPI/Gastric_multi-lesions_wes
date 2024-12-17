# Gastric_MultipleRegion_WES

## Somatic variant calling pipeline codes for the Paper  
**"Evolving Genomic Characteristics of Intestinal Metaplasia and Gastric Cancer"**

This repository contains the core scripts used for preprocessing, mapping, and somatic variant calling in the study.

---

### **Step 1: Trim Adapters and Low-Quality Sequences**  
**Tool:** Trimmomatic (v0.36)

- **Code File**:  
  - `GATK4_pipeline_trim.sh`: Removes adapter sequences and filters out low-quality reads.

---

### **Step 2: Map Read Sequences to the Human Reference Genome**  
**Tools:** BWA-MEM (v0.7.15), GATK (v4.1.7)

#### **2.1 Map Reads and Mark Duplicates**  
- **Code File**:  
  - `GATK4_pipeline_bwamap.sh`: Maps reads to the reference genome and marks duplicate reads.

#### **2.2 Local Realignment and Base Quality Score Recalibration (BQSR)**  
- **Code File**:  
  - `GATK4_pipeline_bqsr.sh`: Performs local realignment around indels and recalibrates base quality scores.

---

### **Step 3: Call Somatic Variants**  
**Tool:** GATK (v4.1.7)

#### **3.1 Create Panel of Normals (PON)**  
- **Code Files**:  
  - `GATK4_pipeline_Step_Mutect2_1.1.sh`: Generate individual normal variant calls.  
  - `GATK4_pipeline_Step_Mutect2_1.2.sh`: Combine individual calls to create the PON.

#### **3.2 Call Raw Variants**  
- **Code File**:  
  - `GATK4_pipeline_Step_Mutect2_1.3.sh`: Perform raw variant calling using the PON.

#### **3.3 Calculate Contamination**  
- **Code Files**:  
  - `GATK4_pipeline_Step_Mutect2_2.1.sh`: Generate contamination metrics.  
  - `GATK4_pipeline_Step_Mutect2_2.2.sh`: Apply contamination estimates to improve variant calling.

#### **3.4 Filter Somatic Variants**  
- **Code Files**:  
  - `GATK4_pipeline_Step_Mutect2_2.3.sh`: Apply standard filters to raw variants.  
  - `GATK4_pipeline_Step_Mutect2_2.4.sh`: Further refine filtered somatic variants.

#### **3.5 Force-Calls Mode**  
Used to refine variant calls across samples.  
- **Code Files**:  
  - `GGA_1_CombineVcf.sh`: Combine VCF files from multiple samples.  
  - `GGA_2_Recall.sh`: Perform force-calls to recall variants across all samples.  
  - `GGA_3_Filter.sh`: Apply final filters to recalled variants.  
  - `GGA_4_Functator.sh`: Annotate variants with functional information.

#### **3.6 Annotation of Somatic SNVs and Indels**  
- **Code File**:  
  - `GGA_4_Functator.sh`: Annotate somatic single nucleotide variants (SNVs) and insertions/deletions (indels) using GATK Funcotator.

---

## **Usage Instructions**
1. Clone this repository.  
2. Ensure all dependencies are installed (e.g., Trimmomatic, BWA-MEM, GATK).  
3. Run each script in the specified order to replicate the preprocessing and somatic variant calling pipeline.

---

## **Contact**
For inquiries or support, please contact [xuxianfeng@njmu.edu.cn].

---
