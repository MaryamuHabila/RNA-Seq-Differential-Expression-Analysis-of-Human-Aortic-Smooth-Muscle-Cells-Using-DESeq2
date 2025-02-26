# RNA-Seq-Differential-Expression-Analysis-of-Human-Aortic-Smooth-Muscle-Cells-Using-DESeq2

This project aims to analyze RNA sequencing data from human aortic smooth muscle cells subjected to a siRNA-mediated knockdown of a target gene. The goal is to identify differentially expressed genes between Control and Treatment groups using the Salmon quantification tool for transcript-level expression estimation, followed by DESeq2 for differential expression analysis.

The analysis involves:
Preprocessing and organizing metadata for samples
Importing transcript quantifications from Salmon
Mapping transcript-level data to gene-level data using GTF annotations
Performing differential expression analysis using DESeq2
Generating a results table with gene annotations. This pipeline allows for reproducible transcriptomic analysis and can be adapted for other RNA-seq datasets.

1. Installation of Required Packages
Several R packages are installed to facilitate transcript quantification, data handling, and statistical analysis:

tximport: To import transcript abundance estimates from Salmon
readr, data.table, tidyverse: For data handling
DESeq2: To perform differential gene expression analysis
apeglm: For shrinkage of log fold changes
rtracklayer: To parse GTF gene annotations

2. Sample Metadata Preparation
A metadata table (samples) is created to define sample groupings: Sample IDs uniquely identify each RNA-seq sample. Condition Labels specify whether a sample belongs to the "Control" or "Treatment" group. This metadata table is crucial for downstream analysis, ensuring correct sample assignments in statistical modeling.

3. Linking Sample IDs to Salmon Quantification Data
File paths for Salmon .quant.sf files are constructed using sample names. Each sample's file path is stored in a named vector (files). The existence of these files is verified to prevent errors in later steps.

4. Annotation Data: Mapping Transcripts to Genes.
To convert transcript-level quantification to gene-level: A GTF annotation file (hg38.refGene.gtf) is imported using rtracklayer. A tx2gene table is created, linking transcript IDs to gene IDs. Metadata such as genomic coordinates, gene names, and exon information are extracted. A gene metadata table is derived by removing duplicate gene entries, ensuring gene-level uniqueness. This mapping is essential for aggregating transcript expression into gene expression.

5. Importing Transcript-Level Data with tximport
tximport is used to aggregate transcript-level quantification from Salmon into gene-level expression data. The tx2gene table is applied to correctly map transcripts to genes. Transcript version numbers are ignored (ignoreTxVersion = TRUE) to ensure compatibility across different database versions. This step prepares the dataset for DESeq2-based differential expression analysis.

6. Differential Expression Analysis with DESeq2: A subset of three Control and three Treatment samples is selected for comparison.
The DESeq2 workflow consists of:
Constructing a DESeq2 dataset (DESeqDataSetFromTximport) using the gene-level counts and sample metadata.
Running DESeq() to perform normalization, dispersion estimation, and statistical testing.
Extracting differential expression results, focusing on log fold change, adjusted p-values, and raw counts.
Merging results with gene metadata to include gene names and annotations.
Exporting the results to a CSV file for further downstream analysis or visualization.

7. Output and Final Data Export
The final differential expression table is saved as "DESeq2_results_subset4.csv", containing:
Gene IDs
Log2 fold changes
Statistical significance (adjusted p-values)
Gene annotations (from the GTF file)
This output serves as the basis for further functional enrichment analysis and biological interpretation.
