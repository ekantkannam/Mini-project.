# Mini-project.
Transcriptomic Analysis of Type 2 Diabetes using RNA-Seq

This repository contains the bioinformatics workflow for identifying differentially expressed genes (DEGs) in Type 2 Diabetes (T2D) using RNA sequencing (RNA-Seq) data. The study aims to explore transcriptomic alterations in T2D and uncover key molecular signatures that contribute to disease progression.

Project Overview

Objective
To analyze RNA-seq data and identify differentially expressed genes (DEGs) in T2D samples.
To perform pathway enrichment analysis to understand biological processes associated with DEGs.
To explore potential biomarkers and therapeutic targets for T2D.
Data Sources
RNA-seq datasets retrieved from NCBI Gene Expression Omnibus (GEO).
Blood and pancreatic islet transcriptomic profiles from healthy and T2D individuals.
Tools and Software Used
FastQC – Quality control of raw sequencing reads.
Trim Galore – Adapter trimming and read filtering.
HISAT2 – Read alignment to the reference genome.
FeatureCounts – Gene expression quantification.
DESeq2 – Differential gene expression analysis.
DAVID & KEGG – Functional enrichment and pathway analysis.
Python & R – Data processing, statistical analysis, and visualization.
Pipeline Workflow

Data Acquisition – Download RNA-seq datasets from GEO.
Preprocessing & Quality Control – Remove low-quality reads and adapters.
Read Alignment & Quantification – Map reads and compute gene expression levels.
Differential Expression Analysis – Identify upregulated/downregulated genes in T2D.
Functional Enrichment Analysis – Determine affected biological pathways and molecular functions.
Visualization & Interpretation – Generate PCA plots, volcano plots, heatmaps, and pathway enrichment charts.
Expected Outcomes

Identification of key DEGs associated with T2D.
Insights into dysregulated pathways linked to insulin resistance, glucose metabolism, and inflammation.
Potential biomarkers for early diagnosis and therapeutic interventions in T2D.
This repository provides all relevant scripts, workflows, and processed data files to ensure reproducibility and scalability for similar transcriptomic studies.
