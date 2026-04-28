# Glycosylation-Associated Prognostic Modeling Across Female Cancers

## Overview
This repository presents a multi-omics framework integrating bulk RNA-seq and single-cell RNA-seq data to identify glycosylation-associated prognostic signatures across female malignancies.

A robust five-gene signature (**TFF3, DAPL1, CP, SLC16A3, MUC20**) was developed using a multi-algorithm machine learning strategy. Among these, **SLC16A3** was prioritised for experimental validation and structural analysis through molecular docking and molecular dynamics simulations.

---

## Repository Structure

| Category | File | Description |
|----------|------|------------|
|  Survival and Prognostic Modeling | **survival_analysis.R** | End-to-end survival analysis pipeline |
|  | **Cox.R** | Univariate and multivariate Cox regression |
|  | **KM_curves.R** | Kaplan–Meier survival curve generation |
|  | **pan_cancer_survival_all_models.R** | Pan-cancer survival modeling |
|  Machine Learning and Feature Selection | **feature_selection.py** | Feature selection framework |
|  | **deep_learning_models.py** | Deep learning-based models |
|  | **RF_GBM_Gene_Signatures.r** | Random Forest and GBM-based gene prioritisation |
|  Functional and Pathway Analysis | **GO_KEGG_analysis.R** | GO and KEGG enrichment analysis |
|  | **ssGSEA.R** | Single-sample GSEA scoring |
|  | **AUCcell.R** | Gene set activity scoring (single-cell) |
|  Single-Cell RNA-seq Analysis | **scRNA-seq analysis.R** | Clustering, annotation, and marker identification |
|  | **CellChat.R** | Cell–cell communication analysis |
|  Visualization | **box_plots.R** | Gene expression boxplots |
---

## Requirements

### R Packages
- survival  
- survminer  
- Seurat  
- GSVA  
- clusterProfiler  
- CellChat  

### Python Libraries
- scikit-learn  
- lifelines  
- pandas  
- numpy  
- tensorflow / pytorch  

---

## Key Highlights
- Multi-omics integration (bulk + single-cell)  
- Cross-cancer robust prognostic signature  
- Multi-algorithm validation (RF, GBM, DL)  
- Experimental validation of **SLC16A3**  
- Structure-based drug candidate identification  

---

## Citation
If you use this repository, please cite the associated manuscript (to be updated).
