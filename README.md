# DV-home-assigment
This repository contains the work for **Home Assignment 2** of the course **Data visualization (2025-2026)** in the Bachelor's Degree in Bioinformatics.

---
## Authors
- Alicia Mañas
- Paula Artiz
- Lidia Sánchez
---
## Repository contents

The dataset used in this project consists of gene expression counts from four count matrices:
  - 'rep1.csv'
  - 'rep2.csv'
  - 'rep3.csv'
  - 'rep4.csv'
 
Each file contains:
  - Rows : genes
  - Columns : samples
 
The data was merged into a single expression matrix and combined with metadata extracted from sample names (wine type, time point, and replicate). 

---


## Key findings

- Strong batch effects were present in raw data due to differences in library size
- Normalisation was essential to recover biological signal
- After correction, PCA revealed clearer separation of biological conditions
- t-SNE was sensitive to parameter choice but consistent in global structure
- Replicates influenced raw structure more than biological conditions

---

## Tools used

- R
- ggplot2
- factoextra
- ggfortify
- Rtsne
- dplyr
- gtsummary


