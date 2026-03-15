# The Developmental and Transcriptomic Identities of GLP-1R Cells

## File Description

- `README.md`: Short project overview and guide to the files in this repository
- `Project_Code_Final.Rmd`: Seurat workflow for loading the DVC dataset, checking QC metrics, subsetting Glp1r/Gfral-positive cells, clustering them, and identifying marker genes
- `project_clean_GSEA_clean.Rmd`: Cleaned rmd analysis that prepares receptor positive cell objects and runs the project GSEA workflow across receptor, diet, neuron, and glial contrasts
- `cellcomp_feeding_analysis_final.R`: R script for converting the source `.h5ad` data to Seurat, annotating receptor positive cells, creating neuron only subsets, and generating feeding state composition plots
- `sig_fgsea.py`: Python helper script that scans existing FGSEA result files and exports statistically significant pathways into summary tables
- `requirements_glp1r_pilot.txt`: Python package list for the pilot GLP1R analysis notebooks
- `nat_analysis_final.ipynb`: Notebook for GLP1R+ cell analysis, including dataset curation, downsampling, and feature modeling with LASSO and XGBoost
- `FINAL_module1.ipynb`: Notebook for pseudobulk-style developmental expression profiling of incretin-axis genes across cell types and time points
- `FINAL_module2.ipynb`: Notebook for developmental trajectory-style analysis using CellxGene Census data, including UMAP, PAGA, and CytoTRACE2


## Data
Trajectory Data:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186069

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186069

(from this paper: https://elifesciences.org/articles/106217)

DVC data: https://singlecell.broadinstitute.org/single_cell/study/SCP2773/a-unified-rodent-atlas-reveals-the-cellular-complexity-and-evolutionary-divergence-of-the-dorsal-vagal-complex?genes=Lmx1b&tab=dotplot

(from this paper: https://www.nature.com/articles/s41586-024-07069-w)
