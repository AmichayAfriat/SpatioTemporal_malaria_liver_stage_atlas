# SpatioTemporal_malaria_liver_stage_atlas

This code is associated with the following manuscript:
(biorxiv-link)
and allows you to go from raw data into a usable Seurat/Scanpy structure for further analysis.

__NOTE! This code is based on Raw mRNA sequencing data__

Files have been deposited in the GenBank GEO database under accession code GSEXXXXXXX
and should be downloaded prior to running this script
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXXXXXX

__Files__

P01_raw_data_preprocessing.ipynb - Jupyter notebook; takes raw data and formats into Scanpy Anndata
P02_Scanpy_processing.ipynb - Jupyter notebook; further processing and filtering of the Anndata structure
R01_Scanpy2Seurat.R - R; transfers Anndata to a Seurat object, then runnning a basic seurat normalization, scaling and dimensionality reduction

/input_data - folder containing various metadata tables needed for the processing of the raw data

__The output of the files was the basis to all scRNAseq analysis done in the manuscript. For further details see Methods section in manuscript.__
