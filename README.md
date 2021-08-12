# SpatioTemporal Malaria Liver-stage Atlas

This code is associated with the following manuscript:
(biorxiv-link)
and allows you to go from raw data into a usable Seurat/Scanpy structure for further analysis.

__NOTE! This code is based on Raw mRNA sequencing data__

Files have been deposited in the GenBank GEO database under accession code GSE181725
and should be downloaded prior to running this script
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181725



__Files:__

* P01_raw_data_preprocessing.ipynb: Jupyter notebook; takes raw data and formats into [Scanpy](https://github.com/theislab/scanpy) Anndata
* P02_Scanpy_processing.ipynb: Jupyter notebook; further processing and filtering of the Anndata structure
* R01_Scanpy2Seurat.R: R; transfers Anndata to a [Seurat](https://github.com/satijalab/seurat) object, then runnning a basic seurat normalization, scaling and dimensionality reduction
* /input_data: folder containing various metadata tables needed for the processing of the raw data

__The output of the files was the basis to all scRNAseq analysis done in the manuscript. For further details see Methods section in manuscript.__
