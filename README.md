# SpatioTemporal Malaria Liver-stage Atlas

This code is associated with the following manuscript:
https://www.nature.com/articles/s41586-022-05406-5
and allows you to go from raw data into a usable Seurat/Scanpy structure for further analysis.

__Files:__

* P01_raw_data_preprocessing.ipynb: Jupyter notebook; takes raw data and formats into [Scanpy](https://github.com/theislab/scanpy) Anndata
* P02_Scanpy_processing.ipynb: Jupyter notebook; further processing and filtering of the Anndata structure
* R01_Scanpy2Seurat.R: R; transfers Anndata to a [Seurat](https://github.com/satijalab/seurat) object, then runnning a basic seurat normalization, scaling and dimensionality reduction
* /input_data: folder containing various metadata tables needed for the processing of the raw data

__The output of the files was the basis to all scRNAseq analysis done in the manuscript. For further details see Methods section in manuscript.__
