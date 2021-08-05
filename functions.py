import scanpy as sc
import numpy as np
import pandas as pd
import os
import re
from anndata import AnnData
from scipy.spatial.distance import correlation, cosine, euclidean
import matplotlib.pyplot as plt

os.environ['R_HOME'] = r'C:\Program Files\R\R-4.0.2'
os.environ['R_USER'] = r'D:\anaconda3\envs\malaria\Lib\site-packages\rpy2'
from rpy2.robjects import pandas2ri
import rpy2.robjects as robjects
pandas2ri.activate()
from rpy2.robjects.packages import STAP

def grouped_obs_mean(adata, groupby):
    """ Given an AnnData object with obs annotations in column 'groupby', groups
    cells and calculates mean expression values (raw counts) which are written
    into a dataframe with groupby categories as columns and genes as index.
    
    Source
    ------
    adapted from https://github.com/theislab/scanpy/issues/181, ivirshup
    """    
    # group obs df
    grouped = adata.obs.groupby(groupby)
    
    # prepare pandas dataframe for mean gene values in groups
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )
    
    # for each group, get the data of member cells and average it,
    # then write to dataframe columns
    for group, idx in grouped.indices.items():
        X = adata[idx].X
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))
        
    return out


def stretch_rotate(
        x_original,
        y_original,
        start_point,
        end_point):
    """ Takes a single point specified by x_original, y_original and transforms
    it onto a new linear axis between start_point and end_point."""
    start_point = np.array(start_point)
    end_point = np.array(end_point)
    # first, calculate the length of the distance between start and end
    length = euclidean(start_point, end_point)
    # we need the angle between start and end
    # because of cosine degenracy we need to check the angle with (1,0) or (-1,0)
    # depending on whether the x val of start-end is negative or positive
    cos = 1 - cosine(end_point-start_point, np.array([1, 0]))
    angle = np.arccos(np.clip(cos, -1, 1))
    
    # rotate the original point
    x_rot = np.cos(angle) * x_original - np.sin(angle) * y_original
    y_rot = np.sin(angle) * x_original + np.cos(angle) * y_original 
    
    # next, the vector still needs to be scaled
    k_x = np.max([1, (end_point[0] - start_point[0])/np.cos(angle)])
    k_y = (end_point[1] - start_point[1])/np.sin(angle)
    
    x_scaled = x_rot*k_x
    y_scaled = y_rot*k_y
        
    # last, we shift the new point to the correct position related to the start point
    x_shifted = x_scaled + start_point[0]
    y_shifted = y_scaled + start_point[1]   
    
    return x_shifted, y_shifted


def rank_genes_groups_to_df(adata, group, pval_cutoff : float =None, logfc_cutoff=None): 
    d = pd.DataFrame() 
    for k in ['scores', 'names', 'logfoldchanges', 'pvals', 'pvals_adj']: 
        d[k] = adata.uns["rank_genes_groups"][k][group] 
    if pval_cutoff is not None: 
        d = d[d["pvals_adj"] < pval_cutoff] 
    if logfc_cutoff is not None: 
        d = d[d["logfoldchanges"].abs() > logfc_cutoff] 
    return d

from textwrap import wrap
# plot function for profiles
def plot_30_profiles(
                df_mean, 
                df_std,
                gene_symbols,
                xlabel=None,
                xpos1=None):
    # make figure grid 5x4
    fig, axes = plt.subplots(5, 6, figsize=(20, 13), sharex=True, sharey=False)
    # get x values
    ticks = [0, 1]
    colnames = ['0', '1']
    if xpos1==None:
        xpos1 = np.linspace(0, 1, len(df_mean.columns))

    # go through gene symbol list supplied and plot mean +- sem onto axes
    for i, ax in enumerate(axes.flat):
        try:
            # data1
            ax.plot(xpos1, df_mean.loc[gene_symbols[i]].values, '-', color='black')
            ax.fill_between(xpos1,
                            df_mean.loc[gene_symbols[i]].values - df_std.loc[gene_symbols[i]].values,
                            df_mean.loc[gene_symbols[i]].values + df_std.loc[gene_symbols[i]].values,
                            alpha=0.3, color='palegreen')

            ax.set_xticks(ticks)
            ax.set_xticklabels(colnames)
            ax.set_title("\n".join(wrap(gene_symbols[i], 25)))
            
        except IndexError:
            ax.set_visible(False)
        except KeyError:
            ax.set_xticks(ticks)
            ax.set_xticklabels(colnames)
            ax.set_title("\n".join(wrap(gene_symbols[i], 25)))
            ax.set_xlabel(xlabel)
    for i in [1,3]:
        axes[i][0].set_ylabel('normalized mean counts +/- sem')
    for i in [0,1,2,3,4,5]:
        axes[4][i].set_xlabel(xlabel)
    
    plt.subplots_adjust(hspace=0.7)

    return fig, axes

#############################
## Preprocessing functions ##
#############################


def zumis_output2pandas(filename):
    """ UMIs is the previous step in the pipeline ran independently on the sequencing data, 
        it takes the original raw-data reads and maps them to a genome.
        code for reading zUMIS output to python modified from: https://github.com/sdparekh/zUMIs/wiki/Output """
    
    tag = re.search(r"(?<=\\IMM-).*?(?=.dgecounts.rds)", filename).group(0)
    
    mfunc = 'to_df <- function(dobj){return(as.data.frame(as.matrix(dobj)))}'
    rsparse2pandas = STAP(mfunc, "to_df")

    readRDS = robjects.r['readRDS']
    asMat = robjects.r['as.matrix']
    
    zumis_data = readRDS(filename)
    zd = dict(zip(zumis_data.names, list(zumis_data)))
    zd_exon=dict(zip(zd['exons'].names, list(zd['exons'])))
    counts = asMat(zd_exon['umicounts'])
    seq_data = robjects.conversion.ri2py(counts)
    
    seq_data = pd.DataFrame(   data = seq_data[:,:],
                               index = counts.rownames,
                               columns = [s + "_IMM-"+ tag for s in counts.colnames])
    
    return seq_data
