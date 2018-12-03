
"""
Created on Mon Dec  3 14:45:33 2018

@author: Murat Cem Köse
"""
import numpy as np
import scanpy.api as sc
import pandas as pd

def getDEgenes(refDataset,annot=None,n=None):
    """ Creates a dictionary for differentially highly expressed genes for all pairwise cell types in the a reference data set.
    
    Parameters
    ----------
    refDataset : DataFrame
        The reference dataset gene expression matrix.
        
    annot : DataFrame
        Annotations for each column in ref_data.
        
    n : Int
        Number of top differential genes to select.
        
    Returns
    -------
    deGenes : dict with multiple index
        Dictionary containing differentially highly expressed genes for each combination of cell types. 
    """
    
    if n is None:
        n=int(500*np.power(2/3,np.log2(len(np.unique(annot.cellType)))))
        
    types=np.unique(annot.cellType)
    median=refDataset.groupby(annot.cellType.values,axis=1).apply(np.median,axis=1)
    deGenes={}
    [deGenes.update({(i,j):median[i]-median[j]}) for i in median.index  for j in median.index if i!=j]
    for i in deGenes.keys():
        deGenes[i]=refDataset.iloc[deGenes.get(i).argsort()[-n:]].index.values.tolist()
    return deGenes

def readData_SingleR(path,min_genes):
    """ Reads, precesses and returns single cell data as a matrix. 
    
    Parameters
    ----------
    path : Str
        Diractory path, the location of single cell data.
        
    min_genes : Int
        The minimum number of genes for a cell to have in order to participate the analysis.
        
    Returns
    -------
    sc_data : DataFrame
        Single cell data matrix. 
        
    """
    data = sc.read(path + 'matrix.mtx').transpose() #, cache=True
    data.var_names = np.genfromtxt(path + 'genes.tsv', dtype=str)[:, 1]
    data.obs_names = np.genfromtxt(path + 'barcodes.tsv', dtype=str)
    data.var_names_make_unique()
    data.obs['n_counts'] = np.sum(data.X, axis=1).A1
    sc.pp.filter_cells(data, min_genes=min_genes)

    sc_data=pd.DataFrame(data.X.toarray())
    sc_data.index=data.obs_names.values
    sc_data.columns=data.var_names.values
    sc_data=sc_data.T
    return sc_data
