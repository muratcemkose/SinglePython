
"""
Created on Mon Dec  3 14:45:33 2018

@author: Murat Cem Köse
"""
import numpy as np
import scanpy.api as sc
import pandas as pd
import os
import scipy
import seaborn as sns
from matplotlib import pyplot as plt

def getDEgenes(refDataset,annot,n=200):
    """ Creates a dictionary for differentially highly expressed genes for all pairwise cell types in the a reference data set.
    
    Parameters
    ----------
    refDataset : DataFrame
        The reference dataset gene expression matrix.
        
    annot : DataFrame
        Annotations for each column in ref_data.
        
    n : Int
        Number of top DE genes to be kept.
    
    Returns
    -------
    deGenes : dict with multiple index
        Dictionary containing differentially highly expressed genes for each combination of cell types. 
    """
    types = annot.groupby("cellType").groups.keys()
    median = refDataset.groupby(annot.cellType.values,axis=1).median()
    medDiff = {}
    [medDiff.update({(i,j):median[i]-median[j]}) for i in types  for j in types if i!=j]
    deGenes={}
    [deGenes.update({i: median.index[medDiff.get(i).argsort()[0:200]]}) for i in medDiff.keys()]
    return deGenes
    
def readCountMartix(path,min_genes):
    """ 
    Reads, precesses and returns single cell data from a count matrix. 
    
    Parameters
    ----------
    path : Str
        Diractory path, the location of single cell data.
        
    min_genes : Int
        The minimum number of genes for a cell to have in order to participate the analysis.
        
    Returns
    -------
    scData : AnnData
        Single cell data. 
        
    """
    
    result = sc.read_csv(path).transpose() #, cache=True
    result.var_names_make_unique()
    n_counts = np.sum(result.X, axis=1)
    result.obs['n_counts'] = n_counts
    sc.pp.filter_cells(result, min_genes=min_genes)
    
    return result

def read10xData(path,min_genes):
    """ 
    Reads, precesses and returns single cell data from 10X mtx file. 
    
    Parameters
    ----------
    path : Str
        Diractory path, the location of single cell data.
        
    min_genes : Int
        The minimum number of genes for a cell to have in order to participate the analysis.
        
    Returns
    -------
    scData : AnnData
        Single cell data. 
        
    """
    
    result = sc.read(path + 'matrix.mtx').transpose() #, cache=True
    result.var_names = np.genfromtxt(path + 'genes.tsv', dtype=str)[:, 1]
    result.obs_names = np.genfromtxt(path + 'barcodes.tsv', dtype=str)
    result.var_names_make_unique()
    result.obs['n_counts'] = np.sum(result.X, axis=1).A1
    sc.pp.filter_cells(result, min_genes=min_genes)
    
    return result

def readSCData(path,min_genes):
    """ 
    The wrapper function for the program to read single cell dataset. 
    
    Parameters
    ----------
    path : Str
        Diractory path, the location of single cell data.
        
    min_genes : Int
        The minimum number of genes for a cell to have in order to participate the analysis.
        
    Returns
    -------
    scData : AnnData
        Single cell data. 
        
    """

    if "matrix.mtx" in os.listdir(path):
        result = read10xData(path, min_genes)
        return result
    else:
        try:
            data_file = [i for i in os.listdir(path) if "expression" in i][0]
            result = readCountMartix(path+data_file, min_genes)
            return result
        except:
            print("Please end your single cell data count matrix file with word 'expression'.")


def convertAnnDataToDf(scData):    
    """ 
    Converts AnnData object into  a DataFrame. 
    
    Parameters
    ----------
    scData : AnnData
        Single cell data. 
        
    Returns
    -------
    scData : DataFrame
        Single cell data. 
        
    """
    try:
        result = pd.DataFrame(scData.X.toarray()) # If data is 10x data
    except:
        result = pd.DataFrame(scData.X[:]) # If data is Digital Gene expression matrix
        
    result.index = scData.obs_names.values
    result.columns = scData.var_names.values
    return result.T


    
