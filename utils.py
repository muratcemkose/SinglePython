
"""
Created on Mon Dec  3 14:45:33 2018

@author: Murat Cem Köse
"""
import numpy as np
import scanpy.api as sc
import pandas as pd
import scipy

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

def readDataDGE(path):
    data = sc.read_text(path,
                delimiter = "\t", first_column_names = True).transpose()
    
    data.var_names_make_unique()
    
    sc.pp.filter_cells(data, min_genes=1) # Hack to generate the n_genes column
    sc.pp.filter_genes(data, min_cells=1)
    
    return data
    
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
    result = sc.read(path + 'matrix.mtx').transpose() #, cache=True
    result.var_names = np.genfromtxt(path + 'genes.tsv', dtype=str)[:, 1]
    result.obs_names = np.genfromtxt(path + 'barcodes.tsv', dtype=str)
    result.var_names_make_unique()
    result.obs['n_counts'] = np.sum(result.X, axis=1).A1
    sc.pp.filter_cells(result, min_genes=min_genes)
    
    return result

def convertAnnDataToDf(scData):    
    try:
        result = pd.DataFrame(scData.X.toarray()) # If data is 10x data
    except:
        result = pd.DataFrame(scData.X[:]) # If data is Digital Gene expression matrix
        
    result.index = scData.obs_names.values
    result.columns = scData.var_names.values
    return result.T

def majorCellTypeCount(annotationResult):
    
    def splitAndReturnFirst(x):
        return x.split(":")[0]

    majorCellTypes = list(set(map(splitAndReturnFirst, annotationResult.final_annotations.transpose()["final_annotations"].tolist())))
    majorCellTypes.sort()

    result = pd.DataFrame()

    for i in majorCellTypes:
        totalCells = len(annotationResult.final_annotations.columns)
        cellTypeCount = len(annotationResult.final_annotations.transpose()["final_annotations"][annotationResult.final_annotations.transpose()["final_annotations"].str.startswith(i)])

        result = result.append([{
            "cellType" : i,
            "cellTypeCount" : cellTypeCount,
            "percentageOfCellType" : cellTypeCount / totalCells
        }])


    return result.set_index("cellType").sort_values("cellTypeCount", ascending = False)
