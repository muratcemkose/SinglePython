"""
Created on Mon Dec  3 14:22:31 2018

@author: Murat Cem Köse
"""


import SinglePython
import pandas as pd
import numpy as np
import utils
import scipy
import tuningMulti

def multiTuning_n(sc_data,refDataset,annot,top_labels,de,d):
    """ The wrapper function for multi processing.
    Parameters
    ----------
    sc_data : DataFrame
        Sc-RNAseq data.
        
    refDataset : DataFrame
        The reference dataset gene expression matrix.
        
    annot : DataFrame
        Annotations for each column in ref_data.
        
    top_labels: List
        Most correlated cell types from the previous round.
        
    de : Dict
        Differentially expressed genes for each combination of cell types.
        
    d : Dict
        The dictionary, containing the top n scoring cell types.
        
    Returns
    -------
    top_annotations : List
        A list of cell types that are associated with given cells.
        
    """
    cols=[j for j in d.keys() if list(d.get(j))==top_labels]
    return tuningMulti._FineTuneRoundByN(sc_data,refDataset,annot,top_labels,de,cols)

def multiTuning_t(sc_data,refDataset,annot,de,i,threshold,scores):
    """ The wrapper function for multi processing.
    Parameters
    ----------
    sc_data : DataFrame
        Sc-RNAseq data.
        
    refDataset : DataFrame
        The reference dataset gene expression matrix.
        
    annot : DataFrame
        Annotations for each column in ref_data.
        
    de : Dict
        Differentially expressed genes for each combination of cell types.
        
    i : String
        The cell name to calculate correlations for.
    
    threshold : Float
        The cutoff value for correlation scores of cell types to choose for the first fine tuning.
        
    scores : DataFrame
        Correlation scores from the first round.
        
    Returns
    -------
    top_annotations : List
        A list of cell types that are associated with given cells.
        
    """
    top_labels=scores[scores[i]>max(scores[i])-threshold].dropna().index.values
    while(len(top_labels)>1):
        top_annotations=tuningMulti._FineTuneRoundByT(sc_data,refDataset,annot,top_labels,de,i,threshold)
    return top_annotations

