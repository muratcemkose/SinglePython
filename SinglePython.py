
"""
Created on Mon Dec  3 14:22:31 2018

@author: Murat Cem Köse
"""

import Result
import utils
import tuning
import scipy
import numpy as np
import pandas as pd

class SinglePythonObject:
    def __init__(self, sc_location,refDataset,annot=None,fine_tuning=False,tuning_by="top_n",tuning_threshold=0.05,tuning_top_n=7,min_gene_th=500,de_genes_n=None):
        """Contructor function for SinglePython class.
    
        Parameters
        ----------
        sc_location : Str
            the location of sc-RNAseq data.
            
        refDataset : DataFrame
            The reference dataset gene expression matrix.
            
        annot : DataFrame
            Annotations for each column in ref_data.
            
        fine_tuning : Boolean
            If fine tuning will be applied.
            
        tuning_by : Str
            There are two tuning options. It can be ither done by taking top n most correlated cell types after
            the initial annotation (top_n) or setting a minimum threshold value to select cell types for tuning 
            (threshod). The default value is top_n.
            
        tuning_threshold : Float
            The threshold value for the selection of cell types for "threshold" tuning. Default value is 0.05.
            
        tuning_top_n : Int
            The number of top correlated cell types to be selected for "top_n" tuning. Default value is 7. 
            
        min_gene_th : Int
            Minimum gene threshold for genes that have expression more than zero. It is necessary to determine
            which cells to contribute to the analysis.
            
        de_genes_n : Int
            The number of top differentially expressed genes to be chosen for the analysis.
            
        """
        
        self.sc_data = utils.readData_SingleR(sc_location,min_gene_th)
        if annot is not None:
            self.refDataset = refDataset.astype(float)
            self.annot = annot
        else:
            self.annot=refDataset.loc["cellType"].T
            refDataset=refDataset.drop("cellType")
            self.refDataset=refDataset.astype(float)
        self.de_genes_n=de_genes_n
        self.fine_tuning=fine_tuning
        self.tuning_by=tuning_by
        self.tuning_threshold=tuning_threshold
        self.tuning_top_n=tuning_top_n
        
    def AnnotateCellTypes(self):
        """ Finds best annotation for single cells and Returns a Result object.
    
        """
        print(self.refDataset.index)
        intersect=np.intersect1d(self.refDataset.index.values,self.sc_data.index)
        sc_data=self.sc_data.loc[intersect]
        refDataset=self.refDataset.loc[intersect]
    
        de=utils.getDEgenes(refDataset,self.annot)
        de_merged=[]
        [de_merged.extend(i) for i in  de.values()]
        de_merged=np.unique(de_merged)
    
        cor=scipy.stats.spearmanr(sc_data.loc[de_merged],refDataset.loc[de_merged])
        cor=pd.DataFrame(cor[0]).iloc[:,0:len(sc_data.columns)][-len(refDataset.columns):]
        cor.columns=sc_data.columns
        cor.index=refDataset.columns
        cor["cellType"]=self.annot["cellType"].values
        scores=cor.groupby("cellType").quantile(q=0.8)
        if self.fine_tuning==True:
            if (self.tuning_by=="top_n"):
                final_annotations=tuning._FineTuneByN(sc_data,refDataset,self.annot,de,scores,self.tuning_top_n)
                return Result.ResultObject(final_annotations,scores,cor,de_merged)
            elif (self.tuning_by=="threshold"):
                final_annotations=tuning._FineTuneByT(sc_data,refDataset,self.annot,de,scores,self.tuning_threshold)
                return Result.ResultObject(final_annotations,scores,cor,de_merged)
            else:
                print("Undefined tuning method.")
        else:
            return Result.ResultObject(pd.DataFrame(scores.idxmax(),columns=["final_annotations"]),scores,cor,de_merged)
        
