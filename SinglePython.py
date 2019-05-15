
"""
Created on Mon Dec  3 14:22:31 2018

@author: Murat Cem Köse
"""

import utils
import tuning
import tuningMulti
import scipy
import numpy as np
import pandas as pd

class SinglePythonObject:
    def __init__(self, scData, refDataset, refAnnot, fine_tuning=False, tuning_by="top_n", tuning_threshold=0.05, tuning_top_n=8, multiProcess = True):
        """Contructor function for SinglePython class.
    
        Parameters
        ----------
        scData : AnnData
            the location of sc-RNAseq data.
            
        refDataset : DataFrame
            The reference dataset gene expression matrix.
            
        refAnnot : DataFrame
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
            
        multiProcess : Boolean
            If the program will use multiple processing.
            
        """
        
        self.scData = scData
        self.refDataset = refDataset.astype(float)
        self.refAnnot = refAnnot
        self.refAnnot.columns = ["cellType"]
        self.method_parameters = {}
        self.method_parameters.update({"fine_tuning":fine_tuning})        
        self.method_parameters.update({"tuning_by":tuning_by})
        self.method_parameters.update({"tuning_top_n":tuning_top_n}) 
        self.method_parameters.update({"tuning_threshold":tuning_threshold}) 
        self.method_parameters.update({"multiProcess":multiProcess}) 
    def annotateCellTypes(self):
        """ Finds best annotation for single cells and adds it on the  Result object.
    
        """
        de=utils.getDEgenes(self.refDataset,self.refAnnot)
        sc_data = utils.convertAnnDataToDf(self.scData)
        intersect=np.intersect1d(self.refDataset.index.values,sc_data.index)
        sc_data=sc_data.loc[intersect]
        refDataset=self.refDataset.loc[intersect]
        
        n=int(500*np.power(2/3,np.log2(len(np.unique(self.refAnnot.cellType)))))
       
        de_merged=[]
        [de_merged.extend(i[:n]) for i in  de.values()]
        de_merged = np.intersect1d(de_merged,sc_data.index)
        de_merged=np.unique(de_merged)
    
        cor=scipy.stats.spearmanr(sc_data.loc[de_merged],refDataset.loc[de_merged])
        cor=pd.DataFrame(cor[0]).iloc[:,0:len(sc_data.columns)][-len(refDataset.columns):]
        cor.columns=sc_data.columns
        cor.index=refDataset.columns
        cor["cellType"]=self.refAnnot["cellType"].values
        scores=cor.groupby("cellType").quantile(q=0.8)
        self.initial_scores = scores
        self.de_dict = de
        if self.method_parameters.get("fine_tuning")==True:
            if (self.method_parameters.get("tuning_by")=="top_n"):
                if self.method_parameters.get("multiProcess") == True:
                    self.scData.obs["cell_type"] = tuningMulti._FineTuneByN(sc_data,refDataset,self.refAnnot,de,scores,self.method_parameters.get("tuning_top_n")).iloc[0]
                else:
                    self.scData.obs["cell_type"] = tuning._FineTuneByN(sc_data,refDataset,self.refAnnot,de,scores,self.method_parameters.get("tuning_top_n")).iloc[0]
            elif (self.method_parameters.get("tuning_by")=="threshold"):
                if self.method_parameters.get("multiProcess") == True:
                    self.scData.obs["cell_type"] = tuningMulti._FineTuneByT(sc_data, refDataset, self.refAnnot, de, scores, self.method_parameters.get("tuning_threshold")).iloc[0]
                else:
                    self.scData.obs["cell_type"] = tuning._FineTuneByT(sc_data,refDataset,self.refAnnot,de,scores,self.method_parameters.get("tuning_threshold")).iloc[0]
            else:
                print("Undefined tuning method.")
        else:
            self.scData.obs["cell_type"] = pd.DataFrame(scores.idxmax(),columns=["final_annotations"]).iloc[0]
        
