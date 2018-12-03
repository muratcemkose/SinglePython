
"""
Created on Mon Dec  3 14:49:07 2018

@author: Murat Cem Köse
"""

class ResultObject:
    def __init__(self,final_annotations,initial_scores,corelation_matrix,common_de_genes):
        """Contructor function for SinglePython class.
    
        Parameters
        ----------
        final_annotations : DataFrame
            Final cell type annotations for each cell.
        
        initial_scores : DataFrame
            Initial scores of cell types for each cell.
            
        corelation_matrix : DataFrame
            Correlation matrix of cell types.
            
        common_de_genes : List
            Differentially highly expressed genes between each pair of cell types.
            
        """
        
        self.final_annotations=final_annotations
        self.initial_scores=initial_scores
        self.corelation_matrix=corelation_matrix
        self.common_de_genes=common_de_genes
        