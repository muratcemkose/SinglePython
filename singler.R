library(SingleR)
library(Seurat)
start.time <- Sys.time()
ref_data=read.csv("C:/Users/murat_gga8ya6/Desktop/Thesis/ref_sample_filtered.csv",header = T,row.names = "genes")
types=read.csv("C:/Users/murat_gga8ya6/Desktop/Thesis/anot_sample.csv",header = T,row.names = "X")$x

#singler = CreateSinglerSeuratObject("C:/Users/murat_gga8ya6/Desktop/Thesis/Datasets/b_cells_filtered_gene_bc_matrice/" , NULL, "sample",  technology = "10x", species = "Human", ref.list = CreateVariableGeneSet(ref_data,types,n=20), fine.tune = F)

#the steps for extracting correlations.

sc.data = ReadSingleCellData("C:/Users/murat_gga8ya6/Desktop/Thesis/Datasets/b_cells_filtered_gene_bc_matrice/" , NULL)

N = colSums(sc.data$counts>0)
min.genes=500
sc_data= sc.data$counts[,N>=min.genes]

#sc.data.gl =read.csv("C:/Users/murat_gga8ya6/Desktop/Thesis/sample_sc_filtered.csv",header = T,row.names = "genes")
#res=SingleR("single",sc_data,as.matrix(ref_data),types,fine.tune = F,method = "de")

ref_data=as.matrix(ref_data)
sc_data=as.matrix(sc_data)
rownames(ref_data) = tolower(rownames(ref_data))
rownames(sc_data) = tolower(rownames(sc_data))
A = intersect(rownames(ref_data),rownames(sc_data))
sc_data = as.matrix(sc_data[A,])
ref_data = ref_data[A,]

not.use = rowSums(is.na(ref_data))>0 | rowSums(is.na(sc_data))>0 | rowSums(ref_data)==0
ref_data = ref_data[!not.use,]
sc_data = sc_data[!not.use,]
mat = medianMatrix(ref_data,types)

n = round(500*(2/3)^(log2(c(ncol(mat)))))
genes.filtered = unique(unlist(unlist(lapply(1:ncol(mat), function(j) {
  lapply(1:ncol(mat), function(i) {
    s=sort(mat[,j]-mat[,i],decreasing=T);
    s=s[s>0];
    names(s)[1:min(n,length(s))]
  })}))))[-1]
print(paste0("Number of DE genes:", length(genes.filtered)))
cell.names = colnames(sc_data)
genes=genes.filtered
sc_data = as.matrix(sc_data[genes,])
ref_data = as.matrix(ref_data[genes,])
step=10000
quantile.use=0.8
r=cor(sc_data,ref_data,method='spearman')
agg_scores = quantileMatrix(r,types,quantile.use)
labels = colnames(agg_scores)[max.col(agg_scores)]
output = list()

names(labels)=t(colnames(sc_data))

output$scores = as.matrix(t(agg_scores))

output$labels = as.matrix(labels)
output$r = r
output$scores = t(output$scores)

end.time <- Sys.time()
print(end.time - start.time)