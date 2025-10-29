library(Seurat)
#import the result of RCTD deconvolution
readRDS(lung_combine,'~/lung_combine.rds')
#perform single cell clustering
sc_pro<-function(object,harmony,resolution,group.by){
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object)
  object <- ScaleData(object)
  set.seed(100)
  object <- RunPCA(object)
  if (harmony==TRUE){
    set.seed(100)
    object <- RunHarmony(object,group.by.vars = group.by)
    set.seed(100)
    object <- FindNeighbors(object, reduction = 'harmony', dims = 1:30)
    set.seed(100)
    object <- FindClusters(object,resolution = resolution, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'harmony',dims = 1:30)
  }
  else {
    set.seed(100)
    object <- FindNeighbors(object, reduction = 'pca', dims = 1:30)
    set.seed(100)
    object <- FindClusters(object,resolution = resolution, verbose = FALSE)
    set.seed(100)
    object <- RunUMAP(object, reduction = 'pca',dims = 1:30)
  }
  return(object)
}
library(harmony)
#uncorrect batch
sc_pro(lung_combine,harmony = FALSE,resolution = 1)->lung_combine_before_batch
lung_combine$orig.ident<-substr(colnames(lung_combine),1,4)
#correct batch
sc_pro(lung_combine,harmony = TRUE,resolution = 1,group.by = 'orig.ident')->lung_combine_after_batch

#reorder the seurat clusters
from <- c(0,2,1,3,11,8,10,5,9,4,6,7)
to <- c(1:12)
library(plyr)
mapvalues(lung_combine_after_batch$seurat_clusters, from = from, to = to)->lung_combine_after_batch$seurat_clusters
lung_combine_after_batch$seurat_clusters<-factor(lung_combine_after_batch$seurat_clusters,levels = 1:12)





