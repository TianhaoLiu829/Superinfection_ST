library(Seurat)

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
readRDS('/ix1/wchen/liutianhao/result/lung_ST/alcorn/lung_combine.rds')->lung_combine

sc_pro(lung_combine,harmony = FALSE,resolution = 1)->lung_combine_before_batch
lung_combine$orig.ident<-substr(colnames(lung_combine),1,4)
sc_pro(lung_combine,harmony = TRUE,resolution = 1,group.by = 'orig.ident')->lung_combine_after_batch

saveRDS(lung_combine_before_batch,'/ix1/wchen/liutianhao/result/lung_ST/alcorn/lung_combine_before_batch.rds')
saveRDS(lung_combine_after_batch,'/ix1/wchen/liutianhao/result/lung_ST/alcorn/lung_combine_after_batch.rds')






