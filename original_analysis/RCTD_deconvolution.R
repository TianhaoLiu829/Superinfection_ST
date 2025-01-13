library(spacexr)
library(Seurat)

ref <- readRDS("/ix1/wchen/liutianhao/result/lung_ST/GSE202325_ref.rds")
cluster <- as.factor(ref$celltype2)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(ref@assays[["RNA"]]@counts, cluster, nUMI)



# set up query with the RCTD function SpatialRNA
N_1<-Load10X_Spatial('/ix1/wchen/Shiyue/Projects/2023_06_Influ_Mouse_Lung_ST/RawData/SpaceRanger_output/01_WT_Naive/outs')
N_2<-Load10X_Spatial('/ix1/wchen/Shiyue/Projects/2023_06_Influ_Mouse_Lung_ST/RawData/SpaceRanger_output/02_WT_Naive/outs')
S_1<-Load10X_Spatial('/ix1/wchen/liutianhao/result/03_WT_SA_Only_2/outs')
S_2<-Load10X_Spatial('/ix1/wchen/liutianhao/result/05_WT_SA_Only_2/outs')
F_1<-Load10X_Spatial('/ix1/wchen/Shiyue/Projects/2023_06_Influ_Mouse_Lung_ST/RawData/SpaceRanger_output/06_WT_Flu_Only/outs')
F_2<-Load10X_Spatial('/ix1/wchen/Shiyue/Projects/2023_06_Influ_Mouse_Lung_ST/RawData/SpaceRanger_output/07_WT_Flu_Only/outs')
FS_1<-Load10X_Spatial('/ix1/wchen/Shiyue/Projects/2023_06_Influ_Mouse_Lung_ST/RawData/SpaceRanger_output/08_WT_F_S/outs')
FS_2<-Load10X_Spatial('/ix1/wchen/Shiyue/Projects/2023_06_Influ_Mouse_Lung_ST/RawData/SpaceRanger_output/09_WT_F_S/outs')

RCTD<-list()
j<-1
for (i in c(N_1,N_2,S_1,S_2,F_1,F_2,FS_1,FS_2)) {
  counts <- i@assays[["Spatial"]]@counts
  coords <- GetTissueCoordinates(i)
  colnames(coords) <- c("x", "y")
  coords[is.na(colnames(coords))] <- NULL
  query <- SpatialRNA(coords, counts, colSums(counts))
  RCTD[[j]] <- create.RCTD(query, reference, max_cores = 29)
  RCTD[[j]] <- run.RCTD(RCTD[[j]], doublet_mode = "doublet")
  j<-j+1
}
save.image('/ix1/wchen/liutianhao/Work.RData')


#incooperate the RCTD result to seurat object
process_RCTD<-function(object,order){
  object<-object[,rownames(RCTD[[order]]@results[["weights"]])]
  object[['conv']]<-CreateAssayObject(t(as.matrix(RCTD[[order]]@results[["weights"]])))
  AddMetaData(object,RCTD[[order]]@results[["results_df"]])->object
  return(object)
}
c(N_1,N_2,S_1,S_2,F_1,F_2,FS_1,FS_2)->union
lapply(1:8,function(x){return(process_RCTD(union[[x]],x))})->union
#combine samples
library(SeuratDisk)
lung_combine<-merge(union[[1]],unlist(union[2:8]),add.cell.ids = c("N_1","N_2","S_1","S_2","F_1","F_2","FS_1","FS_2"))
dim(lung_combine)
saveRDS(lung_combine,'~/lung_combine.rds')





