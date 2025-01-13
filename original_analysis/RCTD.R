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


library(SingleCellExperiment)
library(BayesSpace)
readVisium('/ix1/wchen/Shiyue/Projects/2023_06_Influ_Mouse_Lung_ST/RawData/SpaceRanger_output/02_WT_Naive/outs')->sce
naive_1->sce
as.data.frame(sce@colData@listData)->spatial_info
cell_ID_vec = spatial_info$spot
names(cell_ID_vec) = c(1:nrow(spatial_info))
spatial_locations_matrix = as.matrix(spatial_info[, c("imagerow", "imagecol")])
matrix(nrow = 4024,ncol = 2)->k
k[,1]<-as.numeric(spatial_locations_matrix[,1])
k[,2]<-as.numeric(spatial_locations_matrix[,2])
knn_spatial <- dbscan::kNN(x = k,
                           k = 6)
spatial_locations_matrix<-as.data.frame(spatial_locations_matrix,row.names = spatial_info$spot)
as.data.frame(RCTD[[3]]@results[["weights"]])->celltype_weight
knn_spatial[['id']]->neigh
neigh_fi<-data.frame(neigh,row.names = spatial_info$spot)
for (i in 1:4024) {
  for (j in 1:6) {
    rownames(spatial_locations_matrix)[neigh[i,j]]->neigh_fi[i,j]
  }
}
merge(neigh_fi,celltype_weight, ,by='row.names',all=FALSE)->merge
rownames(merge)<-merge$Row.names
merge<-merge[,-1]
neigh_weight<-matrix(nrow = 4023,ncol=38)
for (i in 1:4023) {
  merge[which(row.names(merge)%in%merge[i,1:6]),]->s
  for (j in 7:36) {
    max(s[,j])->neigh_weight[i,j]
  }
}
neigh_weight<-as.data.frame(neigh_weight[,7:36],row.names = rownames(merge))
colnames(neigh_weight)<-paste(colnames(merge)[7:36],'neighbor',sep='_')
cor_mtx<-merge(merge[,7:36],neigh_weight, by='row.names', all=FALSE)
rownames(cor_mtx)<-cor_mtx$Row.names
cor_mtx<-cor_mtx[,-1]
pheatmap::pheatmap(cor(cor_mtx))->plot
pheatmap::pheatmap(cor(celltype_weight))->plot2
save_pheatmap_pdf <- function(x, filename, width=8, height=13) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(plot10, "NMF_plot10.pdf")

#NMF
library(NMF)
as.data.frame(t(F_1[['conv']]@counts))->celltype_weight
as.data.frame(lung_combine[['conv']]@counts)[,which(lung_combine@meta.data$experiment=='N_1')]->celltype_weight
NMF::nmf(celltype_weight,rank=3:10,seed=500)->nmr_res
nmfEstimateRank(celltype_weight, 6:10,seed=100)->nmf_es

as.data.frame(RCTD[[3]]@results[["weights"]])->celltype_weight
nmf(celltype_weight,rank=8,seed=100)->nmf_res_8
V.hat<-fitted(nmf_res_8) 
h <- coef(nmf_res_8)
w<-basis(nmf_res_8)
rownames(h)<-1:8
pheatmap::pheatmap(t(h),cluster_cols = FALSE)
F_1[['colo']]<-CreateAssayObject(counts = as.data.frame(t(w),row.names = as.character(1:8)))
SpatialFeaturePlot(F_1,'6')
colnames(spa[,which(as.data.frame(spa@assays[["colo"]]@counts)[3,]>6)])

library(SingleCellExperiment)
library(BayesSpace)


  dirname<-'/ix1/wchen/Shiyue/Projects/2023_06_Influ_Mouse_Lung_ST/RawData/SpaceRanger_output/09_WT_F_S/outs'
  spatial_dir <- file.path(dirname, "spatial")
  matrix_dir <- file.path(dirname, "filtered_feature_bc_matrix")
  
  if (!dir.exists(matrix_dir))
    stop("Matrix directory does not exist:\n  ", matrix_dir)
  if (!dir.exists(spatial_dir))
    stop("Spatial directory does not exist:\n  ", spatial_dir)
  
  colData <- read.csv(file.path(spatial_dir, "tissue_positions.csv"), header=FALSE)
  
  colnames(colData) <- c("spot", "in_tissue", "row", "col", "imagerow", "imagecol")
  rownames(colData) <- colData$spot
  colData <- colData[colData$in_tissue > 0, ]
  
  rowData <- read.table(file.path(matrix_dir, "features.tsv.gz"), header=FALSE)
  colnames(rowData) <- c("gene_id", "gene_name", "feature_type")
  rowData <- rowData[, c("gene_id", "gene_name")]
  rownames(rowData) <- scater::uniquifyFeatureNames(rowData$gene_id, rowData$gene_name)
  
  counts <- Matrix::readMM(file.path(matrix_dir, "matrix.mtx.gz"))
  barcodes <- read.table(file.path(matrix_dir, "barcodes.tsv.gz"), header=FALSE)
  colnames(counts) <- barcodes$V1
  rownames(counts) <- rownames(rowData)
  colData<-colData[-1,]
  counts <- counts[, rownames(colData)]
  
  sce <- SingleCellExperiment(assays=list(counts=counts),
                              rowData=rowData,
                              colData=colData)
  
  metadata(sce)$BayesSpace.data <- list()
  metadata(sce)$BayesSpace.data$platform <- "Visium"
  metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
  sce->F_S_2



