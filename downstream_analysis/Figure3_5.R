
make_cell_propor_box<-function(target,object){
  df->df2
  # only include spots with proportion > 0.0005
  if (object%in%c('Alveolar-macrophage','NK','cDC','Plasma-cell')){
    df2<-df2[which(df[,object]>0.001&df2[,object]<max(df2[,object])  ),]
    lung_combine_after_batch2<-lung_combine_after_batch[,which(df[,object]>0.001&df[,object]<max(df[,object]) )]
  } else if (object=='Neutrophil') {
    df2<-df2[which(df[,object]>0.0005&df[,object]<max(df[,object])  ),]
    lung_combine_after_batch2<-lung_combine_after_batch[,which(df[,object]>0.0005&df[,object]<max(df[,object]) )]
  } else {
    lung_combine_after_batch2<-lung_combine_after_batch
  }
  df2[which(lung_combine_after_batch2@meta.data[,target]==TRUE&lung_combine_after_batch2$orig.ident=='N_1'),object]->prop_macro_airway_N_1
  df2[which(lung_combine_after_batch2@meta.data[,target]==TRUE&lung_combine_after_batch2$orig.ident=='N_2'),object]->prop_macro_airway_N_2
  df2[which(lung_combine_after_batch2@meta.data[,target]==TRUE&lung_combine_after_batch2$orig.ident=='S_1'),object]->prop_macro_airway_S_1
  df2[which(lung_combine_after_batch2@meta.data[,target]==TRUE&lung_combine_after_batch2$orig.ident=='S_2'),object]->prop_macro_airway_S_2
  df2[which(lung_combine_after_batch2@meta.data[,target]==TRUE&lung_combine_after_batch2$orig.ident=='F_1'),object]->prop_macro_airway_F_1
  df2[which(lung_combine_after_batch2@meta.data[,target]==TRUE&lung_combine_after_batch2$orig.ident=='F_2'),object]->prop_macro_airway_F_2
  df2[which(lung_combine_after_batch2@meta.data[,target]==TRUE&lung_combine_after_batch2$orig.ident=='FS_1'),object]->prop_macro_airway_FS_1
  df2[which(lung_combine_after_batch2@meta.data[,target]==TRUE&lung_combine_after_batch2$orig.ident=='FS_2'),object]->prop_macro_airway_FS_2
  data <- data.frame(
    value = c(prop_macro_airway_N_1,prop_macro_airway_N_2,prop_macro_airway_S_1,prop_macro_airway_S_2,prop_macro_airway_F_1, prop_macro_airway_F_2, prop_macro_airway_FS_1, prop_macro_airway_FS_2),
    group = factor(rep(c(paste0("prop_",object,"_airway_N_1"),paste0("prop_",object,"_airway_N_2"),paste0("prop_",object,"_airway_S_1"),paste0("prop_",object,"_airway_S_2"),paste0("prop_",object,"_airway_F_1"), paste0("prop_",object,"_airway_F_2"), paste0("prop_",object,"_airway_FS_1"),paste0("prop_",object,"_airway_FS_2") ), 
                       times = c(length(prop_macro_airway_N_1),length(prop_macro_airway_N_2),length(prop_macro_airway_S_1),length(prop_macro_airway_S_2),length(prop_macro_airway_F_1), length(prop_macro_airway_F_2), length(prop_macro_airway_FS_1), length(prop_macro_airway_FS_2))))
  )
  data$condition<-str_sub(data$group,-4)
  data$condition<-gsub('^_','',data$condition)
  data$condition2<-gsub('_.*','',data$condition)
  data$group<-factor(data$group,levels = c(paste0("prop_",object,"_airway_N_1"),paste0("prop_",object,"_airway_N_2"),paste0("prop_",object,"_airway_S_1"),paste0("prop_",object,"_airway_S_2"),paste0("prop_",object,"_airway_F_1"), paste0("prop_",object,"_airway_F_2"), paste0("prop_",object,"_airway_FS_1"),paste0("prop_",object,"_airway_FS_2")))
  ggboxplot(data, x = "group", y = "value",fill = 'condition2') +
    stat_compare_means(comparisons = list(c("prop_macro_airway_F_1", "prop_macro_airway_F_2"), c("prop_macro_airway_F_1", "prop_macro_airway_FS_1"), c("prop_macro_airway_F_1", "prop_macro_airway_FS_2"),
                                          c(paste0("prop_",object,"_airway_F_2"), "prop_macro_airway_FS_1"), c("prop_macro_airway_F_2", "prop_macro_airway_FS_2"), c("prop_macro_airway_FS_1", "prop_macro_airway_FS_2")), 
                       method = "wilcox.test") +
    theme_bw() +
    scale_x_discrete(labels =unique(data$condition))+
    theme(
      axis.line = element_line(size = 1.3), # Add axis lines
      axis.text.x = element_text(size = 18),  # Adjust x-axis text size
      axis.text.y = element_text(size = 18),  # Adjust y-axis text size
      axis.title.x = element_text(size = 20), # Adjust x-axis title size
      axis.title.y = element_text(size = 20),  # Adjust y-axis title size
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()  
    )+
    labs(title = "Boxplot with Significance Annotations", x = "Group", y = object)->p
  return(p)
}


#load data
lung_combine_after_batch <- readRDS("/ix1/wchen/liutianhao/result/lung_ST/alcorn/script/object/lung_combine_after_batch.rds")
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(Seurat)
df<-as.data.frame(t(as.data.frame(lung_combine_after_batch[['conv']]@counts)))
df$all_airway<-rowMeans(df[,c('Ciliated-cell','Club-cell','Goblst-cell')])
immune_manual_list<-list()
for (sample in list.files('/ix1/wchen/liutianhao/result/lung_ST/alcorn/script/revision/immune_manual/')){
  read.csv(paste0('/ix1/wchen/liutianhao/result/lung_ST/alcorn/script/revision/immune_manual/',sample))->a
  a$Barcode<-paste0(gsub('immune_','',gsub('.csv','',sample)) ,'_',a$Barcode)
  a->immune_manual_list[[sample]]
}
do.call(rbind,immune_manual_list)->immune_manual_list
lung_combine_after_batch$immune_manual<-FALSE
lung_combine_after_batch$immune_manual[which(colnames(lung_combine_after_batch)%in%immune_manual_list$Barcode)]<-TRUE
immune_mtx<-df[,c('B-cell','Alveolar-macrophage','Interstitial-macrophage','Cd8-tcell','Cd4-tcell','cDC','IgA-plasma-cell','Macrophage','Mast-cell','mitotic-tcell','Monocyte','NK','Plasma-cell','Neutrophil')]
#immune_mtx$`Interstitial-macrophage`<-immune_mtx$`Interstitial-macrophage`-0.1
#immune_mtx$Macrophage<-immune_mtx$Macrophage-0.02
#df$all_immune<-log1p(rowMeans(apply(immune_mtx,2,function(x) x/mean(x)) ))
#df$all_immune<-log(rowMeans(immune_mtx)+1,base = 10)
df$all_immune<-rowMeans(immune_mtx)
df$immune_manual<-lung_combine_after_batch$immune_manual
df$airway<-lung_combine_after_batch$airway
df$region2<-'Parenchyla_non_immune'
df$region2[which(lung_combine_after_batch$near_airway)]<-'near_airway'
df$region2[which(!lung_combine_after_batch$near_airway&lung_combine_after_batch$immune_manual)]<-'Parenchyma_immune'
df$region2<-factor(df$region2,levels = unique(df$region2)[c(1,3,2)])
ggboxplot(df, x = "region2", y = "all_airway",fill = 'region2',palette=scales::hue_pal()(3)[c(3,1,2)])+
  theme_bw() +
  theme(
    axis.line = element_line(size = 1.3), # Add axis lines
    axis.text.x = element_text(size = 18,angle = 45, hjust = 1),  # Adjust x-axis text size
    axis.text.y = element_text(size = 18),  # Adjust y-axis text size
    axis.title.x = element_text(size = 20), # Adjust x-axis title size
    axis.title.y = element_text(size = 20),  # Adjust y-axis title size
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )+
  labs(title = "Boxplot with Significance Annotations", x = "Group", y = 'all_immune')->p
ggsave('p1.png',p,width = 7,height = 8,dpi = 300)
ggboxplot(df[which(!lung_combine_after_batch$airway),], x = "region2", y = "all_airway",fill = 'region2',palette=scales::hue_pal()(3)[c(3,1,2)])+
  theme_bw() +
  theme(
    axis.line = element_line(size = 1.3), # Add axis lines
    axis.text.x = element_text(size = 18,angle = 45, hjust = 1),  # Adjust x-axis text size
    axis.text.y = element_text(size = 18),  # Adjust y-axis text size
    axis.title.x = element_text(size = 20), # Adjust x-axis title size
    axis.title.y = element_text(size = 20),  # Adjust y-axis title size
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )+
  
  labs(title = "Boxplot with Significance Annotations", x = "Group", y = 'all_immune')->p
ggsave('p1.png',p,width = 7,height = 8,dpi = 300)
lung_combine_after_batch$region2<-df$region2
lung_combine_after_batch$region_final<-'Parenchyma_non_immune'
lung_combine_after_batch$region_final[which(lung_combine_after_batch$immune)]<-'Parenchyma_immune'
lung_combine_after_batch$region_final[which(lung_combine_after_batch$airway)]<-'airway'
lung_combine_after_batch$region_final<-factor(lung_combine_after_batch$region_final,levels = unique(lung_combine_after_batch$region_final)[c(3,1,2)])
SpatialDimPlot(lung_combine_after_batch,group.by = 'region_final',images = 'slice1.5')
SpatialDimPlot(lung_combine_after_batch,group.by = 'region2',images = 'slice1.5')
SpatialDimPlot(lung_combine_after_batch,group.by = 'region_final',pt.size.factor = 1.3)->p
ggsave('p.png',p,width = 30,height = 17,dpi = 300)

DimPlot(lung_combine_after_batch,group.by = c('region2','region_final'),cols=scales::hue_pal()(3)[c(3,1,2)])

lung_combine_after_batch$near_airway<-FALSE
lung_combine_after_batch$near_airway[which(lung_combine_after_batch$distance_min<200)]<-TRUE
lung_combine_after_batch$distal<-FALSE
lung_combine_after_batch$distal[which(lung_combine_after_batch$distance_min>200)]<-TRUE
lung_combine_after_batch$whole<-TRUE
t.test(data$value[which(data$condition2=='S')],data$value[which(data$condition2=='F')])

#find some differential LR pairs
lung_combine_after_batch[,which(lung_combine_after_batch$airway==TRUE)]->airway
DefaultAssay(airway)<-'sender'
airway$IM_prop<-airway@assays$conv@counts['Interstitial-macrophage',]
FindMarkers(airway,group.by = 'orig.ident',ident.1 = c('FS_1','FS_2'),ident.2 = c('F_1','F_2'),test.use = 'LR',latent.vars = c('immune_prop','IM_prop'))->double_marker_send_flu
FindMarkers(airway,group.by = 'orig.ident',ident.1 = c('FS_1','FS_2'),ident.2 = c('S_1','S_2'),test.use = 'LR',latent.vars = c('immune_prop','IM_prop'))->double_marker_send_bac
DefaultAssay(airway)<-'reciever'
FindMarkers(airway,group.by = 'orig.ident',ident.1 = c('FS_1','FS_2'),ident.2 = c('F_1','F_2'),test.use = 'LR',latent.vars = c('immune_prop','IM_prop'))->double_marker_recei_flu
FindMarkers(airway,group.by = 'orig.ident',ident.1 = c('FS_1','FS_2'),ident.2 = c('S_1','S_2'),test.use = 'LR',latent.vars = c('immune_prop','IM_prop'))->double_marker_recei_bac

#Find the enrichemnt of cell for LR pairs
df_lr1<-lung_combine_after_batch@assays$sender[c('s.Cxcl5.Cxcr2','s.Ccl5.Ccr5'),]
df_lr2<-lung_combine_after_batch@assays$reciever[c('r.Cxcl1.Cxcr2','r.TNF','r.Ccl2.Ccr2','r.VCAM'),]
rbind(df_lr1,df_lr2)->df_lr_all
df_lr_all<-as.data.frame(lung_combine_after_batch@assays$sender$counts)+as.data.frame(lung_combine_after_batch@assays$reciever$counts)
df_lr_all<-df_lr_all[c('s.Cxcl13.Cxcr5','s.Ccl7.Ccr2','s.Ccl5.Ccr5','s.Cxcl10.Cxcr3','s.Ccl8.Ccr5','s.Ccl8.Ccr2','s.Cxcl9.Cxcr3','s.Ccl2.Ccr2','s.Cxcl12.Ackr3','s.Cxcl12.Cxcr4'),which(lung_combine_after_batch$ibalt)]
gsub('^.[.]','',rownames(df_lr_all))->rownames(df_lr_all)
df_cell<-lung_combine_after_batch@assays$conv[c('Neutrophil','Interstitial-macrophage','Alveolar-macrophage','Monocyte','NK','cDC','Endothelial-cell','Cd274-endothelial-cell','Lymphatic-endothelial-cell'),]
#df_cell<-lung_combine_after_batch@assays$conv[c('Neutrophil','Interstitial-macrophage','Alveolar-macrophage','Monocyte','NK','cDC','Endothelial-cell','Cd274-endothelial-cell','Lymphatic-endothelial-cell','Cd4-tcell','Cd8-tcell','B-cell'),which(lung_combine_after_batch$ibalt)]
#df_cell<-lung_combine_after_batch@assays$conv[c('NK','cDC','Endothelial-cell','Cd274-endothelial-cell','Lymphatic-endothelial-cell','Cd4-tcell','Cd8-tcell','B-cell'),which(lung_combine_after_batch$ibalt)]
#rownames(df_cell)[4]<-'CD274-endothelial-cell'
rownames(df_cell)[8]<-'CD274-endothelial-cell'

cor_matrix <- apply(df_lr_all, 1, function(x) apply(df_cell, 1, function(y) cor(x, y)))
pheatmap::pheatmap(cor_matrix,scale = 'row')->p
library(ggplot2)
ggsave('p.png',p,width = 6,height = 5)

#Find the enrichemnt of cell for LR pairs
df_lr1<-lung_combine_after_batch@assays$sender[c('s.Cxcl5.Cxcr2','s.Ccl5.Ccr5'),]
df_lr2<-lung_combine_after_batch@assays$reciever[c('r.Cxcl1.Cxcr2','r.TNF','r.Ccl2.Ccr2'),]
rbind(df_lr1,df_lr2)->df_lr_all
gsub('^.[.]','',rownames(df_lr_all))->rownames(df_lr_all)
df_cell<-lung_combine_after_batch@assays$conv[c('Neutrophil','Interstitial-macrophage','Alveolar-macrophage','Monocyte','NK','cDC','Endothelial-cell','Cd274-endothelial-cell','Lymphatic-endothelial-cell'),]
rownames(df_cell)[8]<-'CD274-endothelial-cell'

cor_matrix <- apply(df_lr_all, 1, function(x) apply(df_cell, 1, function(y) cor(x, y)))
pheatmap::pheatmap(cor_matrix,scale = 'column')->p
library(ggplot2)
ggsave('p.png',p,width = 7,height = 5)


#Figure3A
make_cell_propor_box('near_airway','Neutrophil')->p
#Figure3B
make_cell_propor_box('immune','Neutrophil')->p
#Figure3C
SpatialFeaturePlot(lung_combine_after_batch[,which(lung_combine_after_batch$airway=='TRUE')],'s.Cxcl5.Cxcr2',images = c('slice1.3','slice1.4','slice1.7','slice1.8'),pt.size.factor = 2.2)
#Figure3D
SpatialFeaturePlot(lung_combine_after_batch[,which(lung_combine_after_batch$immune=='TRUE')],'r.Cxcl1.Cxcr2',images = c('slice1.3','slice1.4','slice1.7','slice1.8'),pt.size.factor = 2.2)
#Figure3E
SpatialFeaturePlot(lung_combine_after_batch[,which(lung_combine_after_batch$immune=='TRUE')],'r.TNF',images = c('slice1.3','slice1.4','slice1.7','slice1.8'),pt.size.factor = 2.2)


#Figure5A
make_cell_propor_box('near_airway','Interstitial-macrophage')->p
#Figure5B
make_cell_propor_box('immune','Interstitial-macrophage')->p
#Figure5C
SpatialFeaturePlot(lung_combine_after_batch[,which(lung_combine_after_batch$airway=='TRUE')],'s.Ccl5.Ccr5',images = c('slice1.5','slice1.6','slice1.7','slice1.8'),pt.size.factor = 2.2)
#Figure5D
SpatialFeaturePlot(lung_combine_after_batch[,which(lung_combine_after_batch$immune=='TRUE')],'r.Ccl2.Ccr2',images = c('slice1.5','slice1.6','slice1.7','slice1.8'),pt.size.factor = 2.2)
#Figure5E
SpatialFeaturePlot(lung_combine_after_batch[,which(lung_combine_after_batch$immune=='TRUE')],'r.VCAM',images = c('slice1.5','slice1.6','slice1.7','slice1.8'),pt.size.factor = 2.2)


#AMs
ggsave('p_AM_airway.png',make_cell_propor_box('near_airway','Alveolar-macrophage'),width = 8,height = 8)
ggsave('p_AM_distal.png',make_cell_propor_box('distal','Alveolar-macrophage'),width = 8,height = 8)

ggsave('p_NK_airway.png',make_cell_propor_box('near_airway','NK'),width = 8,height = 8)
ggsave('p_NK_distal.png',make_cell_propor_box('distal','NK'),width = 8,height = 8)

ggsave('p_DC_airway.png',make_cell_propor_box('near_airway','cDC'),width = 8,height = 8)
ggsave('p_DC_distal.png',make_cell_propor_box('distal','cDC'),width = 8,height = 8)

ggsave('p_Plasma_airway.png',make_cell_propor_box('near_airway','Plasma-cell'),width = 8,height = 8)
ggsave('p_Plasma_distal.png',make_cell_propor_box('distal','Plasma-cell'),width = 8,height = 8)

ggsave('p_Monocyte_airway.png',make_cell_propor_box('near_airway','Monocyte'),width = 8,height = 8)
ggsave('p_Monocyte_distal.png',make_cell_propor_box('distal','Monocyte'),width = 8,height = 8)

make_cell_propor_box('near_airway','NK')
make_cell_propor_box('immune','NK')
make_cell_propor_box('distal','NK')

make_cell_propor_box('near_airway','cDC')
make_cell_propor_box('immune','cDC')
make_cell_propor_box('distal','cDC')

make_cell_propor_box('near_airway','Monocyte')
make_cell_propor_box('immune','Monocyte')
make_cell_propor_box('distal','Monocyte')

make_cell_propor_box('near_airway','Cd4-tcell')
make_cell_propor_box('immune','Cd4-tcell')
make_cell_propor_box('distal','Cd4-tcell')

make_cell_propor_box('near_airway','B-cell')
make_cell_propor_box('immune','Cd4-tcell')
make_cell_propor_box('distal','Cd4-tcell')



#deconvolution cluster
scaled_mat <- t(scale(t(df)))
dist_mat <- dist(scaled_mat, method = "euclidean")
hc <- hclust(dist_mat, method = "ward.D2")  # or use cor_dist
plot(hc)
clusters <- cutree(hc, k = 29)  # Adjust number of clusters as needed

library(NMF)
set.seed(123)
nmf_res <- nmf(t(df[,-dim(df)[2]]), rank = 29, method = "brunet",nrun=30)
clusters <- predict(nmf_res)  # Assign spot to cluster with highest contribution
table(clusters)
pheatmap

kmeans(df[,-dim(df)[2]],10)->result
table(result$cluster)
pheatmap::pheatmap(df[order(result$cluster),],cluster_rows = FALSE,show_rownames = FALSE,scale='row',annotation_row = data.frame(row.names = rownames(df),cluster=as.character(result$cluster)))
lung_combine_after_batch$cluster<-result$cluster
DimPlot(lung_combine_after_batch,group.by = 'cluster')
pheatmap::pheatmap(df[order(lung_combine_after_batch$seurat_clusters),],cluster_rows = FALSE,show_rownames = FALSE,scale='row',annotation_row = data.frame(row.names = rownames(df),cluster=as.character(lung_combine_after_batch$seurat_clusters)))
pheatmap::pheatmap(df[order(lung_combine_after_batch$seurat_clusters),],show_rownames = FALSE,scale='row',annotation_row = data.frame(row.names = rownames(df),cluster=as.character(lung_combine_after_batch$seurat_clusters)))

inte<-intersect(rownames(lung_combine_after_batch@assays$Spatial@data),marker)
A <- as.matrix(lung_combine_after_batch@assays$Spatial@data[inte,])
B <- as.matrix(lung_combine_after_batch@assays$conv@counts[selected_cells,])
# Center the rows (subtract row means)
A_centered <- A - rowMeans(A)
B_centered <- B - rowMeans(B)

# Compute the norms (standard deviation for each row)
A_norms <- sqrt(rowSums(A_centered^2))
B_norms <- sqrt(rowSums(B_centered^2))

# Compute correlation using matrix multiplication
cor_matrix <- (A_centered %*% t(B_centered)) / (A_norms %*% t(B_norms))
cor_matrix[as.array(unlist(lapply(as.data.frame(cor_matrix), function(x) rownames(cor_matrix)[order(0-x)] ))),]->a
pheatmap::pheatmap(cor_matrix,scale = 'column')

#plot the cell type proportion in each cluster

#plot the figure 2
colors <- c(
  "#E41A1C",  # red
  "#377EB8",  # blue
  "#4DAF4A",  # green
  "#984EA3",  # purple
  "#FF7F00",  # orange
  "#A65628",  # brown
  "#F781BF",  # pink
  '#61ffe9'
)
names(colors)<-c(1,3,0,5,2,4,6,7)
index<-0:7
index<-data.frame(index,row.names = c(0,1,2,3,6,5,7,4))
lung_combine_after_batch_harmony$seurat_all<-lung_combine_after_batch_harmony$Spatial_snn_res.0.3
lung_combine_after_batch_harmony$seurat_all<-as.character(lung_combine_after_batch_harmony$seurat_all)
lung_combine_after_batch_harmony$seurat_all[which(lung_combine_after_batch_harmony$Spatial_snn_res.0.3==6)]<-4
lung_combine_after_batch_harmony$seurat_all[which(lung_combine_after_batch_harmony$Spatial_snn_res.0.3==7)]<-6
lung_combine_after_batch_harmony$seurat_all[which(lung_combine_after_batch_harmony$Spatial_snn_res.0.3==8)]<-7
lung_combine_after_batch_harmony$seurat_all<-as.array(index[lung_combine_after_batch_harmony$seurat_all,1])
SpatialDimPlot(lung_combine_after_batch_harmony,group.by = 'seurat_all',images = 'slice1.8',cols=colors)
lung_combine_after_batch_harmony$all_airway<-df$all_airway
SpatialFeaturePlot(lung_combine_after_batch_harmony,'all_airway',images = 'slice1.8')

SpatialDimPlot(lung_combine_after_batch_harmony,group.by = 'seurat_all',images = 'slice1.5',cols=colors,pt.size.factor = 2)->p
ggsave('p.png',p,width = 8,height = 8,dpi = 300)
DimPlot(lung_combine_after_batch_harmony,group.by = 'seurat_all',label = F,cols=colors)->p
ggsave('p1.png',p,width = 6,height = 6,dpi = 300)
SpatialDimPlot(lung_combine_after_batch_harmony,group.by = 'seurat_all',cols=colors,pt.size.factor = 2)->p
ggsave('p.png',p,width = 30,height = 15,dpi = 300)
