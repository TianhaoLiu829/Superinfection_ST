
make_cell_propor_box<-function(target,object){
  df->df2
  lung_combine_after_batch->lung_combine_after_batch
  # set spots with neutrophil proportion of 0.0005 as 0
  if (object=='Neutrophil'){
    df2<-df2[which(df$Neutrophil>0.0005),]
    lung_combine_after_batch<-lung_combine_after_batch[,which(df$Neutrophil>0.0005)]
  }
  df2[which(lung_combine_after_batch@meta.data[,target]==TRUE&lung_combine_after_batch$orig.ident=='N_1'),object]->prop_macro_airway_N_1
  df2[which(lung_combine_after_batch@meta.data[,target]==TRUE&lung_combine_after_batch$orig.ident=='N_2'),object]->prop_macro_airway_N_2
  df2[which(lung_combine_after_batch@meta.data[,target]==TRUE&lung_combine_after_batch$orig.ident=='S_1'),object]->prop_macro_airway_S_1
  df2[which(lung_combine_after_batch@meta.data[,target]==TRUE&lung_combine_after_batch$orig.ident=='S_2'),object]->prop_macro_airway_S_2
  df2[which(lung_combine_after_batch@meta.data[,target]==TRUE&lung_combine_after_batch$orig.ident=='F_1'),object]->prop_macro_airway_F_1
  df2[which(lung_combine_after_batch@meta.data[,target]==TRUE&lung_combine_after_batch$orig.ident=='F_2'),object]->prop_macro_airway_F_2
  df2[which(lung_combine_after_batch@meta.data[,target]==TRUE&lung_combine_after_batch$orig.ident=='FS_1'),object]->prop_macro_airway_FS_1
  df2[which(lung_combine_after_batch@meta.data[,target]==TRUE&lung_combine_after_batch$orig.ident=='FS_2'),object]->prop_macro_airway_FS_2
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
lung_combine_after_batch$near_airway<-FALSE
lung_combine_after_batch$near_airway[which(lung_combine_after_batch$distance_min<200)]<-TRUE


#find some differential LR pairs
lung_combine_after_batch[,which(lung_combine_after_batch$airway==TRUE)]->airway
DefaultAssay(airway)<-'sender'
airway$IM_prop<-airway@assays$conv@counts['Interstitial-macrophage',]
FindMarkers(airway,group.by = 'orig.ident',ident.1 = c('FS_1','FS_2'),ident.2 = c('F_1','F_2'),test.use = 'LR',latent.vars = c('immune_prop','IM_prop'))->double_marker_send_flu
FindMarkers(airway,group.by = 'orig.ident',ident.1 = c('FS_1','FS_2'),ident.2 = c('S_1','S_2'),test.use = 'LR',latent.vars = c('immune_prop','IM_prop'))->double_marker_send_bac
DefaultAssay(airway)<-'reciever'
FindMarkers(airway,group.by = 'orig.ident',ident.1 = c('FS_1','FS_2'),ident.2 = c('F_1','F_2'),test.use = 'LR',latent.vars = c('immune_prop','IM_prop'))->double_marker_recei_flu
FindMarkers(airway,group.by = 'orig.ident',ident.1 = c('FS_1','FS_2'),ident.2 = c('S_1','S_2'),test.use = 'LR',latent.vars = c('immune_prop','IM_prop'))->double_marker_recei_bac


#Figure3A
make_cell_propor_box('near_airway','Neutrophil')
#Figure3B
make_cell_propor_box('immune','Neutrophil')
#Figure3C
SpatialFeaturePlot(lung_combine_after_batch[,which(lung_combine_after_batch$airway=='TRUE')],'s.Cxcl5.Cxcr2',images = c('slice1.3','slice1.4','slice1.7','slice1.8'),pt.size.factor = 2.2)
#Figure3D
SpatialFeaturePlot(lung_combine_after_batch[,which(lung_combine_after_batch$immune=='TRUE')],'r.Cxcl1.Cxcr2',images = c('slice1.3','slice1.4','slice1.7','slice1.8'),pt.size.factor = 2.2)
#Figure3E
SpatialFeaturePlot(lung_combine_after_batch[,which(lung_combine_after_batch$immune=='TRUE')],'r.TNF',images = c('slice1.3','slice1.4','slice1.7','slice1.8'),pt.size.factor = 2.2)


#Figure4A
make_cell_propor_box('near_airway','Interstitial-macrophage')
#Figure4B
make_cell_propor_box('immune','Interstitial-macrophage')
#Figure4C
SpatialFeaturePlot(lung_combine_after_batch[,which(lung_combine_after_batch$airway=='TRUE')],'s.Ccl5.Ccr5',images = c('slice1.5','slice1.6','slice1.7','slice1.8'),pt.size.factor = 2.2)
#Figure4D
SpatialFeaturePlot(lung_combine_after_batch[,which(lung_combine_after_batch$immune=='TRUE')],'r.Ccl2.Ccr2',images = c('slice1.5','slice1.6','slice1.7','slice1.8'),pt.size.factor = 2.2)
#Figure4E
SpatialFeaturePlot(lung_combine_after_batch[,which(lung_combine_after_batch$immune=='TRUE')],'r.VCAM',images = c('slice1.5','slice1.6','slice1.7','slice1.8'),pt.size.factor = 2.2)





