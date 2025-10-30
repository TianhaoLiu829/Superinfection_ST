c("#E5C494","#00CED1","#4DAF4A" ,"#984EA3","#FFB4A2" ,"#FFD700" ,"#A65628","#3357FF","#33FFF5","#FF33A1","#FF8333","#A133FF")->color
c("grey","#FFD700")->color2
c("grey","#00CED1")->color2
SpatialDimPlot(lung_combine_after_batch,group.by = 'immune',label = FALSE,label.box=FALSE,images = 'slice1',pt.size.factor = 2.1)+scale_fill_manual(values=color2)+ggplot2::theme(legend.key.size = unit(5,'mm'))->p1
SpatialDimPlot(lung_combine_after_batch,group.by = 'immune',label = FALSE,label.box=FALSE,images = 'slice1.2',pt.size.factor = 2.1)+scale_fill_manual(values=color2)+ggplot2::theme(legend.key.size = unit(5,'mm'))->p2
SpatialDimPlot(lung_combine_after_batch,group.by = 'immune',label = FALSE,label.box=FALSE,images = 'slice1.3',pt.size.factor = 2.1)+scale_fill_manual(values=color2)+ggplot2::theme(legend.key.size = unit(5,'mm'))->p3
SpatialDimPlot(lung_combine_after_batch,group.by = 'immune',label = FALSE,label.box=FALSE,images = 'slice1.4',pt.size.factor = 2.1)+scale_fill_manual(values=color2)+ggplot2::theme(legend.key.size = unit(5,'mm'))->p4
SpatialDimPlot(lung_combine_after_batch,group.by = 'immune',label = FALSE,label.box=FALSE,images = 'slice1.5',pt.size.factor = 2.1)+scale_fill_manual(values=color2)+ggplot2::theme(legend.key.size = unit(5,'mm'))->p5
SpatialDimPlot(lung_combine_after_batch,group.by = 'immune',label = FALSE,label.box=FALSE,images = 'slice1.6',pt.size.factor = 2.1)+scale_fill_manual(values=color2)+ggplot2::theme(legend.key.size = unit(5,'mm'))->p6
SpatialDimPlot(lung_combine_after_batch,group.by = 'immune',label = FALSE,label.box=FALSE,images = 'slice1.7',pt.size.factor = 2.1)+scale_fill_manual(values=color2)+ggplot2::theme(legend.key.size = unit(5,'mm'))->p7
SpatialDimPlot(lung_combine_after_batch,group.by = 'immune',label = FALSE,label.box=FALSE,images = 'slice1.8',pt.size.factor = 1.8)+scale_fill_manual(values=color2)+ggplot2::theme(legend.key.size = unit(5,'mm'))->p8



lung_combine_after_batch_archieve$seurat_clusters2<-factor(lung_combine_after_batch_archieve$seurat_clusters2,levels = seq(1,12,1))
SpatialDimPlot(lung_combine_after_batch_archieve,group.by = 'seurat_clusters2',label = FALSE,label.box=FALSE,images = 'slice1',pt.size.factor = 2.1)+scale_fill_manual(values=color)+ggplot2::theme(legend.key.size = unit(5,'mm'))->p1
SpatialDimPlot(lung_combine_after_batch_archieve,group.by = 'seurat_clusters2',label = FALSE,label.box=FALSE,images = 'slice1.2',pt.size.factor = 2.1)+scale_fill_manual(values=color)+ggplot2::theme(legend.key.size = unit(5,'mm'))->p2
SpatialDimPlot(lung_combine_after_batch_archieve,group.by = 'seurat_clusters2',label = FALSE,label.box=FALSE,images = 'slice1.3',pt.size.factor = 2.1)+scale_fill_manual(values=color)+ggplot2::theme(legend.key.size = unit(5,'mm'))->p3
SpatialDimPlot(lung_combine_after_batch_archieve,group.by = 'seurat_clusters2',label = FALSE,label.box=FALSE,images = 'slice1.4',pt.size.factor = 2.1)+scale_fill_manual(values=color)+ggplot2::theme(legend.key.size = unit(5,'mm'))->p4
SpatialDimPlot(lung_combine_after_batch_archieve,group.by = 'seurat_clusters2',label = FALSE,label.box=FALSE,images = 'slice1.5',pt.size.factor = 2.1)+scale_fill_manual(values=color)+ggplot2::theme(legend.key.size = unit(5,'mm'))->p5
SpatialDimPlot(lung_combine_after_batch_archieve,group.by = 'seurat_clusters2',label = FALSE,label.box=FALSE,images = 'slice1.6',pt.size.factor = 2.1)+scale_fill_manual(values=color)+ggplot2::theme(legend.key.size = unit(5,'mm'))->p6
SpatialDimPlot(lung_combine_after_batch_archieve,group.by = 'seurat_clusters2',label = FALSE,label.box=FALSE,images = 'slice1.7',pt.size.factor = 2.1)+scale_fill_manual(values=color)+ggplot2::theme(legend.key.size = unit(5,'mm'))->p7
SpatialDimPlot(lung_combine_after_batch_archieve,group.by = 'seurat_clusters2',label = FALSE,label.box=FALSE,images = 'slice1.8',pt.size.factor = 1.7)+scale_fill_manual(values=color)+ggplot2::theme(legend.key.size = unit(5,'mm'))->p8
ggsave('~/ST_Superinfection/Figure1/p1.png',p1,width = 10,height = 9,dpi = 400)
ggsave('~/ST_Superinfection/Figure1/p2.png',p2,width = 10,height = 9,dpi = 400)
ggsave('~/ST_Superinfection/Figure1/p3.png',p3,width = 10,height = 9,dpi = 400)
ggsave('~/ST_Superinfection/Figure1/p4.png',p4,width = 10,height = 9,dpi = 400)
ggsave('~/ST_Superinfection/Figure1/p5.png',p5,width = 10,height = 9,dpi = 400)
ggsave('~/ST_Superinfection/Figure1/p6.png',p6,width = 10,height = 9,dpi = 400)
ggsave('~/ST_Superinfection/Figure1/p7.png',p7,width = 10,height = 9,dpi = 400)
ggsave('~/ST_Superinfection/Figure1/p8.png',p8,width = 10,height = 9,dpi = 400)

DimPlot(lung_combine_after_batch_archieve,group.by = 'seurat_clusters2',cols = color)->p_all
ggsave('~/ST_Superinfection/Figure1/p_all.png',p_all,width = 7,height = 6,dpi = 400)


read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/send4_N_1.csv')->send_N_1
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/send4_N_2.csv')->send_N_2
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/send4_S_1.csv')->send_S_1
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/send4_S_2.csv')->send_S_2
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/send4_F_1.csv')->send_F_1
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/send4_F_2.csv')->send_F_2
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/send4_FS_1.csv')->send_FS_1
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/send4_FS_2.csv')->send_FS_2

send_N_1$X<-paste0('N_1_',send_N_1$X)
send_N_2$X<-paste0('N_2_',send_N_2$X)
send_S_1$X<-paste0('S_1_',send_S_1$X)
send_S_2$X<-paste0('S_2_',send_S_2$X)
send_F_1$X<-paste0('F_1_',send_F_1$X)
send_F_2$X<-paste0('F_2_',send_F_2$X)
send_FS_1$X<-paste0('FS_1_',send_FS_1$X)
send_FS_2$X<-paste0('FS_2_',send_FS_2$X)


rbind(send_N_1,send_N_2,send_S_1,send_S_2,send_F_1,send_F_2,send_FS_1,send_FS_2)->send


read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/recie4_N_1.csv')->recie_N_1
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/recie4_N_2.csv')->recie_N_2
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/recie4_S_1.csv')->recie_S_1
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/recie4_S_2.csv')->recie_S_2
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/recie4_F_1.csv')->recie_F_1
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/recie4_F_2.csv')->recie_F_2
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/recie4_FS_1.csv')->recie_FS_1
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/recie4_FS_2.csv')->recie_FS_2

recie_N_1$X<-paste0('N_1_',recie_N_1$X)
recie_N_2$X<-paste0('N_2_',recie_N_2$X)
recie_S_1$X<-paste0('S_1_',recie_S_1$X)
recie_S_2$X<-paste0('S_2_',recie_S_2$X)
recie_F_1$X<-paste0('F_1_',recie_F_1$X)
recie_F_2$X<-paste0('F_2_',recie_F_2$X)
recie_FS_1$X<-paste0('FS_1_',recie_FS_1$X)
recie_FS_2$X<-paste0('FS_2_',recie_FS_2$X)

rbind(recie_N_1,recie_N_2,recie_S_1,recie_S_2,recie_F_1,recie_F_2,recie_FS_1,recie_FS_2)->recie

rownames(recie)<-recie$X
recie<-recie[,-1]
rownames(send)<-send$X
send<-send[,-1]

inter<-intersect(rownames(recie),colnames(lung_combine_after_batch))
lung_combine_after_batch<-lung_combine_after_batch[,inter]
recie<-recie[inter,]
send<-send[inter,]

lung_combine_after_batch[['reciever']]<-CreateAssayObject(as.data.frame(t(recie[match(colnames(lung_combine_after_batch),rownames(recie)),])))
lung_combine_after_batch[['sender']]<-CreateAssayObject(as.data.frame(t(send[match(colnames(lung_combine_after_batch),rownames(send)),])))

library(parallel)
library(Seurat)
# Define a function to find markers for each cluster
find_markers_for_cluster <- function(cluster_id) {
  FindMarkers(lung_combine_after_batch, group.by = 'seurat_clusters',ident.1 = cluster_id,only.pos = TRUE)
}
cluster_ids <- levels(lung_combine_after_batch$seurat_clusters)
# Use future_lapply to parallelize
markers_list <- mclapply(cluster_ids, find_markers_for_cluster,mc.cores = 14)

# Combine results
markers <- do.call(rbind, markers_list)

findall_sc<-function(object,test){
  alllist<-data.frame()
  DefaultAssay(object)<-test
  for (i in 1:length(unique(object$orig.ident))) {
    FindMarkers(object,group.by = 'orig.ident',ident.1 = unique(object$orig.ident)[i],only.pos = TRUE)->markerlist
    markerlist$cluster<-unique(object$orig.ident)[i]
    rownames(markerlist)->markerlist$gene
    rbind(alllist,markerlist)->alllist
  }
  return(alllist)
}
findall_sc(lung_combine_after_batch,test = 'sender')->marker_sender
findall_sc(lung_combine_after_batch,test = 'reciever')->marker_reciever
findall_sc(lung_combine_after_batch[,which(lung_combine_after_batch$seurat_clusters==1)],test = 'sender')->marker_sender_neu
findall_sc(lung_combine_after_batch[,which(lung_combine_after_batch$seurat_clusters==1)],test = 'reciever')->marker_reciever_neu
findall_sc(lung_combine_after_batch[,which(lung_combine_after_batch$seurat_clusters==3)],test = 'sender')->marker_sender_macro
findall_sc(lung_combine_after_batch[,which(lung_combine_after_batch$seurat_clusters==3)],test = 'reciever')->marker_reciever_macro


lung_combine_after_batch[['conv']]['Neutrophil',]->lung_combine_after_batch$neutrophil
lung_combine_after_batch[['conv']]['Interstitial-macrophage',]->lung_combine_after_batch$IM
DefaultAssay(lung_combine_after_batch)<-'sender'
FindMarkers(lung_combine_after_batch[,which(lung_combine_after_batch$seurat_clusters==1)],group.by='orig.ident', ident.1 = c('S_1_','S_2_'),ident.2 = c('FS_1','FS_2'),latent.vars = 'neutrophil',test.use = 'LR')->marker_sender_neu_fvsfs
DefaultAssay(lung_combine_after_batch)<-'reciever'
FindMarkers(lung_combine_after_batch[,which(lung_combine_after_batch$seurat_clusters==1)],group.by='orig.ident', ident.1 = c('S_1_','S_2_'),ident.2 = c('FS_1','FS_2'),latent.vars = 'neutrophil',test.use = 'LR')->marker_reciever_neu_fvsfs

DefaultAssay(lung_combine_after_batch)<-'sender'
FindMarkers(lung_combine_after_batch[,which(lung_combine_after_batch$seurat_clusters==3)],group.by='orig.ident', ident.1 = c('F_1_','F_2_'),ident.2 = c('FS_1','FS_2'),latent.vars = 'IM',test.use = 'LR')->marker_sender_macro_fvsfs
DefaultAssay(lung_combine_after_batch)<-'reciever'
FindMarkers(lung_combine_after_batch[,which(lung_combine_after_batch$seurat_clusters==3)],group.by='orig.ident', ident.1 = c('F_1_','F_2_'),ident.2 = c('FS_1','FS_2'),latent.vars = 'IM',test.use = 'LR')->marker_reciever_macro_fvsfs


# Function to perform regression on a pair of columns
# Initialize matrices for beta values and p-values
df<-as.data.frame(t(as.data.frame(lung_combine_after_batch[['conv']]@counts)))
n<-ncol(df)

# Loop through each pair of columns
lm_test<-function(df){
  n<-ncol(df)
  beta_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(colnames(df), colnames(df)))
  pval_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(colnames(df), colnames(df)))
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Perform linear regression
      model <- lm(df[, j] ~ df[, i])
      summary_model <- summary(model)
      
      # Extract beta (slope) and p-value
      beta <- summary_model$coefficients[2, "Estimate"]
      p_value <- summary_model$coefficients[2, "Pr(>|t|)"]
      
      # Fill the matrices
      beta_matrix[i, j] <- beta
      beta_matrix[j, i] <- beta  # symmetric matrix
      
      pval_matrix[i, j] <- p_value
      pval_matrix[j, i] <- p_value  # symmetric matrix
    }
  }
  return(list(beta_matrix,pval_matrix))
}
lm_test(df[which(lung_combine_after_batch$orig.ident=='N_1_'),])->N_1_cor
lm_test(df[which(lung_combine_after_batch$orig.ident=='N_2_'),])->N_2_cor
lm_test(df[which(lung_combine_after_batch$orig.ident=='S_1_'),])->S_1_cor
lm_test(df[which(lung_combine_after_batch$orig.ident=='S_2_'),])->S_2_cor
lm_test(df[which(lung_combine_after_batch$orig.ident=='F_1_'),])->F_1_cor
lm_test(df[which(lung_combine_after_batch$orig.ident=='F_2_'),])->F_2_cor
lm_test(df[which(lung_combine_after_batch$orig.ident=='FS_1'),])->FS_1_cor
lm_test(df[which(lung_combine_after_batch$orig.ident=='FS_2'),])->FS_2_cor

N_1_cor[[1]]['B-cell','Cd4-tcell']
N_2_cor[[1]]['B-cell','Cd4-tcell']
S_1_cor[[1]]['B-cell','Cd4-tcell']
S_2_cor[[1]]['B-cell','Cd4-tcell']
F_1_cor[[1]]['B-cell','Cd4-tcell']
F_2_cor[[1]]['B-cell','Cd4-tcell']
FS_1_cor[[1]]['B-cell','Cd4-tcell']
FS_2_cor[[1]]['B-cell','Cd4-tcell']

N_1_cor[[2]]['B-cell','Cd4-tcell']
N_2_cor[[2]]['B-cell','Cd4-tcell']
S_1_cor[[2]]['B-cell','Cd4-tcell']
S_2_cor[[2]]['B-cell','Cd4-tcell']
F_1_cor[[2]]['B-cell','Cd4-tcell']
F_2_cor[[2]]['B-cell','Cd4-tcell']
FS_1_cor[[2]]['B-cell','Cd4-tcell']
FS_2_cor[[2]]['B-cell','Cd4-tcell']


F_1_cor[[2]]['Interstitial-macrophage','Cd4-tcell']
F_2_cor[[2]]['Interstitial-macrophage','Cd4-tcell']
FS_1_cor[[2]]['Interstitial-macrophage','Cd4-tcell']
FS_2_cor[[2]]['Interstitial-macrophage','Cd4-tcell']
F_1_cor[[1]]['Interstitial-macrophage','Cd4-tcell']
F_2_cor[[1]]['Interstitial-macrophage','Cd4-tcell']
FS_1_cor[[1]]['Interstitial-macrophage','Cd4-tcell']
FS_2_cor[[1]]['Interstitial-macrophage','Cd4-tcell']


F_1_cor[[2]]['Interstitial-macrophage','Ciliated-cell']
F_2_cor[[2]]['Interstitial-macrophage','Ciliated-cell']
FS_1_cor[[2]]['Interstitial-macrophage','Ciliated-cell']
FS_2_cor[[2]]['Interstitial-macrophage','Ciliated-cell']
F_1_cor[[1]]['Interstitial-macrophage','Ciliated-cell']
F_2_cor[[1]]['Interstitial-macrophage','Ciliated-cell']
FS_1_cor[[1]]['Interstitial-macrophage','Ciliated-cell']
FS_2_cor[[1]]['Interstitial-macrophage','Ciliated-cell']
library(ggpubr)
library(ggplot2)
library(stringr)
for (object in rownames(lung_combine_after_batch@assays$conv)) {
  df[which(lung_combine_after_batch$first_type%in%c('Ciliated_cell','Goblst_cell','Club_cell')&lung_combine_after_batch$orig.ident=='N_1_'),object]->prop_macro_airway_N_1
  df[which(lung_combine_after_batch$first_type%in%c('Ciliated_cell','Goblst_cell','Club_cell')&lung_combine_after_batch$orig.ident=='N_2_'),object]->prop_macro_airway_N_2
  df[which(lung_combine_after_batch$first_type%in%c('Ciliated_cell','Goblst_cell','Club_cell')&lung_combine_after_batch$orig.ident=='S_1_'),object]->prop_macro_airway_S_1
  df[which(lung_combine_after_batch$first_type%in%c('Ciliated_cell','Goblst_cell','Club_cell')&lung_combine_after_batch$orig.ident=='S_2_'),object]->prop_macro_airway_S_2
  df[which(lung_combine_after_batch$first_type%in%c('Ciliated_cell','Goblst_cell','Club_cell')&lung_combine_after_batch$orig.ident=='F_1_'),object]->prop_macro_airway_F_1
  df[which(lung_combine_after_batch$first_type%in%c('Ciliated_cell','Goblst_cell','Club_cell')&lung_combine_after_batch$orig.ident=='F_2_'),object]->prop_macro_airway_F_2
  df[which(lung_combine_after_batch$first_type%in%c('Ciliated_cell','Goblst_cell','Club_cell')&lung_combine_after_batch$orig.ident=='FS_1'),object]->prop_macro_airway_FS_1
  df[which(lung_combine_after_batch$first_type%in%c('Ciliated_cell','Goblst_cell','Club_cell')&lung_combine_after_batch$orig.ident=='FS_2'),object]->prop_macro_airway_FS_2
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
    theme_minimal() +
    scale_x_discrete(labels =unique(data$condition))+
    theme(
      axis.line = element_line(), # Add axis lines
      axis.text.x = element_text(size = 18),  # Adjust x-axis text size
      axis.text.y = element_text(size = 18),  # Adjust y-axis text size
      axis.title.x = element_text(size = 20), # Adjust x-axis title size
      axis.title.y = element_text(size = 20),  # Adjust y-axis title size
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()  
    )+
    labs(title = "Boxplot with Significance Annotations", x = "Group", y = "Value")->p
  print(p)
  
}


# Combine them into a data frame
data <- data.frame(
  value = c(prop_macro_airway_F_1, prop_macro_airway_F_2, prop_macro_airway_FS_1, prop_macro_airway_FS_2),
  group = factor(rep(c("prop_macro_airway_F_1", "prop_macro_airway_F_2", "prop_macro_airway_FS_1", "prop_macro_airway_FS_2"), 
                     times = c(length(prop_macro_airway_F_1), length(prop_macro_airway_F_2), length(prop_macro_airway_FS_1), length(prop_macro_airway_FS_2))))
)
library(ggplot2)
library(ggpubr)
ggboxplot(data, x = "group", y = "value") +
  stat_compare_means(comparisons = list(c("prop_macro_airway_F_1", "prop_macro_airway_F_2"), c("prop_macro_airway_F_1", "prop_macro_airway_FS_1"), c("prop_macro_airway_F_1", "prop_macro_airway_FS_2"),
                                        c("prop_macro_airway_F_2", "prop_macro_airway_FS_1"), c("prop_macro_airway_F_2", "prop_macro_airway_FS_2"), c("prop_macro_airway_FS_1", "prop_macro_airway_FS_2")), 
                     method = "wilcox.test") +
  theme_minimal() +
  labs(title = "Boxplot with Significance Annotations", x = "Group", y = "Value")

# Print the matrices
print("Beta Values Matrix:")
print(beta_matrix)

print("P-Values Matrix:")
print(pval_matrix)
# Apply the regression_test function to all pairs of columns
apply(combn(ncol(df), 2), 2, regression_test)

df$B_zone<-'B_cell_absent'
df$B_zone[which(df[,'B-cell']>0.02)]<-'B_cell_enriched'
t.test(df[which(df$B_zone=='B_cell_enriched'&lung_combine_after_batch$orig.ident=='F_1_'),'Cd4-tcell'],df[which(df$B_zone=='B_cell_absent'&lung_combine_after_batch$orig.ident=='F_1_'),'Cd4-tcell'])
t.test(df[which(df$B_zone=='B_cell_enriched'&lung_combine_after_batch$orig.ident=='F_2_'),'Cd4-tcell'],df[which(df$B_zone=='B_cell_absent'&lung_combine_after_batch$orig.ident=='F_2_'),'Cd4-tcell'])
t.test(df[which(df$B_zone=='B_cell_enriched'&lung_combine_after_batch$orig.ident=='FS_2'),'Cd4-tcell'],df[which(df$B_zone=='B_cell_absent'&lung_combine_after_batch$orig.ident=='FS_2'),'Cd4-tcell'])
t.test(df[which(df$B_zone=='B_cell_enriched'&lung_combine_after_batch$orig.ident=='FS_2'),'Cd4-tcell'],df[which(df$B_zone=='B_cell_absent'&lung_combine_after_batch$orig.ident=='FS_2'),'Cd4-tcell'])

# Load the necessary library
library(ggplot2)

# Sample data frame
df <- data.frame(
  Value = rnorm(100),
  Group1 = rep(c("A", "B"), each = 50),
  Group2 = rep(c("X", "Y"), times = 50)
)

# Create the boxplot
df$orig<-lung_combine_after_batch$orig.ident[match(rownames(df),colnames(lung_combine_after_batch))]
ggplot(df[,c('Cd4-tcell','B_zone','orig')], aes(x = interaction(B_zone,orig), y = `Cd4-tcell`)) +
  geom_boxplot() +
  labs(x = "Group1 and Group2", y = "Value") +
  theme_minimal()


as.data.frame(lung_combine_after_batch[['conv']]@counts)->conv
conv_infla<-conv[c('Alveolar-macrophage','Interstitial-macrophage','Cd4-tcell','Cd8-tcell','B-cell','IgA-plasma-cell','Macrophage','cDC','Mast-cell','mitotic-tcell','Monocyte','Neutrophil','NK','Plasma-cell'),]
hist(colSums(conv_infla))
lung_combine_after_batch$immune_prop<-colSums(conv_infla)
lung_combine_after_batch$immune<-FALSE
lung_combine_after_batch$immune[which(lung_combine_after_batch$immune_prop>0.3)]<-TRUE
SpatialDimPlot(lung_combine_after_batch,group.by = 'immune',images = 'slice1.6')
table(lung_combine_after_batch$immune)
conv_infla2<-conv_infla[,which(lung_combine_after_batch$immune)]
conv_infla2 <- data.frame(t(conv_infla2), label = lung_combine_after_batch$orig.ident[which(lung_combine_after_batch$immune)])
colnames(conv_infla2)[dim(conv_infla2)[2]]<-'Sample'

# Reshape data to long format
library(dplyr)
library(tidyverse)
conv_infla2 <- gather(conv_infla2, key = "cell_type", value = "proportion", -Sample)
conv_infla2$proportion <- as.numeric(conv_infla2$proportion )  # Ensure values are numeric
unique_rows <- unique(conv_infla2$cell_type)
conv_infla2$Sample<-gsub('_$','',conv_infla2$Sample)
conv_infla2$Sample<-factor(conv_infla2$Sample,levels = c('N_1','N_2','S_1','S_2','F_1','F_2','FS_1','FS_2'))

conv_infla2$condition<-substr(conv_infla2$Sample,1,(nchar(as.character(conv_infla2$Sample))-2))
for (r in unique_rows) {
  ggplot(subset(conv_infla2, cell_type == r), aes(x = Sample, y = proportion, fill = condition)) +
    geom_boxplot() +
    labs(title = paste("Boxplot for", r), x = "Sample", y = r) +
    theme_minimal() +
    theme(      axis.line = element_line(), # Add axis lines
                axis.text.x = element_text(size = 18),  # Adjust x-axis text size
                axis.text.y = element_text(size = 18),  # Adjust y-axis text size
                axis.title.x = element_text(size = 20), # Adjust x-axis title size
                axis.title.y = element_text(size = 20),  # Adjust y-axis title size
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank()  
    )+
    theme(legend.position = "none") -> p
  
  print(p)  # This will display the plot in the loop
}

#CD8 we can show
#Interstitial immune regions can be retained
#Also talk about the adaptive immunity (T cell and B cell)
#Indeed, IL-17R−/− and IL-22−/− mice displayed impaired bacterial clearance of S. aureus compared with that of wild-type mice. Mice challenged with influenza A PR/8/34 H1N1 and subsequently with S. aureus had increased inflammation and decreased clearance of both virus and bacteria. Coinfection resulted in greater type I and II IFN production in the lung compared with that with virus infection alone. Importantly, influenza A coinfection resulted in substantially decreased IL-17, IL-22, and IL-23 production after S. aureus infection. The decrease in S. aureus-induced IL-17, IL-22, and IL-23 was independent of type II IFN but required type I IFN production in influenza A-infected mice. Furthermore, overexpression of IL-23 in influenza A, S. aureus-coinfected mice rescued the induction of IL-17 and IL-22 and markedly improved bacterial clearance. These data indicate a novel mechanism by which influenza A-induced type I IFNs inhibit Th17 immunity and increase susceptibility to secondary bacterial pneumonia.
#test the association between those two types of cells





