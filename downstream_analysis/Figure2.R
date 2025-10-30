#figure 2

DimPlot(lung_combine_after_batch,group.by = "orig.ident")
DimPlot(lung_combine_after_batch,group.by = "seurat_clusters")


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
