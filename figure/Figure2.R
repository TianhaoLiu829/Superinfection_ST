library(Seurat)
#Figure 2A
DimPlot(lung_combine_before_batch,group.by = "orig.ident")
DimPlot(lung_combine_after_batch,group.by = "orig.ident")
#Figure 2B
FeaturePlot(lung_combine_before_batch,'Scgb3a2')
FeaturePlot(lung_combine_after_batch,'Scgb3a2')
#Figure 2C
c("#E5C494","#00CED1","#4DAF4A" ,"#984EA3","#FFB4A2" ,"#FFD700" ,"#A65628","#3357FF","#33FFF5","#FF33A1","#FF8333","#A133FF")->color
DimPlot(lung_combine_after_batch,group.by = 'seurat_clusters',cols = color)
#Figure 2D
df<-as.data.frame(lung_combine_after_batch@meta.data[,c('orig.ident','seurat_clusters')])
detach("package:tidyverse", unload = TRUE)
detach("package:dplyr", unload = TRUE)
detach("package:tidyverse", unload = TRUE)
library(dplyr)
df_proportion <- df %>%
  group_by(orig.ident, seurat_clusters) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))
library(ggplot2)
df_proportion$orig.ident<-factor(df_proportion$orig.ident,levels = c('N_1','N_2','S_1','S_2','F_1','F_2','FS_1','FS_2'))
ggplot(df_proportion, aes(x = orig.ident, y = proportion, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(
    x = "Label1",
    y = "Proportion",
    fill = "Label2",
    title = "Proportion of Label2 within each Label1"
  ) +
  theme_bw() +
  scale_fill_manual(values = color) +
  theme(legend.title = element_blank(),panel.grid=element_blank())
#Figure 2E
SpatialDimPlot(lung_combine_after_batch,group.by = 'seurat_clusters',label = FALSE,label.box=FALSE,images = 'slice1.7',pt.size.factor = 2.1)+scale_fill_manual(values=color)+ggplot2::theme(legend.key.size = unit(5,'mm'))
SpatialFeaturePlot(lung_combine_after_batch,'immune_prop',images = 'slice1.7',pt.size.factor = 2)
as.data.frame(lung_combine_after_batch[['conv']]@counts)->conv
conv_infla<-conv[c('Alveolar-macrophage','Interstitial-macrophage','Cd4-tcell','Cd8-tcell','B-cell','IgA-plasma-cell','Macrophage','cDC','Mast-cell','mitotic-tcell','Monocyte','Neutrophil','NK','Plasma-cell'),]
hist(colSums(conv_infla))
lung_combine_after_batch$immune_prop<-colSums(conv_infla)
SpatialFeaturePlot(lung_combine_after_batch,'immune_prop',images = 'slice1.7',pt.size.factor = 2)


