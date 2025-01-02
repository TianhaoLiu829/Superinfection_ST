
lung_combine_after_batch[,which(lung_combine_after_batch$airway==TRUE)]->airway
DefaultAssay(airway)<-'sender'
airway$macrophage<-airway[['conv']]@counts['Interstitial-macrophage',]
FindMarkers(airway,group.by = 'orig.ident',ident.1 = c('FS_1','FS_2'),ident.2 = c('F_1','F_2'),test.use = 'LR',latent.vars = c('immune_prop','macrophage'))->double_marker_send
FindMarkers(airway,group.by = 'orig.ident',ident.1 = c('FS_1','FS_2'),ident.2 = c('S_1','S_2'),test.use = 'LR',latent.vars = c('immune_prop','macrophage'))->double_marker_send_bac

#Ccl5-Ccr5
#Ccl7-Ccr2
SpatialFeaturePlot(airway,'s.Ccl7.Ccr2',pt.size.factor = 2.2)->p
ggsave('~/ST_Superinfection/s_Ccl7_Ccr2.png',p,width = 25,height = 18,dpi = 400)
SpatialFeaturePlot(airway,'s.Ccl5.Ccr5',pt.size.factor = 2.1)->p
ggsave('~/ST_Superinfection/s_Ccl5_Ccr5.png',p,width = 25,height = 14,dpi = 400)


SpatialFeaturePlot(airway,'s.Cxcl5.Cxcr2',pt.size.factor = 2.2)->p
ggsave('~/ST_Superinfection/s_Cxcl5_Cxcr2.png',p,width = 25,height = 18,dpi = 400)

SpatialFeaturePlot(airway,'r.MHC.I',pt.size.factor = 2.1)->p
ggsave('~/ST_Superinfection/r_MHC_I.png',p,width = 25,height = 18,dpi = 400)



DefaultAssay(airway)<-'reciever'
FindMarkers(airway,group.by = 'orig.ident',ident.1 = c('FS_1','FS_2'),ident.2 = c('F_1','F_2'))->double_marker_recieve
#MHC.I
#Ifng.ifngr1-ifngr2
SpatialFeaturePlot(airway,'s.Ifng.Ifngr1-Ifngr2',images = c('slice1.5','slice1.6','slice1.7','slice1.8'))
SpatialFeaturePlot(airway,'s.Ccl5.Ccr5',images = c('slice1.5','slice1.6','slice1.7','slice1.8'))








library(ggplot2)
library(ggpubr)
df<-as.data.frame(t(as.data.frame(lung_combine_after_batch[['conv']]@counts)))
n<-ncol(df)
lung_combine_after_batch$near_airway<-FALSE
lung_combine_after_batch$near_airway[which(lung_combine_after_batch$distance_min<200)]<-TRUE
lung_combine_after_batch$airway_auto<-FALSE
lung_combine_after_batch$airway_auto[which(lung_combine_after_batch$first_type%in%c('Ciliated_cell','Club_cell','Goblst_cell')|lung_combine_after_batch$second_type%in%c('Ciliated_cell','Club_cell','Goblst_cell'))]<-TRUE
target<-'ibalt'
#df<-as.data.frame(t(as.data.frame(lung_combine_after_batch[['sender']]@counts)))
for (object in 'r_total') {
  df->df2
  lung_combine_after_batch->lung_combine_after_batch2
  # set neutrophil as 0.0005
  if (object=='Neutrophil'){
    df2<-df2[which(df$Neutrophil>0.0005),]
    lung_combine_after_batch2<-lung_combine_after_batch2[,which(df$Neutrophil>0.0005)]
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
  print(p)
}
ggsave('~/ST_Superinfection/macrophage_inf_region.png',p,width = 10,height = 8,dpi = 300)

# Load libraries
library(ggplot2)
library(dplyr)



# Sample data
as.data.frame(lung_combine_after_batch@assays$reciever@counts+lung_combine_after_batch@assays$sender@counts)->df4
df4<-as.data.frame(t(df4))
df3 <- cbind(df4$r.Icosl.Icos,lung_combine_after_batch$orig.ident)
colnames(df3)<-c('value','label')
as.data.frame(df3)->df3
df3$value<-as.numeric(df3$value)
df3<-df3[which(lung_combine_after_batch$immune==TRUE),]
# Calculate count of non-zero values and log mean for each label
summary_df <- df3 %>%
  filter(value != 0) %>%
  group_by(label) %>%
  summarize(
    count_non_zero = n(),
    log_mean_value = mean(value)
  )
summary_df$count_non_zero<-as.numeric(summary_df$count_non_zero/table(df3$label)[match(summary_df$label,names(table(df3$label)))])
summary_df$label<-factor(summary_df$label,levels = c('N_1','N_2','S_1','S_2','F_1','F_2','FS_1','FS_2'))
# Plot
ggplot(summary_df[which(summary_df$label%in%c('F_1','F_2','FS_1','FS_2')),], aes(x = label, y = 1, size = count_non_zero, color = log_mean_value)) +
  geom_point() +
  scale_size_continuous(range = c(3, 15)) +  # Adjust size range as needed
  scale_color_viridis_c() +                  # Color scale for log mean
  labs(x = "Sample Label", y = "", size = "Non-Zero Count", color = "Log Mean Value") +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.grid=element_blank(),
    axis.text.x = element_text(size = 14),  # Adjust x-axis text size
    axis.text.y = element_text(size = 14),   # Adjust y-axis text size
    axis.title.x = element_text(size = 25),  # Adjust x-axis text size
    axis.title.y = element_text(size = 25)   # Adjust y-axis text size
  )+
  theme(axis.text.y = element_blank(),       # Remove y-axis text
        axis.ticks.y = element_blank(),      # Remove y-axis ticks
        legend.position = "right") +    # Remove y-axis title
  ggtitle("Dot Plot of Non-Zero Counts and Log Mean Value by Label")->p
ggsave('~/ST_Superinfection/ICOS_re_nar.png',p,width = 7,height = 2.7,dpi = 300)
SpatialFeaturePlot(lung_combine_after_batch,'r.ICOS',pt.size.factor = 2.1)->p
ggsave('~/ST_Superinfection/ICOS_re.png',p,width = 25,height = 17,dpi = 300)

# Sample data
as.data.frame(lung_combine_after_batch@assays$reciever@counts)->df4
df4<-as.data.frame(t(df4))
df3 <- cbind(df4$r.CD40,lung_combine_after_batch$orig.ident)
colnames(df3)<-c('value','label')
as.data.frame(df3)->df3
df3$value<-as.numeric(df3$value)
df3<-df3[which(lung_combine_after_batch$immune==TRUE),]
# Calculate count of non-zero values and log mean for each label
summary_df <- df3 %>%
  filter(value != 0) %>%
  group_by(label) %>%
  summarize(
    count_non_zero = n(),
    log_mean_value = median(value)
  )
summary_df$count_non_zero<-as.numeric(summary_df$count_non_zero/table(df3$label)[match(summary_df$label,names(table(df3$label)))])
summary_df$label<-factor(summary_df$label,levels = c('N_1','N_2','S_1','S_2','F_1','F_2','FS_1','FS_2'))
summarize(group_by(df3,label), log_mean_value=mean(value))$log_mean_value->summary_df$log_mean_value
aggregate(df3$value,by=list(df3$label),mean)$x[match(summary_df$label,aggregate(df3$value,by=list(df3$label),mean)$Group.1)]->summary_df$log_mean_value
# Plot
ggplot(summary_df[which(summary_df$label%in%c('F_1','F_2','FS_1','FS_2')),], aes(x = label, y = 1, size = count_non_zero, color = log_mean_value)) +
  geom_point() +
  scale_size_continuous(range = c(3, 15)) +  # Adjust size range as needed
  scale_color_viridis_c() +                  # Color scale for log mean
  labs(x = "Sample Label", y = "", size = "Non-Zero Count", color = "Log Mean Value") +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.grid=element_blank(),
    axis.text.x = element_text(size = 14),  # Adjust x-axis text size
    axis.text.y = element_text(size = 14),   # Adjust y-axis text size
    axis.title.x = element_text(size = 25),  # Adjust x-axis text size
    axis.title.y = element_text(size = 25)   # Adjust y-axis text size
  )+
  theme(axis.text.y = element_blank(),       # Remove y-axis text
        axis.ticks.y = element_blank(),      # Remove y-axis ticks
        axis.title.y = element_blank(),
        legend.position = "right") +    # Remove y-axis title
  ggtitle("Dot Plot of Non-Zero Counts and Log Mean Value by Label")->p
ggsave('~/ST_Superinfection/IFNI_se.png',p,width = 7,height = 2.7,dpi = 300)

# Sample data
as.data.frame(lung_combine_after_batch@assays$recie@counts)->df4
df4<-as.data.frame(t(df4))
df3 <- cbind(df4[,c('r.Icosl.Icos','r.Cd40lg.Cd40')],lung_combine_after_batch_archieve$orig.ident)
colnames(df3)<-c('r.Icosl.Icos','r.Cd40lg.Cd40','label')
as.data.frame(df3)->df3
df3<-df3[which(lung_combine_after_batch$immune),]
# Calculate count of non-zero values and log mean for each label
result <- df3 %>%
  pivot_longer(cols = starts_with("r"), names_to = "row", values_to = "value") %>%
  group_by(label, row) %>%
  summarise(
    proportion = mean(value != 0), # Proportion of non-zero values
    log_mean = log(mean(value[value != 0] + 1)) # Log mean value (avoid log(0))
  ) %>%
  ungroup()
# Bubble plot
result<-result[which(result$label%in%c('F_1','F_2','FS_1','FS_2')),]
ggplot(result[which(result$row=='r.Cd40lg.Cd40'),], aes(x = label, y = row, size = proportion, color = log_mean)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red", name = "Log Mean") +
  scale_size_continuous(name = "Proportion") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(x = "Group", y = "Features", title = "Bubble Plot")
ggsave('~/ST_Superinfection/CXCL13_CXCR5.png',p,width = 10,height = 4,dpi = 300)



for (object in 'CD8_score1') {
  df[which(lung_combine_after_batch@meta.data[,target]==TRUE&lung_combine_after_batch$orig.ident=='N_1'),object]->prop_macro_airway_N_1
  df[which(lung_combine_after_batch@meta.data[,target]==TRUE&lung_combine_after_batch$orig.ident=='N_2'),object]->prop_macro_airway_N_2
  df[which(lung_combine_after_batch@meta.data[,target]==TRUE&lung_combine_after_batch$orig.ident=='S_1'),object]->prop_macro_airway_S_1
  df[which(lung_combine_after_batch@meta.data[,target]==TRUE&lung_combine_after_batch$orig.ident=='S_2'),object]->prop_macro_airway_S_2
  df[which(lung_combine_after_batch@meta.data[,target]==TRUE&lung_combine_after_batch$orig.ident=='F_1'),object]->prop_macro_airway_F_1
  df[which(lung_combine_after_batch@meta.data[,target]==TRUE&lung_combine_after_batch$orig.ident=='F_2'),object]->prop_macro_airway_F_2
  df[which(lung_combine_after_batch@meta.data[,target]==TRUE&lung_combine_after_batch$orig.ident=='FS_1'),object]->prop_macro_airway_FS_1
  df[which(lung_combine_after_batch@meta.data[,target]==TRUE&lung_combine_after_batch$orig.ident=='FS_2'),object]->prop_macro_airway_FS_2
  data <- data.frame(
    value = c(prop_macro_airway_N_1,prop_macro_airway_N_2,prop_macro_airway_S_1,prop_macro_airway_S_2,prop_macro_airway_F_1, prop_macro_airway_F_2, prop_macro_airway_FS_1, prop_macro_airway_FS_2),
    group = factor(rep(c(paste0("prop_",object,"_airway_N_1"),paste0("prop_",object,"_airway_N_2"),paste0("prop_",object,"_airway_S_1"),paste0("prop_",object,"_airway_S_2"),paste0("prop_",object,"_airway_F_1"), paste0("prop_",object,"_airway_F_2"), paste0("prop_",object,"_airway_FS_1"),paste0("prop_",object,"_airway_FS_2") ), 
                       times = c(length(prop_macro_airway_N_1),length(prop_macro_airway_N_2),length(prop_macro_airway_S_1),length(prop_macro_airway_S_2),length(prop_macro_airway_F_1), length(prop_macro_airway_F_2), length(prop_macro_airway_FS_1), length(prop_macro_airway_FS_2))))
  )
  data$condition<-str_sub(data$group,-4)
  data$condition<-gsub('^_','',data$condition)
  data$condition2<-gsub('_.*','',data$condition)
  data$group<-factor(data$group,levels = c(paste0("prop_",object,"_airway_N_1"),paste0("prop_",object,"_airway_N_2"),paste0("prop_",object,"_airway_S_1"),paste0("prop_",object,"_airway_S_2"),paste0("prop_",object,"_airway_F_1"), paste0("prop_",object,"_airway_F_2"), paste0("prop_",object,"_airway_FS_1"),paste0("prop_",object,"_airway_FS_2")))
  ggplot(data, aes(x = group, y = value, fill = condition2)) +
    geom_violin(trim = FALSE, alpha = 0.4,scale = "width") +  # Add violin plot to show distribution, adjust transparency with alpha
    geom_boxplot(width = 0.3, position = position_dodge(0.9)) +  # Add boxplot, narrower width for visibility over violin
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
    labs(title = "Boxplot with Significance Annotations", x = "Group", y = object)->p
  print(p)
}


AddModuleScore(lung_combine_after_batch,list(c('Ccl5','Ms4a4b','Nkg7','Cd8a','Cd3d')),name = 'CD8_score')->lung_combine_after_batch
AddModuleScore(lung_combine_after_batch,features = list(c('Ccr7','Tcf7','Cd4','Il7r','Cd3d')),name = 'CD4_score')->lung_combine_after_batch
AddModuleScore(lung_combine_after_batch,features = list(c('Ccr7','Tcf7','Cd4','Il7r','Cd3d','Ccl5','Ms4a4b','Cxcr3','Ifng',"Tbx21","Gzmk")),name = 'Th1_score')->lung_combine_after_batch
AddModuleScore(lung_combine_after_batch,features = list(c('Ms4a1','Cd19','Ighm','Bank1','Pax5')),name = 'B_score')->lung_combine_after_batch
AddModuleScore(lung_combine_after_batch,features = list(c('Itgae','Xcr1','Clec9a','Itgax','Cd74','H2-Eb1','H2-Ab1','H2-Aa')),name = 'cDC_marker_1')->lung_combine_after_batch
AddModuleScore(lung_combine_after_batch,features = list(rownames(DC_marker)[1:6]),name = 'cDC_marker_2')->lung_combine_after_batch
AddModuleScore(lung_combine_after_batch,features = list(c('Csf3r','Ly6g','S100a9','S100a8','Retnlg')),name = 'neu_score')->lung_combine_after_batch

#ibalt centric analysis
lung_combine_after_batch$ibalt<-FALSE
hist(colSums(lung_combine_after_batch[['conv']]$counts[c("Cd4-tcell","Cd8-tcell","B-cell","cDC"),]))
hist(lung_combine_after_batch[['Spatial']]$data['Cxcl13',])
lung_combine_after_batch$ibalt[which(lung_combine_after_batch[['Spatial']]$data['Cxcl13',]>0.6&colSums(as.matrix(lung_combine_after_batch[['conv']]$counts[c("Cd4-tcell","Cd8-tcell","B-cell","cDC"),]))>0.1)]<-TRUE

table(lung_combine_after_batch$ibalt,lung_combine_after_batch$orig.ident)
ibalt<-lung_combine_after_batch[,which(lung_combine_after_batch$ibalt)]
Idents(ibalt)<-'orig.ident'
DefaultAssay(ibalt)<-'Spatial'
immune_name<-c('Alveolar-macrophage','Cd4-tcell','Cd8-tcell','cDC','B-cell','IgA-plasma-cell','Interstitial-macrophage','Macrophage','Mast-cell','mitotic-tcell','Monocyte','Neutrophil','NK','Plasma-cell')
as.data.frame(t(as.data.frame(ibalt[['conv']]$counts[immune_name,])))->immune_name
colnames(immune_name)<-gsub('-','.',colnames(immune_name))
ibalt@meta.data<-cbind(as.data.frame(ibalt@meta.data),immune_name)

#cxcl13_cxcr5
ibalt2<-ibalt[,c(which(ibalt$orig.ident=='F_1'),which(ibalt$orig.ident=='F_2'),which(ibalt$orig.ident=='FS_1'),which(ibalt$orig.ident=='FS_2'))]
pheatmap::pheatmap(as.data.frame(ibalt2@assays$reciever@counts['r.Cxcl13.Cxcr5',]),cluster_rows = FALSE,cluster_cols = FALSE,show_rownames = FALSE,show_colnames = FALSE)


FindMarkers(ibalt,ident.1 = c('FS_1','FS_2'),ident.2=c('F_1','F_2'),latent.vars = colnames(immune_name),test.use = 'LR')->Flu_ibalt_marker
DefaultAssay(ibalt)<-'reciever'
FindMarkers(ibalt,ident.1 = c('FS_1','FS_2'),ident.2=c('F_1','F_2'),latent.vars = colnames(immune_name),test.use = 'LR')->Flu_ibalt_marker_recie
DefaultAssay(ibalt)<-'sender'
FindMarkers(ibalt,ident.1 = c('FS_1','FS_2'),ident.2=c('F_1','F_2'),latent.vars = colnames(immune_name),test.use = 'LR')->Flu_ibalt_marker_send
Idents(lung_combine_after_batch)<-"ibalt"
SpatialDimPlot(lung_combine_after_batch,pt.size.factor = 2.1,group.by = 'seurat_clusters')+scale_fill_manual(values= c("lightblue", "orange"))->p
ggsave('~/ST_Superinfection/ibalt.png',p,width = 25,height = 17,dpi = 400)

#Th17 score 


lung_combine_after_batch$Th17_score1<-colMeans(as.data.frame(lung_combine_after_batch@assays$Spatial@data[c("Rorc","Il17a",'Ccr6','Il22','Il17f'),]))
lung_combine_after_batch$Th17_score1<-colMeans(as.data.frame(lung_combine_after_batch@assays$Spatial@data[c('Ccl5','Ms4a4b','Cxcr3','Ifng',"Tbx21","Gzmk"),]))

df$Th17_score<-lung_combine_after_batch$Th17_score1
df$Th17_score<-log((df$Th17_score+1),base = 2)

# Sample data
as.data.frame(lung_combine_after_batch$Th17_score1)->df4
df3 <- cbind(df4[,1],lung_combine_after_batch$orig.ident)
colnames(df3)<-c('value','label')
as.data.frame(df3)->df3
df3$value<-as.numeric(df3$value)
df3<-df3[which(lung_combine_after_batch$airway==TRUE),]
# Calculate count of non-zero values and log mean for each label
summary_df <- df3 %>%
  filter(value >0) %>%
  group_by(label) %>%
  summarize(
    count_non_zero = n(),
    log_mean_value = mean(value)
  )
summary_df$count_non_zero<-as.numeric(summary_df$count_non_zero/table(df3$label)[match(summary_df$label,names(table(df3$label)))])
summary_df$label<-factor(summary_df$label,levels = c('N_1','N_2','S_1','S_2','F_1','F_2','FS_1','FS_2'))
summary_df$log_mean_value<-summarize(group_by(df3,label),log=log(mean(value)))$log
# Plot
ggplot(summary_df, aes(x = label, y = 1, size = count_non_zero, color = log_mean_value)) +
  geom_point() +
  scale_size_continuous(range = c(3, 15)) +  # Adjust size range as needed
  scale_color_viridis_c() +                  # Color scale for log mean
  labs(x = "Sample Label", y = "", size = "Non-Zero Count", color = "Log Mean Value") +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.grid=element_blank(),
    axis.text.x = element_text(size = 14),  # Adjust x-axis text size
    axis.text.y = element_text(size = 14),   # Adjust y-axis text size
    axis.title.x = element_text(size = 25),  # Adjust x-axis text size
    axis.title.y = element_text(size = 25)   # Adjust y-axis text size
  )+
  theme(axis.text.y = element_blank(),       # Remove y-axis text
        axis.ticks.y = element_blank(),      # Remove y-axis ticks
        axis.title.y = element_blank(),
        legend.position = "right") +    # Remove y-axis title
  ggtitle("Dot Plot of Non-Zero Counts and Log Mean Value by Label")->p
ggsave('~/ST_Superinfection/Th17_pathway_airway.png',p,width = 10,height = 2.7,dpi = 300)









