library(ggplot2)
library(ggpubr)
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
