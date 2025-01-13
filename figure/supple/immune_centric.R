#

DefaultAssay(lung_combine_after_batch)<-'sender'
Idents(lung_combine_after_batch)<-'orig.ident'
FindMarkers(lung_combine_after_batch[,which(lung_combine_after_batch$immune_exclude_airway==TRUE)],ident.1 = c('FS_1','FS_2'),ident.2 = c('F_1','F_2'))->marker_sender_macro
FindMarkers(lung_combine_after_batch[,which(lung_combine_after_batch$immune_exclude_airway==TRUE)],ident.1 = c('FS_1','FS_2'),ident.2 = c('S_1','S_2'))->marker_sender_neu

DefaultAssay(lung_combine_after_batch)<-'reciever'
FindMarkers(lung_combine_after_batch[,which(lung_combine_after_batch$immune_exclude_airway==TRUE)],ident.1 = c('FS_1','FS_2'),ident.2 = c('F_1','F_2'))->marker_recei_macro
FindMarkers(lung_combine_after_batch[,which(lung_combine_after_batch$immune_exclude_airway==TRUE)],ident.1 = c('FS_1','FS_2'),ident.2 = c('S_1','S_2'))->marker_recei_neu



#define immune regions
conv<-as.data.frame(lung_combine_after_batch[['conv']]$counts)
conv_infla<-conv[c('Alveolar-macrophage','Cd4-tcell','Cd8-tcell','cDC','B-cell','IgA-plasma-cell','Interstitial-macrophage','Macrophage','Mast-cell','mitotic-tcell','Monocyte','Neutrophil','NK','Plasma-cell'),]
lung_combine_after_batch$immune_prop<-colSums(conv_infla)[match(colnames(lung_combine_after_batch),colnames(conv_infla))]
lung_combine_after_batch$immune<-FALSE
lung_combine_after_batch$immune[which(lung_combine_after_batch$immune_prop>0.35)]<-TRUE
lung_combine_after_batch$ALL<-TRUE
SpatialDimPlot(lung_combine_after_batch,group.by = 'immune')->p
ggsave('p.png',p,width = 35,height = 22,dpi = 400)

#correlation of CD4 T cell and B cell 
immune<-lung_combine_after_batch[,which(lung_combine_after_batch$immune==TRUE)]
conv_infla_only<-as.data.frame(t(as.matrix(immune@assays[['conv']]@counts)))
cbind(conv_infla_only,immune@meta.data[,39:47])->conv_infla_only
summary(lm(`Cd4-tcell`~`cDC`,conv_infla_only[which(immune$orig.ident=='F_1'),]))
summary(lm(`Cd4-tcell`~`cDC`,conv_infla_only[which(immune$orig.ident=='F_2'),]))
summary(lm(`Cd4-tcell`~`cDC`,conv_infla_only[which(immune$orig.ident=='FS_1'),]))
summary(lm(`Cd4-tcell`~`cDC`,conv_infla_only[which(immune$orig.ident=='FS_2'),]))
set.seed(100)
summary(lm(`CD4_score1`~`B_score1`,conv_infla_only[which(immune$orig.ident=='F_1'),][sample(1:(dim(conv_infla_only[which(immune$orig.ident=='F_1'),])[1]),1000),]))
set.seed(100)
summary(lm(`CD4_score1`~`B_score1`,conv_infla_only[which(immune$orig.ident=='F_2'),][sample(1:(dim(conv_infla_only[which(immune$orig.ident=='F_2'),])[1]),1000),]))
set.seed(100)
summary(lm(`CD4_score1`~`B_score1`,conv_infla_only[which(immune$orig.ident=='FS_1'),][sample(1:(dim(conv_infla_only[which(immune$orig.ident=='FS_1'),])[1]),1000),]))
set.seed(100)
summary(lm(`CD4_score1`~`B_score1`,conv_infla_only[which(immune$orig.ident=='FS_2'),][sample(1:(dim(conv_infla_only[which(immune$orig.ident=='FS_2'),])[1]),1000),]))

cor(conv_infla_only[which(immune$orig.ident=='F_1'),c('cDC_marker_21','CD4_score1')],method = "kendall")
cor(conv_infla_only[which(immune$orig.ident=='F_2'),c('Cd4-tcell','B-cell')],method = "kendall")
cor(conv_infla_only[which(immune$orig.ident=='FS_1'),c('Cd4-tcell','B-cell')],method = "kendall")
cor(conv_infla_only[which(immune$orig.ident=='FS_2'),c('Cd4-tcell','B-cell')],method = "kendall")
conv_infla_only$`B-cell`<-0-log((conv_infla_only$`B-cell`+0.00001),base = 10)
conv_infla_only$`Cd4-tcell`<-0-log((conv_infla_only$`Cd4-tcell`+0.00001),base = 10)
conv_infla_only$`B-cell`<-sqrt(conv_infla_only$`B-cell`)
conv_infla_only$`Cd4-tcell`<-sqrt(conv_infla_only$`Cd4-tcell`)

conv_infla_only$`B-cell` <- asin(sqrt(conv_infla_only$`B-cell`))
conv_infla_only$`Cd4-tcell` <- asin(sqrt(conv_infla_only$`Cd4-tcell`))

conv_infla_only$`Cd4-tcell`<-1 / (conv_infla_only$`Cd4-tcell` + 1e-6)
conv_infla_only$`B-cell`<-1 / (conv_infla_only$`B-cell` + 1e-6)

conv_infla_only$CD4_bi<-'low_CD4'
conv_infla_only$CD4_bi[which(conv_infla_only$`Cd4-tcell`>0.02)]<-'high_CD4'
conv_infla_only$B_bi<-'low_B'
conv_infla_only$B_bi[which(conv_infla_only$`B-cell`>0.02)]<-'high_B'
chisq.test(table(conv_infla_only[which(immune$orig.ident=='FS_2'),"B_bi"],conv_infla_only[which(immune$orig.ident=='FS_2'),"CD4_bi"]))


conv_infla_only<-cbind(conv_infla_only,immune@meta.data[,39:42])

set.seed(100)
ggplot(conv_infla_only[which(immune$orig.ident=='FS_2'),], aes(x = `CD4_score1`, y = `macro_score1`)) +
  geom_point(size=0.5) +                             # Scatter plot (dot plot)
  geom_smooth(method = "lm", col = "red") + # Linear regression line
  theme_minimal() +
  theme(
    panel.grid = element_blank(),            # Remove grid lines
    panel.background = element_blank(),      # Remove background color
    axis.line = element_line(color = "black")# Add x and y axis lines
  )->p
cor(conv_infla_only[which(immune$orig.ident=='FS_2'),c('macro_score1','CD4_score1')])
summary(lm(macro_score1~CD4_score1,conv_infla_only[which(immune$orig.ident=='FS_2'),c('macro_score1','CD4_score1')]))$coefficients[2,4]
ggsave('~/ST_Superinfection/CD4_macro2_FS_2.png',p,width = 7,height = 7)
ggplot(conv[which(lung_combine_after_batch$orig.ident=='FS_1'),], aes(x = `CD4_score1`, y = `macro_score1`)) +
  geom_point(size=0.5) +                             # Scatter plot (dot plot)
  geom_smooth(method = "lm", col = "red") + # Linear regression line
  theme_minimal() +
  theme(
    panel.grid = element_blank(),            # Remove grid lines
    panel.background = element_blank(),      # Remove background color
    axis.line = element_line(color = "black")# Add x and y axis lines
  )->p
ggsave('~/ST_Superinfection/CD4_macro_FS_1.png',p,width = 7,height = 7)
cor(conv[which(lung_combine_after_batch$orig.ident=='FS_1'),c('macro_score1','CD4_score1')])



makeframe<-function(sample,tar1,tar2){
  conv_infla_only[which(immune$orig.ident==sample),]->s1
  set1 <- s1[which(s1[,tar1]<=0.01),tar2]  # 6 values
  set2 <- s1[which(s1[,tar1]>0.01),tar2]  # 8 values
  data <- data.frame(
    value = c(set1, set2),  # All values in one column
    group = c(rep(paste0('non-CD4-Tcell_',sample), length(set1)), rep(paste0('CD4-Tcell_',sample), length(set2)))  # Labels
  )
  return(data)
}
makeframe('S_1','Cd4-tcell','B-cell')->data1
makeframe('S_2','Cd4-tcell','B-cell')->data2
makeframe('F_1','Cd4-tcell','B-cell')->data3
makeframe('F_2','Cd4-tcell','B-cell')->data4
makeframe('FS_1','Cd4-tcell','B-cell')->data5
makeframe('FS_2','Cd4-tcell','B-cell')->data6
rbind(data3,data4)->data
data$group<-factor(data$group,levels = c("CD4-Tcell_F_1", "non-CD4-Tcell_F_1","CD4-Tcell_F_2", "non-CD4-Tcell_F_2"))
# Example data: two sets of values with different lengths
set1<-log((set1+1),2)
set1<-max(set1)-set1
set2<-max(set2)-set2
set2<-log((set2+1),2)
set2<-max(set2)-set2

# Create the box plot
ggplot(data[which(data$value>0.02),], aes(x = group, y = value)) +
  geom_violin(trim = FALSE, fill = "lightblue") +  # Create the violin plot
  geom_boxplot(width = 0.1, fill = "white") +  # Overlay boxplot for better visualization
  labs(title = "Distance from fibroblast to nearest macrophage", x = "Sample", y = "Distance to nearest macrophage") +
  theme_minimal() +
  ylim(c(0,0.35))+
  stat_compare_means(method = "wilcox.test",comparisons = list(c("CD4-Tcell_F_1", "non-CD4-Tcell_F_1"),c("CD4-Tcell_F_2", "non-CD4-Tcell_F_2")))+
  theme_minimal() +              # Use a minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(colour = "black"),
    text =element_text(size = 20), # Add axis lines in black
    axis.title.y = ( element_text(size=23,face = 'bold')),
    axis.title.x = ( element_text(size=23,face = 'bold'))
  )->p

# Display the plot
print(p)


#DC and T cell
cor(conv_infla_only[which(immune$orig.ident=='F_2'),c('cDC','Cd4-tcell','Cd8-tcell')])

FindMarkers(GSE202325_ref,group.by = 'celltype2',ident.1 = 'Cd4_tcell')->Cd4_marker
FindMarkers(GSE202325_ref,group.by = 'celltype2',ident.1 = 'B_cell')->B_marker

AddModuleScore(immune,list(c('Ccl5','Ms4a4b','Nkg7','Cd8a','Cd3d')),name = 'CD8_score')->immune
AddModuleScore(immune,features = list(c('Ccr7','Tcf7','Cd4','Il7r','Cd3d')),name = 'CD4_score')->immune
AddModuleScore(immune,features = list(c('Ccl5','Ms4a4b','Cxcr3','Ifng',"Tbx21","Gzmk")),name = 'Th1_score')->immune
AddModuleScore(immune,features = list(c('Ms4a1','Cd19','Ighm','Bank1','Pax5')),name = 'B_score')->immune
hist(immune$CD4_score1)


#find markers only in immune regions
DefaultAssay(immune)<-'reciever'
immune$macrophage<-immune[['conv']]$counts['Interstitial-macrophage',]
FindMarkers(immune,group.by = 'orig.ident',ident.1 = c('FS_1','FS_2'),ident.2 = c('F_1','F_2'),test.use = 'LR',latent.vars = c('immune_prop','macrophage'))->double_marker_reci_immune_viral
FindMarkers(immune,group.by = 'orig.ident',ident.1 = c('FS_1','FS_2'),ident.2 = c('S_1','S_2'),test.use = 'LR',latent.vars = c('immune_prop','macrophage'))->double_marker_reci_immune_bac
SpatialFeaturePlot(immune,'r.Ccl2.Ccr2',pt.size.factor = 2.1)->p
ggsave('p.png',p,width = 25,height = 15)
SpatialFeaturePlot(immune,'r.VCAM',pt.size.factor = 2.1)->p
ggsave('p2.png',p,width = 25,height = 15)

SpatialFeaturePlot(immune_exclude,'r.Tnf.Tnfrsf1a',pt.size.factor = 2.1)->p
ggsave('~/ST_Superinfection/r_Tnf_Tnfrsf1a.png',p,width = 25,height = 18,dpi = 400)
SpatialFeaturePlot(immune,'s.Tnf.Tnfrsf1a',pt.size.factor = 2.1)->p
ggsave('~/ST_Superinfection/s_Tnf_Tnfrsf1a.png',p,width = 25,height = 18,dpi = 400)

SpatialFeaturePlot(immune,'r.Cxcl5.Cxcr2',pt.size.factor = 1.8)->p
ggsave('~/ST_Superinfection/s_Cxcl1_Cxcr2.png',p,width = 25,height = 18,dpi = 400)
SpatialFeaturePlot(immune,'s.Tnf.Tnfrsf1a',pt.size.factor = 2.1)->p
ggsave('~/ST_Superinfection/s_Tnf_Tnfrsf1a.png',p,width = 25,height = 18,dpi = 400)

#exclude airway near
lung_combine_after_batch$immune_exclude_airway<-FALSE
lung_combine_after_batch$immune_exclude_airway[which(lung_combine_after_batch$immune==TRUE&lung_combine_after_batch$near_airway==FALSE)]<-TRUE


#cell colocalization
target1<-'CD4_score1'
target2<-'cDC_marker_11'

target1<-'CD8_score1'
target2<-'cDC_marker_21'

target1<-'CD4_score1'
target2<-'macro_score1'

target1<-'CD8_score1'
target2<-'macro_score1'

target1<-'Th1_score1'
target2<-'macro_score1'


target1<-'CD4_score1'
target2<-'B_score1'

mtx_CD4_B<-as.data.frame(matrix(nrow = 8,ncol = 3))
mtx_CD4_B$V2<-unique(immune$orig.ident)
mtx_CD4_B$V3<-substr(mtx_CD4_B$V2,1,2)
get_cor<-function(x){return(cor(conv_infla_only[which(immune$orig.ident==x),c(target1,target2)])[1,2])}
j<-1
for (x in unique(immune$orig.ident)) {
  get_cor(x)->mtx_CD4_B[j,1]
  j<-j+1
}

j<-1
for (x in unique(immune$orig.ident)) {
  0-log10(summary(lm(`CD4_score1`~`B_score1`,conv_infla_only[which(immune$orig.ident==x),]))$coefficients[2,4])->mtx_CD4_B$p[j]
  j<-j+1
}

mtx_CD4_B$V2<-factor(mtx_CD4_B$V2,levels = c('N_1','N_2','S_1','S_2','FS_2','FS_1','F_2','F_1'))
ggplot(mtx_CD4_B[5:8,], aes(x = V1, y = V2, color = p)) +
  geom_point(alpha = 0.7, size = 5) + # Adjust size and transparency as needed
  scale_color_gradient(low = "darkblue", high = "red") + # Adjust colors for gradient
  labs(x = "Value 1", y = "Value 2", color = "p") +
  geom_line(aes(group=V3))+
  theme_bw() +
  scale_x_continuous()+
  theme(
    legend.position = "right", # Position legend
    axis.text = element_text(size = 12), # Adjust axis text size
    axis.title = element_text(size = 14), # Adjust axis title size
    panel.grid.major.y = element_line(color = "grey"),
    legend.title = element_blank(),panel.grid=element_blank()
  )->p
ggsave('~/ST_Superinfection/cor_CD4_B.png',p,width = 7,height = 5.8,dpi = 400)


target1<-'CD4_score1'
target2<-'cDC_marker_11'

mtx_CD4_DC<-as.data.frame(matrix(nrow = 8,ncol = 3))
mtx_CD4_DC$V2<-unique(immune$orig.ident)
mtx_CD4_DC$V3<-substr(mtx_CD4_DC$V2,1,2)
get_cor<-function(x){return(cor(conv_infla_only[which(immune$orig.ident==x),c(target1,target2)])[1,2])}
j<-1
for (x in unique(immune$orig.ident)) {
  get_cor(x)->mtx_CD4_DC[j,1]
  j<-j+1
}

j<-1
for (x in unique(immune$orig.ident)) {
  0-log10(summary(lm(`CD4_score1`~`cDC_marker_11`,conv_infla_only[which(immune$orig.ident==x),]))$coefficients[2,4])->mtx_CD4_DC$p[j]
  j<-j+1
}

mtx_CD4_DC$V2<-factor(mtx_CD4_DC$V2,levels = c('N_1','N_2','S_1','S_2','FS_2','FS_1','F_2','F_1'))
ggplot(mtx_CD4_DC[5:8,], aes(x = V1, y = V2, color = p)) +
  geom_point(alpha = 0.7, size = 5) + # Adjust size and transparency as needed
  scale_color_gradient(low = "darkblue", high = "red") + # Adjust colors for gradient
  labs(x = "Value 1", y = "Value 2", color = "p") +
  geom_line(aes(group=V3))+
  theme_bw() +
  scale_x_continuous()+
  theme(
    legend.position = "right", # Position legend
    axis.text = element_text(size = 12), # Adjust axis text size
    axis.title = element_text(size = 14), # Adjust axis title size
    panel.grid.major.y = element_line(color = "grey"),
    legend.title = element_blank(),panel.grid=element_blank()
  )->p
ggsave('~/ST_Superinfection/cor_CD4_DC.png',p,width = 7,height = 5.8,dpi = 400)





target1<-'CD4_score1'
target2<-'neu_score1'

mtx_CD8_macro<-as.data.frame(matrix(nrow = 8,ncol = 3))
mtx_CD8_macro$V2<-unique(immune$orig.ident)
mtx_CD8_macro$V3<-substr(mtx_CD8_macro$V2,1,2)
get_cor<-function(x){return(cor(conv_infla_only[which(immune$orig.ident==x),c(target1,target2)])[1,2])}
j<-1
for (x in unique(immune$orig.ident)) {
  get_cor(x)->mtx_CD8_macro[j,1]
  j<-j+1
}

j<-1
for (x in unique(immune$orig.ident)) {
  0-log10(summary(lm(`CD8_score1`~`macro_score1`,conv_infla_only[which(immune$orig.ident==x),]))$coefficients[2,4])->mtx_CD8_macro$p[j]
  j<-j+1
}

mtx_CD8_macro$V2<-factor(mtx_CD8_macro$V2,levels = c('N_1','N_2','S_1','S_2','FS_2','FS_1','F_2','F_1'))
ggplot(mtx_CD8_macro[5:8,], aes(x = V1, y = V2, color = p)) +
  geom_point(alpha = 0.7, size = 3) + # Adjust size and transparency as needed
  scale_color_gradient(low = "darkblue", high = "red") + # Adjust colors for gradient
  labs(x = "Value 1", y = "Value 2", color = "p") +
  geom_line(aes(group=V3))+
  theme_bw() +
  scale_x_continuous()+
  theme(
    legend.position = "right", # Position legend
    axis.text = element_text(size = 12), # Adjust axis text size
    axis.title = element_text(size = 14), # Adjust axis title size
    panel.grid.major.y = element_line(color = "grey"),
    legend.title = element_blank(),panel.grid=element_blank()
  )->p
ggsave('~/ST_Superinfection/cor_CD4_DC.png',p,width = 7,height = 5.8,dpi = 400)

#correlation plot
# Define breaks to use the same color scale across all heatmaps
cor(conv_infla_only[which(immune$orig.ident=='N_1'),c('Alveolar-macrophage','Cd4-tcell','Cd8-tcell','cDC','B-cell','IgA-plasma-cell','Interstitial-macrophage','Mast-cell','Monocyte','Neutrophil','NK','Plasma-cell')])->mtx_N_1
cor(conv_infla_only[which(immune$orig.ident=='N_2'),c('Alveolar-macrophage','Cd4-tcell','Cd8-tcell','cDC','B-cell','IgA-plasma-cell','Interstitial-macrophage','Mast-cell','Monocyte','Neutrophil','NK','Plasma-cell')])->mtx_N_2
cor(conv_infla_only[which(immune$orig.ident=='S_1'),c('Alveolar-macrophage','Cd4-tcell','Cd8-tcell','cDC','B-cell','IgA-plasma-cell','Interstitial-macrophage','Mast-cell','Monocyte','Neutrophil','NK','Plasma-cell')])->mtx_S_1
cor(conv_infla_only[which(immune$orig.ident=='S_2'),c('Alveolar-macrophage','Cd4-tcell','Cd8-tcell','cDC','B-cell','IgA-plasma-cell','Interstitial-macrophage','Mast-cell','Monocyte','Neutrophil','NK','Plasma-cell')])->mtx_S_2
cor(conv_infla_only[which(immune$orig.ident=='F_1'),c('Alveolar-macrophage','Cd4-tcell','Cd8-tcell','cDC','B-cell','IgA-plasma-cell','Interstitial-macrophage','Mast-cell','Monocyte','Neutrophil','NK','Plasma-cell')])->mtx_F_1
cor(conv_infla_only[which(immune$orig.ident=='F_2'),c('Alveolar-macrophage','Cd4-tcell','Cd8-tcell','cDC','B-cell','IgA-plasma-cell','Interstitial-macrophage','Mast-cell','Monocyte','Neutrophil','NK','Plasma-cell')])->mtx_F_2
cor(conv_infla_only[which(immune$orig.ident=='FS_1'),c('Alveolar-macrophage','Cd4-tcell','Cd8-tcell','cDC','B-cell','IgA-plasma-cell','Interstitial-macrophage','Mast-cell','Monocyte','Neutrophil','NK','Plasma-cell')])->mtx_FS_1
cor(conv_infla_only[which(immune$orig.ident=='FS_2'),c('Alveolar-macrophage','Cd4-tcell','Cd8-tcell','cDC','B-cell','IgA-plasma-cell','Interstitial-macrophage','Mast-cell','Monocyte','Neutrophil','NK','Plasma-cell')])->mtx_FS_2

diag(mtx_N_1)<-0
mtx_N_1[which(abs(mtx_N_1)<0.3)]<-0
chordDiagram(mtx_N_1, transparency = 0.5)
diag(mtx_N_2)<-0
mtx_N_2[which(abs(mtx_N_2)<0.3)]<-0
chordDiagram(mtx_N_2, transparency = 0.5)
chordDiagram(mtx_S_1, transparency = 0.5)
chordDiagram(mtx_S_2, transparency = 0.5)
chordDiagram(mtx_F_1, transparency = 0.5)
chordDiagram(mtx_F_2, transparency = 0.5)
chordDiagram(mtx_FS_1, transparency = 0.5)
chordDiagram(mtx_FS_2, transparency = 0.5)

global_min<-min(c(mtx_N_1,mtx_N_2,mtx_S_1,mtx_S_2,mtx_F_1,mtx_F_2,mtx_FS_1,mtx_FS_2))
breaks <- seq(global_min, 1, length.out = 101)
colorRampPalette(c("navy", "white", "firebrick3"))(100)->color_palette
# Plot each heatmap with the same color scale
pheatmap::pheatmap(mtx_F_1, color = color_palette, breaks = breaks)->p
p$tree_col$order->order
pheatmap::pheatmap(mtx_F_1[order,order], color = color_palette, breaks = breaks,cluster_rows = FALSE,cluster_cols = FALSE,fontsize_row = 20,fontsize_col = 20)->p_F_1
pheatmap::pheatmap(mtx_F_2[order,order], color = color_palette, breaks = breaks,cluster_rows = FALSE,cluster_cols = FALSE,fontsize_row = 20,fontsize_col = 20)->p_F_2
pheatmap::pheatmap(mtx_N_1[order,order], color = color_palette, breaks = breaks,cluster_rows = FALSE,cluster_cols = FALSE,fontsize_row = 20,fontsize_col = 20)->p_N_1
pheatmap::pheatmap(mtx_N_2[order,order], color = color_palette, breaks = breaks,cluster_rows = FALSE,cluster_cols = FALSE,fontsize_row = 20,fontsize_col = 20)->p_N_2
pheatmap::pheatmap(mtx_S_1[order,order], color = color_palette, breaks = breaks,cluster_rows = FALSE,cluster_cols = FALSE,fontsize_row = 20,fontsize_col = 20)->p_S_1
pheatmap::pheatmap(mtx_S_2[order,order], color = color_palette, breaks = breaks,cluster_rows = FALSE,cluster_cols = FALSE,fontsize_row = 20,fontsize_col = 20)->p_S_2
pheatmap::pheatmap(mtx_FS_1[order,order], color = color_palette, breaks = breaks,cluster_rows = FALSE,cluster_cols = FALSE,fontsize_row = 20,fontsize_col = 20)->p_FS_1
pheatmap::pheatmap(mtx_FS_2[order,order], color = color_palette, breaks = breaks,cluster_rows = FALSE,cluster_cols = FALSE,fontsize_row = 20,fontsize_col = 20)->p_FS_2
library(ggplot2)
ggsave('~/ST_Superinfection/heatmap/p_F_1.png',p_F_1,width = 10,height = 10,dpi = 300)
ggsave('~/ST_Superinfection/heatmap/p_F_2.png',p_F_2,width = 10,height = 10,dpi = 300)
ggsave('~/ST_Superinfection/heatmap/p_S_1.png',p_S_1,width = 10,height = 10,dpi = 300)
ggsave('~/ST_Superinfection/heatmap/p_S_2.png',p_S_2,width = 10,height = 10,dpi = 300)
ggsave('~/ST_Superinfection/heatmap/p_N_1.png',p_N_1,width = 10,height = 10,dpi = 300)
ggsave('~/ST_Superinfection/heatmap/p_N_2.png',p_N_2,width = 10,height = 10,dpi = 300)
ggsave('~/ST_Superinfection/heatmap/p_FS_1.png',p_FS_1,width = 10,height = 10,dpi = 300)
ggsave('~/ST_Superinfection/heatmap/p_FS_2.png',p_FS_2,width = 10,height = 10,dpi = 300)




pheatmap(mat2, color = color_palette, breaks = breaks)
pheatmap(mat3, color = color_palette, breaks = breaks)

cor(conv_infla_only[which(immune$orig.ident=='F_2'),c('Alveolar-macrophage','Cd4-tcell','Cd8-tcell','cDC','B-cell','IgA-plasma-cell','Interstitial-macrophage','Mast-cell','Monocyte','Neutrophil','NK','Plasma-cell')])->mtx
pheatmap::pheatmap(mtx)

breaks <- seq(global_min, global_max, length.out = 101)



