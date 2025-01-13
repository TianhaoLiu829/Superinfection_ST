library(Seurat)
library(ggplot2)
#Figure 6A
SpatialFeaturePlot(lung_combine_after_batch,'Cxcl13',images = c('slice1.5','slice1.6','slice1.7','slice1.8'))
#Figure 6B
#see the annotation for immune hubs in define_regions.R
SpatialDimPlot(lung_combine_after_batch,'immune_hub',images = c('slice1.5','slice1.6','slice1.7','slice1.8'))

#Figure 6D
# Sample data
as.data.frame(lung_combine_after_batch_archieve@assays$reciever@counts)->df4
df4<-as.data.frame(t(df4))
df3 <- cbind(df4$r.Cxcl13.Cxcr5,lung_combine_after_batch$orig.ident)
colnames(df3)<-c('value','label')
as.data.frame(df3)->df3
df3$value<-as.numeric(df3$value)
# Calculate count of non-zero values and log mean for each label
summary_df <- df3 %>%
  filter(value != 0) %>%
  group_by(label) %>%
  summarize(
    count_non_zero = n(),
    log_mean_value = mean(value)
  )
summary_df$count_non_zero<-as.numeric(summary_df$count_non_zero/table(df3$label)[match(summary_df$label,names(table(df3$label)))])
summary_df$label<-factor(summary_df$label,levels = c('N_1','N_2','S_1','S_2','F_1','F_2','FS_1','FS_2')[8:1])
# Plot
ggplot(summary_df, aes(x = 1, y = label, size = count_non_zero, color = log_mean_value)) +
  geom_point() +
  scale_x_continuous(breaks = NULL)+
  scale_size_continuous(range = c(3, 15)) +  # Adjust size range as needed
  scale_color_viridis_c() +                  # Color scale for log mean
  labs(x = "Sample Label", y = "", size = "Non-Zero Count", color = "Log Mean Value") +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.grid=element_blank(),
    axis.text.x = element_blank(),  # Adjust x-axis text size
    axis.text.y = element_text(size = 14),   # Adjust y-axis text size
    axis.title.x = element_text(size = 25),  # Adjust x-axis text size
    axis.title.y = element_text(size = 25)   # Adjust y-axis text size
  )+
  theme(     # Remove y-axis ticks
    axis.title.y = element_blank(),
    legend.position = "right") +    # Remove y-axis title
  ggtitle("Dot Plot of Non-Zero Counts and Log Mean Value by Label")->p
p


#Figure 6E
ibalt<-lung_combine_after_batch[,which(lung_combine_after_batch$ibalt)]
as.data.frame(ibalt@assays$sender@counts+ibalt@assays$reciever@counts)->ibalt_all_sender
ibalt_all_sender[,which(ibalt$orig.ident%in%c('F_1','F_2'))]->ibalt_F_sender
ibalt_all_sender[,which(ibalt$orig.ident%in%c('FS_1','FS_2'))]->ibalt_FS_sender
as.data.frame(matrix(ncol = 6,nrow = dim(ibalt_F)[1]))->ibalt_marker_sender
for (i in 1:dim(ibalt_F_sender)[1]) {
  wilcox.test(as.matrix(ibalt_FS_sender[i,]),as.matrix(ibalt_F_sender[i,]))$p.value->ibalt_marker_sender[i,1]
  log(mean(as.matrix(ibalt_FS_sender[i,]))/mean(as.matrix(ibalt_F_sender[i,])),2)->ibalt_marker_sender[i,2]
  mean(as.matrix(ibalt_all_sender[i,which(ibalt$orig.ident=='F_1')]))->ibalt_marker_sender[i,3]
  mean(as.matrix(ibalt_all_sender[i,which(ibalt$orig.ident=='F_2')]))->ibalt_marker_sender[i,4]
  mean(as.matrix(ibalt_all_sender[i,which(ibalt$orig.ident=='FS_1')]))->ibalt_marker_sender[i,5]
  mean(as.matrix(ibalt_all_sender[i,which(ibalt$orig.ident=='FS_2')]))->ibalt_marker_sender[i,6]
}
rownames(ibalt_marker_sender)<-rownames(ibalt_all_sender)
colnames(ibalt_marker_sender)<-c('pvalue','log2FC','F_1','F_2','FS_1','FS_2')
ibalt_marker_sender$padj<-p.adjust(ibalt_marker_sender$pvalue,method = 'fdr')
ibalt_marker_sender$mean_all_sender<-rowMeans(ibalt_all_sender)
as.data.frame(t(ibalt_all_sender))->ibalt_all_sender2 
agg<-function(x){
  aggregate(x,by=list(ibalt$orig.ident),FUN=function(y){return(sum(y>0.1)/length(y))})->m
  m<-as.data.frame(m)
  return(as.data.frame(t(m))[2,1:4])
}
do.call(rbind,apply(ibalt_all_sender2, 2, agg))->ibalt_all_sender2
cbind(ibalt_marker_sender,ibalt_all_sender2)->ibalt_marker_sender
colnames(ibalt_marker_sender)[9:12]<-paste0(colnames(ibalt_marker_sender)[3:6],'_count')
ibalt_marker_sender<-ibalt_marker_sender[order(ibalt_marker_sender$pvalue),]
ibalt_marker_sender_ll<-ibalt_marker_sender[grepl("CXCL|CCL|Cxcl|Ccl",rownames(ibalt_marker_sender)),]
ibalt_marker_sender_ll<-ibalt_marker_sender_ll[order(ibalt_marker_sender_ll$pvalue),]
ibalt_marker_sender_ll$Molecule<-rownames(ibalt_marker_sender_ll)
library(tidyverse)
library(reshape2)
mean_data <- melt(ibalt_marker_sender_ll, id.vars = "Molecule", measure.vars = c("F_1", "F_2", "FS_1", "FS_2"),
                  variable.name = "Group", value.name = "Mean")
prop_data <- melt(ibalt_marker_sender_ll, id.vars = "Molecule", measure.vars = c("F_1_count", "F_2_count", "FS_1_count", "FS_2_count"),
                  variable.name = "Group", value.name = "Proportion")

# Step 2: Combine mean and proportion data
mean_data$Proportion <- prop_data$Proportion
mean_data$Group <- gsub("Mean_", "", mean_data$Group)  # Clean group names
mean_data$Proportion<-as.numeric(mean_data$Proportion)


mean_data_plot<-mean_data[which(mean_data$Molecule%in%c('s.Cxcl13.Cxcr5',rownames(ibalt_marker_sender_ll)[which(ibalt_marker_sender_ll$log2FC<0)][c(3:11)])),]

scaled_data <- mean_data_plot %>%
  group_by(Molecule) %>%                          # Group by the label
  mutate(scaled_value = scale(Mean)) %>%      # Standardize values
  ungroup()                                    # Remove grouping
colnames(scaled_data)

scaled_data$Group<-factor(scaled_data$Group,levels = c('FS_2','FS_1','F_2','F_1'))
scaled_data$Molecule<-factor(scaled_data$Molecule,levels = c('s.Cxcl13.Cxcr5',rownames(ibalt_marker_sender_ll)[c(3:11)]))
ggplot(scaled_data, aes(x = Molecule, y = Group, size = Proportion, color =scaled_value)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(3, 12)) +  # Bubble size range
  scale_color_viridis_c() +  # Color gradient for mean values
  theme_bw() +
  theme( 
    legend.position="none",
    panel.grid=element_blank(),
    axis.text.x = element_text(size = 14),  # Adjust x-axis text size
    axis.text.y = element_text(size = 14),   # Adjust y-axis text size
    axis.title.x = element_text(size = 25),  # Adjust x-axis text size
    axis.title.y = element_text(size = 25)   # Adjust y-axis text size
  )+
  theme(     # Remove y-axis ticks
        axis.title.y = element_blank(),
        legend.position = "right") +    # Remove y-axis title
  ggtitle("Dot Plot of Non-Zero Counts and Log Mean Value by Label")->p
p

