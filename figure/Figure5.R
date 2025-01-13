#
immune<-lung_combine_after_batch[,which(lung_combine_after_batch$immune==TRUE)]
conv_infla_only<-as.data.frame(t(as.matrix(immune@assays[['conv']]@counts)))
AddModuleScore(lung_combine_after_batch,features = list(c('Ccr7','Tcf7','Cd4','Il7r','Cd3d')),name = 'CD4_score')->lung_combine_after_batch
AddModuleScore(lung_combine_after_batch,features = list(c('Ms4a1','Cd19','Ighm','Bank1','Pax5')),name = 'B_score')->lung_combine_after_batch
AddModuleScore(lung_combine_after_batch,features = list(c('Itgae','Xcr1','Clec9a','Itgax','Cd74','H2-Eb1','H2-Ab1','H2-Aa')),name = 'cDC_marker')->lung_combine_after_batch
#Figure 5A
plot_list_CD4_B<-list()
for (i in 1:4) {
  ggplot(conv_infla_only[which(immune$orig.ident==c("F_1","F_2","FS_1","FS_2")[i]),], aes(x = CD4_score1, y = B_score1)) +
    geom_point(size=0.5) +                             # Scatter plot (dot plot)
    geom_smooth(method = "lm", col = "red") + # Linear regression line
    theme_minimal() +
    theme(
      panel.grid = element_blank(),            # Remove grid lines
      panel.background = element_blank(),      # Remove background color
      axis.line = element_line(color = "black")# Add x and y axis lines
    )->plot_list_CD4_B[[i]]
}
library(patchwork)
plot_list_CD4_B[[1]] + plot_list_CD4_B[[2]] + plot_list_CD4_B[[3]] + plot_list_CD4_B[[4]]
get_cor<-function(x){return(cor(conv_infla_only[which(immune$orig.ident==x),c("CD4_score1","B_score1")])[1,2])}
mtx_CD4_B<-data.frame(correlation=NA,sample=unique(immune$orig.ident),condition=substr(unique(immune$orig.ident),1,2))
j<-1
for (x in unique(immune$orig.ident)) {
  get_cor(x)->mtx_CD4_B[j,1]
  j<-j+1
}
mtx_CD4_B

#Figure 5B
plot_list_CD4_cDC<-list()
for (i in 1:4) {
  ggplot(conv_infla_only[which(immune$orig.ident==c("F_1","F_2","FS_1","FS_2")[i]),], aes(x = CD4_score1, y = cDC_marker1)) +
    geom_point(size=0.5) +                             # Scatter plot (dot plot)
    geom_smooth(method = "lm", col = "red") + # Linear regression line
    theme_minimal() +
    theme(
      panel.grid = element_blank(),            # Remove grid lines
      panel.background = element_blank(),      # Remove background color
      axis.line = element_line(color = "black")# Add x and y axis lines
    )->plot_list_CD4_cDC[[i]]
}
library(patchwork)
plot_list_CD4_cDC[[1]] + plot_list_CD4_cDC[[2]] + plot_list_CD4_cDC[[3]] + plot_list_CD4_cDC[[4]]
get_cor<-function(x){return(cor(conv_infla_only[which(immune$orig.ident==x),c("CD4_score1","cDC_marker1")])[1,2])}
mtx_CD4_DC<-data.frame(correlation=NA,sample=unique(immune$orig.ident),condition=substr(unique(immune$orig.ident),1,2))
j<-1
for (x in unique(immune$orig.ident)) {
  get_cor(x)->mtx_CD4_DC[j,1]
  j<-j+1
}
mtx_CD4_DC

#Figure 5C
SpatialFeaturePlot(lung_combine_after_batch,'B-cell',images = c('slice1.5','slice1.6','slice1.7','slice1.8'))
SpatialFeaturePlot(lung_combine_after_batch,'Cd4-tcell',images = c('slice1.5','slice1.6','slice1.7','slice1.8'))

#Figure 5D
as.data.frame(lung_combine_after_batch@assays$reciever@counts+lung_combine_after_batch@assays$sender@counts)->df4
df4<-as.data.frame(t(df4))
df3 <- cbind(df4[,c("r.ICOS","r.CD40")],lung_combine_after_batch$orig.ident)
colnames(df3)<-c('r.ICOS','r.CD40','label')
as.data.frame(df3)->df3
df3<-df3[which(lung_combine_after_batch$immune==TRUE),]
# Summarize the data: calculate proportion of non-zero values and mean for each label
summary_df <- df3 %>%
  group_by(label) %>%
  summarise(
    Prop_NonZero_r.ICOS = mean(r.ICOS != 0),
    Prop_NonZero_r.CD40 = mean(r.CD40 != 0),
    Mean_r.ICOS = mean(r.ICOS, na.rm = TRUE),
    Mean_r.CD40 = mean(r.CD40, na.rm = TRUE)
  )

# View the summarized data
print(summary_df)
# Reshape data to long format for ggplot
long_df <- summary_df %>%
  pivot_longer(cols = starts_with("Prop_NonZero"), 
               names_to = c("Set", "Metric"), 
               names_pattern = "(.*)_(.*)",
               values_to = "Proportion")
long_df<-long_df[which(long_df$label%in%c('F_1','F_2','FS_1','FS_2')),]
# Create a bubble plot
ggplot(long_df[which(long_df$Metric=='r.ICOS'),], aes(y = Metric, x = label, size = Proportion, color = Mean_r.ICOS)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(3, 15)) +  # Adjust the size range of bubbles
  scale_color_viridis_c()+  # Adjust color scale based on mean
  theme_bw() +
  theme( 
    legend.position="none",
    panel.grid=element_blank(),
    axis.text.x = element_text(size = 14),  # Adjust x-axis text size
    axis.text.y = element_text(size = 14),   # Adjust y-axis text size
    axis.title.x = element_text(size = 25),  # Adjust x-axis text size
    axis.title.y = element_text(size = 25)   # Adjust y-axis text size
  )+
  theme(    # Remove y-axis ticks
    axis.title.y = element_blank(),
    legend.position = "right") +    # Remove y-axis title
  ggtitle("Dot Plot of Non-Zero Counts and Log Mean Value by Label")


#Figure 5E
lung_combine_after_batch$Th17_score1<-colMeans (as.data.frame(lung_combine_after_batch@assays$Spatial@data[c("Rorc","Il17a",'Ccr6','Il22','Il17f'),]))
lung_combine_after_batch$Th1_score1<-colMeans(as.data.frame(lung_combine_after_batch@assays$Spatial@data[c('Ccl5','Ms4a4b','Cxcr3','Ifng',"Tbx21","Gzmk"),]))
# Sample data
as.data.frame(lung_combine_after_batch$Th1_score1)->df4
df3 <- cbind(df4[,1],lung_combine_after_batch$orig.ident)
colnames(df3)<-c('value','label')
as.data.frame(df3)->df3
df3$value<-as.numeric(df3$value)
df3<-df3[which(lung_combine_after_batch$near_airway ==TRUE),]
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
  theme(axis.text.y = element_blank(), 
        panel.grid=element_blank(),
        axis.ticks.y = element_blank(),      # Remove y-axis ticks
        axis.title.y = element_blank(),
        legend.position = "right") +    # Remove y-axis title
  ggtitle("Dot Plot of Non-Zero Counts and Log Mean Value by Label")






