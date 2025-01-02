#proportion plot
df<-as.data.frame(lung_combine_after_batch@meta.data[,c(1,20)])
library(dplyr)

# Calculate proportions
df_proportion <- df %>%
  group_by(orig.ident, seurat_clusters) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))

library(ggplot2)
df_proportion$orig.ident<-factor(df_proportion$orig.ident,levels = c('N_1','N_2','S_1','S_2','F_1','F_2','FS_1','FS_2'))
ggplot(df_proportion, aes(x = orig.ident, y = proportion, fill = seurat_clusters2)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(
    x = "Label1",
    y = "Proportion",
    fill = "Label2",
    title = "Proportion of Label2 within each Label1"
  ) +
  theme_bw() +
  theme(legend.title = element_blank(),panel.grid=element_blank())



as.data.frame(lung_combine_after_batch@assays$conv@counts)->s
sweep(s,2, colSums(s),'/')->s
library(tidyr)
gather(s)->s1
s1$celltype<-rep(rownames(s),dim(s)[2])
s1$domain<-lung_combine_after_batch$seurat_clusters[match(s1$key,colnames(lung_combine_after_batch))]
aggregate(s1$value,by=list(s1$domain,s1$celltype),FUN=mean)->k
aggregate(k$x,by=list(k$Group.1),FUN=sum)
colnames(k)<-c('spatial_domain','cell_type','proportion')
library(plyr)
ggplot(data=k,mapping = aes(x=spatial_domain ,y=proportion,fill=cell_type))+geom_bar(stat = 'identity',width = 0.6)+theme_classic()+theme(text = element_text(size = 20),axis.title.y = element_text(size = 20,face = 'bold') ,axis.text.x = element_text(angle = 50,hjust = 1,size = 15,face='bold'))->p
ggsave('~/ST_Superinfection/prop.png',p,width = 20,height = 7,dpi = 300)

a<-array()
b<-array()
c<-array()
for (i in 0:11) {
  k$cell_type[which(k$spatial_domain==i)][order(0-k$proportion[which(k$spatial_domain==i)])[1]]->a[i+1]
}

for (i in 0:11) {
  k$cell_type[which(k$spatial_domain==i)][order(0-k$proportion[which(k$spatial_domain==i)])[2]]->b[i+1]
}
paste0(a,'_',b)

for (i in 1:length(unique(k$cell_type))) {
  k$spatial_domain[which(k$cell_type==unique(k$cell_type)[i])][order(0-k$proportion[which(k$cell_type==unique(k$cell_type)[i])])[2]]->c[i]
}
lapply(1:29,FUN = function(x){return(unique(k$cell_type)[which(c==x)])})->m

from <- c(0,2,1,3,11,8,10,5,9,4,6,7)
to <- c(1:12)

# Mutate labels
mapvalues(lung_combine_after_batch$seurat_clusters, from = from, to = to)->lung_combine_after_batch$seurat_clusters2


#deconvolution_cluster_plot


