#annotation for airway
lung_combine_after_batch$airway<-FALSE
lung_combine_after_batch$airway[which(lung_combine_after_batch$airway_manual==TRUE|lung_combine_after_batch$first_type%in%c("Ciliated_cell","Club_cell","Goblst_cell"))]<-TRUE
table(lung_combine_after_batch$airway)

#define inflammatory regions
conv<-as.data.frame(lung_combine_after_batch[['conv']]$counts)
conv_infla<-conv[c('Alveolar-macrophage','Cd4-tcell','Cd8-tcell','cDC','B-cell','IgA-plasma-cell','Interstitial-macrophage','Macrophage','Mast-cell','mitotic-tcell','Monocyte','Neutrophil','NK','Plasma-cell'),]
lung_combine_after_batch$immune_prop<-colSums(conv_infla)[match(colnames(lung_combine_after_batch),colnames(conv_infla))]
lung_combine_after_batch$immune<-FALSE
lung_combine_after_batch$immune[which(lung_combine_after_batch$immune_prop>0.35)]<-TRUE
lung_combine_after_batch$ALL<-TRUE
SpatialDimPlot(lung_combine_after_batch,group.by = 'immune')->p
ggsave('p.png',p,width = 35,height = 22,dpi = 400)


#define follicle-like immune hubs
lung_combine_after_batch$immune_hub<-FALSE
hist(colSums(lung_combine_after_batch[['conv']]$counts[c("Cd4-tcell","Cd8-tcell","B-cell","cDC"),]))
hist(lung_combine_after_batch[['Spatial']]$data['Cxcl13',])
lung_combine_after_batch$immune_hub[which(lung_combine_after_batch[['Spatial']]$data['Cxcl13',]>0.6&colSums(as.matrix(lung_combine_after_batch[['conv']]$counts[c("Cd4-tcell","Cd8-tcell","B-cell","cDC"),]))>0.1)]<-TRUE
