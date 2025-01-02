AddModuleScore(lung_combine_after_batch,list(c('Ccl5','Ms4a4b','Nkg7','Cd8a','Cd3d')),name = 'CD8_score')->lung_combine_after_batch
AddModuleScore(lung_combine_after_batch,features = list(c('Ccr7','Tcf7','Cd4','Il7r','Cd3d')),name = 'CD4_score')->lung_combine_after_batch
AddModuleScore(lung_combine_after_batch,features = list(c('Ms4a1','Cd19','Ighm','Bank1','Pax5')),name = 'B_score')->lung_combine_after_batch
AddModuleScore(lung_combine_after_batch,features = list(c('Itgae','Xcr1','Clec9a','Itgax','Cd74','H2-Eb1','H2-Ab1','H2-Aa')),name = 'cDC_marker_1')->lung_combine_after_batch
AddModuleScore(lung_combine_after_batch,features = list(rownames(DC_marker)[1:6]),name = 'cDC_marker_2')->lung_combine_after_batch
AddModuleScore(lung_combine_after_batch,features = list(c('Csf3r','Ly6g','S100a9','S100a8','Retnlg')),name = 'neu_score')->lung_combine_after_batch

lung_combine_after_batch$Th17_score1<-colMeans(as.data.frame(lung_combine_after_batch@assays$Spatial@data[c("Rorc","Il17a",'Ccr6','Il22','Il17f'),]))
lung_combine_after_batch$Th1_score1<-colMeans(as.data.frame(lung_combine_after_batch@assays$Spatial@data[c('Ccl5','Ms4a4b','Cxcr3','Ifng',"Tbx21","Gzmk"),]))