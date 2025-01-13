# Functions to calculate distance to airway
#calculation of distance was based on the union annotation (deconvolution + manual histology) of airway regions

get_distance_airway<-function(object,image,celltype){
  object@images[[image]]@coordinates->mtx
  mtx[which(object$airway==TRUE),]->mtx_air
  k<-1
  out<-as.data.frame(matrix(nrow = length(which(object$first_type==celltype)),ncol = dim(mtx_air)[1]))
  for (i in which(object$first_type==celltype)) {
    as.data.frame(slice(mtx,rep(i,dim(mtx_air)[1]))[,c(4,5)]-mtx_air[,c(4,5)])->m
    apply(m,1,FUN = function(x){sqrt(sum(x^2))})->out[k,]
    k<-k+1
  }
  colnames(out)<-rownames(mtx_air)
  rownames(out)<-colnames(object)[which(object$first_type==celltype)]
  return(out)
}

get_distance_airway_loop<-function(sample,image){
  distance_list<-list()
  k<-1
  for (i in unique(lung_combine_after_batch$first_type)) {
    get_distance_airway(lung_combine_after_batch[,which(lung_combine_after_batch$orig.ident==sample)],image,i)->distance_list[[k]]
    k<-k+1
  }
  names(distance_list)<-unique(lung_combine_after_batch$first_type)
  distance_list<-do.call(rbind,distance_list)
  distance_list<-as.data.frame(distance_list)
  distance_list$cell_type<-sub("[.].*", "", rownames(distance_list))
  return(distance_list)
}

#calculate distance
library(tidyverse)
get_distance_airway_loop('N_1','slice1')->distance_list_N_1
get_distance_airway_loop('N_2','slice1.2')->distance_list_N_2
get_distance_airway_loop('S_1','slice1.3')->distance_list_S_1
get_distance_airway_loop('S_2','slice1.4')->distance_list_S_2
get_distance_airway_loop('F_1','slice1.5')->distance_list_F_1
get_distance_airway_loop('F_2','slice1.6')->distance_list_F_2
get_distance_airway_loop('FS_1','slice1.7')->distance_list_FS_1
get_distance_airway_loop('FS_2','slice1.8')->distance_list_FS_2

#minimum distance to airway regions
min2<-function(object){
  return(as.array(apply(object[,-dim(object)[2]], 1, FUN = min)))
}
min_distance<-c(min2(distance_list_N_1),min2(distance_list_N_2),min2(distance_list_S_1),min2(distance_list_S_2),min2(distance_list_F_1),min2(distance_list_F_2),min2(distance_list_FS_1),min2(distance_list_FS_2))

combine_distance<-data.frame(barcode=c(rownames(distance_list_N_1),rownames(distance_list_N_2),rownames(distance_list_S_1),rownames(distance_list_S_2),rownames(distance_list_F_1),rownames(distance_list_F_2),rownames(distance_list_FS_1),rownames(distance_list_FS_2)),distance_min=min_distance)
rownames(combine_distance)<-sub('.*[.]','',combine_distance$barcode)
lung_combine_after_batch$distance_min<-combine_distance$distance_min[match(colnames(lung_combine_after_batch),rownames(combine_distance))]





