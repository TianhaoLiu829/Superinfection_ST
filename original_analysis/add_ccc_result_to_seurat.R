#read sending signal csv file
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/send_N_1.csv')->send_N_1
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/send_N_2.csv')->send_N_2
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/send_S_1.csv')->send_S_1
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/send_S_2.csv')->send_S_2
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/send_F_1.csv')->send_F_1
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/send_F_2.csv')->send_F_2
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/send_FS_1.csv')->send_FS_1
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/send_FS_2.csv')->send_FS_2
send_N_1$X<-paste0('N_1_',send_N_1$X)
send_N_2$X<-paste0('N_2_',send_N_2$X)
send_S_1$X<-paste0('S_1_',send_S_1$X)
send_S_2$X<-paste0('S_2_',send_S_2$X)
send_F_1$X<-paste0('F_1_',send_F_1$X)
send_F_2$X<-paste0('F_2_',send_F_2$X)
send_FS_1$X<-paste0('FS_1_',send_FS_1$X)
send_FS_2$X<-paste0('FS_2_',send_FS_2$X)
rbind(send_N_1,send_N_2,send_S_1,send_S_2,send_F_1,send_F_2,send_FS_1,send_FS_2)->send

#read receiving signal csv file
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/recie4_N_1.csv')->recie_N_1
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/recie4_N_2.csv')->recie_N_2
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/recie4_S_1.csv')->recie_S_1
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/recie4_S_2.csv')->recie_S_2
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/recie4_F_1.csv')->recie_F_1
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/recie4_F_2.csv')->recie_F_2
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/recie4_FS_1.csv')->recie_FS_1
read.csv('/ix1/wchen/liutianhao/result/lung_ST/result/recie4_FS_2.csv')->recie_FS_2
recie_N_1$X<-paste0('N_1_',recie_N_1$X)
recie_N_2$X<-paste0('N_2_',recie_N_2$X)
recie_S_1$X<-paste0('S_1_',recie_S_1$X)
recie_S_2$X<-paste0('S_2_',recie_S_2$X)
recie_F_1$X<-paste0('F_1_',recie_F_1$X)
recie_F_2$X<-paste0('F_2_',recie_F_2$X)
recie_FS_1$X<-paste0('FS_1_',recie_FS_1$X)
recie_FS_2$X<-paste0('FS_2_',recie_FS_2$X)
rbind(recie_N_1,recie_N_2,recie_S_1,recie_S_2,recie_F_1,recie_F_2,recie_FS_1,recie_FS_2)->recie

#clean the data
rownames(recie)<-recie$X
recie<-recie[,-1]
rownames(send)<-send$X
send<-send[,-1]

inter<-intersect(rownames(recie),colnames(lung_combine_after_batch))
lung_combine_after_batch<-lung_combine_after_batch[,inter]
recie<-recie[inter,]
send<-send[inter,]

#incoperate to the seurat object
lung_combine_after_batch[['reciever']]<-CreateAssayObject(as.data.frame(t(recie[match(colnames(lung_combine_after_batch),rownames(recie)),])))
lung_combine_after_batch[['sender']]<-CreateAssayObject(as.data.frame(t(send[match(colnames(lung_combine_after_batch),rownames(send)),])))







