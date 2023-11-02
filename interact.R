rm(list = ls())
###############################4_protein_pearson#####################
setwd("E:\\Testis_IHC\\network\\108")
gene_2<-c("NSUN5","UTP6","CCNT1","PRPF4")
id_map<-read.delim("E:\\Testis_IHC\\network\\id_map.txt",header = T,check.names = F)
gene_lnc<-id_map[which(id_map[,3]=="lncRNA"),2]
gene_mi<-id_map[which(id_map[,3]=="miRNA"),2]
data_tumor<-read.delim("E:\\Testis_IHC\\network\\data_tumor.txt",header = T,sep = "\t")
gene_p1<-setdiff(rownames(data_tumor),gene_mi)
gene_p1<-setdiff(gene_p1,gene_lnc)
gene_lnc1<-intersect(gene_lnc,rownames(data_tumor))
gene_mi1<-intersect(gene_mi,rownames(data_tumor))
data_tumor_lnc<-data_tumor[match(gene_lnc1,rownames(data_tumor)),]
data_tumor_mi<-data_tumor[match(gene_mi1,rownames(data_tumor)),]
data_tumor_p<-data_tumor[match(gene_p1,rownames(data_tumor)),]

data_normal<-read.delim("E:\\Testis_IHC\\network\\data_normal11.txt",header = T,sep = "\t")
gene_lnc2<-intersect(gene_lnc,rownames(data_normal))
gene_mi2<-intersect(gene_mi,rownames(data_normal))
data_normal_lnc<-data_normal[match(gene_lnc2,rownames(data_normal)),]
data_normal_mi<-data_normal[match(gene_mi2,rownames(data_normal)),]
data_normal_p<-data_normal[match(gene_p1,rownames(data_normal)),]

data_exp <- data_tumor[match(gene_2,rownames(data_tumor)),]
p1 <- matrix(NA,nrow(data_tumor_p),4)
R<-matrix(NA,nrow(data_tumor_p),4)
for (j in 1:4) {
  x<-matrix(NA,nrow(data_tumor_p),1)
  y<- matrix(NA,nrow(data_tumor_p),1)
  for(i in 1:nrow(data_tumor_p)){
    p<-cor.test(as.numeric(data_tumor_p[i,]),as.numeric(data_exp[j,]),method = ("pearson"))
    x[i,1] <- p$p.value
    y[i,1] <- p$estimate
    
  }
  p1[,j] <- x
  R[,j]  <- y
  print(j)
}
p_p<-as.data.frame(p1)
colnames(p_p)<-gene_2
rownames(p_p)<-rownames(data_tumor_p)
R_p<-as.data.frame(R)
colnames(R_p)<-gene_2
rownames(R_p)<-rownames(data_tumor_p)


p1 <- matrix(NA,nrow(data_tumor_mi),4)
R<-matrix(NA,nrow(data_tumor_mi),4)
for (j in 1:4) {
  x<-matrix(NA,nrow(data_tumor_mi),1)
  y<- matrix(NA,nrow(data_tumor_mi),1)
  for(i in 1:nrow(data_tumor_mi)){
    p<-cor.test(as.numeric(data_tumor_mi[i,]),as.numeric(data_exp[j,]),method = ("pearson"))
    x[i,1] <- p$p.value
    y[i,1] <- p$estimate
    
  }
  p1[,j] <- x
  R[,j]  <- y
  print(j)
}
p_mi<-as.data.frame(p1)
colnames(p_mi)<-gene_2
rownames(p_mi)<-rownames(data_tumor_mi)
R_mi<-as.data.frame(R)
colnames(R_mi)<-gene_2
rownames(R_mi)<-rownames(data_tumor_mi)


p1 <- matrix(NA,nrow(data_tumor_lnc),4)
R<-matrix(NA,nrow(data_tumor_lnc),4)
for (j in 1:4) {
  x<-matrix(NA,nrow(data_tumor_lnc),1)
  y<- matrix(NA,nrow(data_tumor_lnc),1)
  for(i in 1:nrow(data_tumor_lnc)){
    p<-cor.test(as.numeric(data_tumor_lnc[i,]),as.numeric(data_exp[j,]),method = ("pearson"))
    x[i,1] <- p$p.value
    y[i,1] <- p$estimate
    
  }
  p1[,j] <- x
  R[,j]  <- y
  print(j)
}
p_lnc<-as.data.frame(p1)
colnames(p_lnc)<-gene_2
rownames(p_lnc)<-rownames(data_tumor_lnc)
R_lnc<-as.data.frame(R)
colnames(R_lnc)<-gene_2
rownames(R_lnc)<-rownames(data_tumor_lnc)

gene_p<-list()
for (i in 1:4) {
  p<-p.adjust(p_p[,i],method="BH",length(p_p[,i]))#fdr校正
  g1<-which(p<0.05)
  #g1<-which(p_p[,i]<0.05)
  gene_p[[i]]<-cbind(rownames(data_tumor_p)[g1],R_p[g1,i])
  #gene_p[[i]]<-gene_p[[i]][-which(gene_p[[i]][,2]==1),]
  gene_p[[i]]<-gene_p[[i]][which(gene_p[[i]][,2]>0.3),]
}
names(gene_p)<-gene_2

gene_lnc<-list()
for (i in 1:4) {
  #g1<-which(p_lnc[,i]<0.05)
  p<-p.adjust(p_lnc[,i],method="BH",length(p_lnc[,i]))#fdr校正
  g1<-which(p<0.05)
  gene_lnc[[i]]<-cbind(rownames(data_tumor_lnc)[g1],R_lnc[g1,i])
  gene_lnc[[i]]<-gene_lnc[[i]][which(gene_lnc[[i]][,2]>0.3),]
}
names(gene_lnc)<-gene_2

gene_mi<-list()
for (i in 1:4) {
  #g1<-which(p_mi[,i]<0.05)
  p<-p.adjust(p_mi[,i],method="fdr",length(p_mi[,i]))#fdr校正
  g1<-which(p<0.05)
  gene_mi[[i]]<-cbind(rownames(data_tumor_mi)[g1],R_mi[g1,i])
  gene_mi[[i]]<-gene_mi[[i]][which(abs(R_mi[g1,i])>0.3),]
  gene_mi[[i]]<-gene_mi[[i]][which(gene_mi[[i]][,2]<0),]
}
names(gene_mi)<-gene_2

gene_in_cancer<-list()
for (i in 1:4) {
  gene1<-rbind(gene_p[[i]],gene_mi[[i]],gene_lnc[[i]])
  gene_in_cancer[[i]]<-gene1[,1]
  write(gene_in_cancer[[i]],paste("gene",gene_2[i],"c_in.txt",sep = "_"))
}
names(gene_in_cancer)<-colnames(p_p)

gene_all<-NULL
for (i in 1:4) {
  g1<-gene_2[i]
  gene_in1<-paste("gene",g1,"c_in.txt",sep = "_")
  gene_in1<-read.table(gene_in1,header = F)
  gene_all<-rbind(gene_all,gene_in1)
}
gene_all_cancer<-unique(gene_all[,1])

#####################protein_inter#########################################
###
#string
data_inter_pp<-NULL
for (i in 1:4) {
  g1<-gene_2[i]
  d1<-paste(g1,"0.15.tsv",sep = "_")
  d2<-read.delim(d1,header = T,check.names = F)
  data_inter_pp<-rbind(data_inter_pp,d2)
}

gene_in_pp<-list()
for (i in 1:4) {
  g1<-gene_2[i]
  d1<-which(data_inter_pp[,1]==g1)
  d2<-which(data_inter_pp[,2]==g1)
  gene_in_pp[[i]]<-unique(c(data_inter_pp[d1,2],data_inter_pp[d2,1]))
}
gene_in_cancer_p<-list()
for (i in 1:4) {
  gene_in_cancer_p[[i]]<-intersect(gene_p[[i]][,1],gene_in_pp[[i]])
}

#biogrid
data_inter_pp<-NULL
for (i in 1:4) {
  g1<-gene_2[i]
  d1<-paste(g1,"PP.txt",sep = "_")
  d2<-read.delim(d1,header = T,check.names = F)
  data_inter_pp<-rbind(data_inter_pp,d2)
}

gene_in_pp<-list()
for (i in 1:4) {
  g1<-gene_2[i]
  d1<-which(data_inter_pp$`Official Symbol Interactor A`==g1)
  d2<-which(data_inter_pp$`Official Symbol Interactor B`==g1)
  gene_in_pp[[i]]<-unique(c(data_inter_pp[d1,9],data_inter_pp[d2,8]))
}
gene_in_cancer_p<-list()
for (i in 1:4) {
  gene_in_cancer_p[[i]]<-intersect(gene_p[[i]][,1],gene_in_pp[[i]])
}
gene_in_cancer_p1<-unique(c(gene_in_cancer_p[[1]],gene_in_cancer_p[[2]],
                            gene_in_cancer_p[[3]],gene_in_cancer_p[[4]]))
write(gene_in_cancer_p1,"gene_in_cancer_p.txt")


data_p<-data_tumor[match(gene_in_cancer_p1,rownames(data_tumor)),]
#pearson
p1 <- matrix(NA,nrow(data_p),nrow(data_p))
R<-matrix(NA,nrow(data_p),nrow(data_p))
for (j in 1:nrow(data_p)) {
  x<-matrix(NA,nrow(data_p),1)
  y<- matrix(NA,nrow(data_p),1)
  for(i in 1:nrow(data_p)){
    p<-cor.test(as.numeric(data_p[i,]),as.numeric(data_p[j,]),method = ("pearson"))
    x[i,1] <- p$p.value
    y[i,1] <- p$estimate
    
  }
  p1[,j] <- x
  R[,j]  <- y
  print(j)
}
pp_c<-as.data.frame(p1)
colnames(pp_c)<-rownames(data_p)
rownames(pp_c)<-rownames(data_p)
pR_c<-as.data.frame(R)
colnames(pR_c)<-rownames(data_p)
rownames(pR_c)<-rownames(data_p)

gene_p_c<-list()
for (i in 1:131) {
  g1<-which(pp_c[,i]<0.05)
  gene_p_c[[i]]<-cbind(rownames(data_p)[g1],pR_c[g1,i])
  gene_p_c[[i]]<-gene_p_c[[i]][which(gene_p_c[[i]][,2]>0.3),]
}
names(gene_p_c)<-rownames(data_p)
#inter
data_inter_p_c<-read.delim("cancer_all.tsv",header = T,sep = "\t")
###########cancer_net
gene_in_p<-list()
data_in_p<-list()
for (i in 1:131) {
  g1<-rownames(data_p)[i]
  d1<-which(data_inter_p_c[,1]==g1)
  d2<-which(data_inter_p_c[,2]==g1)
  data_in_p[[i]]<-data_inter_p_c[c(d1,d2),]
  gene_in_p[[i]]<-unique(c(data_in_p[[i]][,1],data_in_p[[i]][,2]))
  gene_in_p[[i]]<-setdiff(gene_in_p[[i]],g1)
}
names(gene_in_p)<-rownames(data_p)
names(data_in_p)<-rownames(data_p)
gene_in_cancer_pp<-list()
for (i in 1:131) {
  gene_in_cancer_pp[[i]]<-intersect(gene_in_p[[i]],gene_p_c[[i]])
}

net<-NULL
for (i in 1:131) {
  g1<-rownames(data_p)[i]
  gene_in1<-gene_in_cancer_pp[[i]]
  net<-rbind(net,cbind(rep(g1,length(unique(gene_in1))),unique(gene_in1)))
  print(i)
}
net_cancer3<-net
colnames(net_cancer3)<-colnames(net_cancer)[1:2]
net_cancer3<-rbind(net_cancer2[,1:2],net_cancer3)
net_cancer3<-unique(as.data.frame(t(apply(net_cancer3,1,sort))))

net<-NULL
for (i in 1:4) {
  g1<-gene_2[i]
  gene_in1<-gene_in_cancer_p[[i]]
  net<-rbind(net,cbind(rep(g1,length(unique(gene_in1))),unique(gene_in1)))
  print(i)
}
colnames(net_cancer4)<-colnames(net_cancer)[1:2]
net_cancer4<-rbind(net_cancer3[,1:2],net_cancer4)
net_cancer4<-unique(as.data.frame(t(apply(net_cancer4,1,sort))))

type_c<-rep("mRNA",nrow(net_cancer4))
for (i in 1:3) {
  type_c[which(net_cancer4[,2]==gene_net_cancer_lnc[i])]<-"lncRNA"
  type_c[which(net_cancer4[,1]==gene_net_cancer_lnc[i])]<-"lncRNA"
}
gene_in_cancer_p2<-setdiff(gene_in_cancer_p1,gene_2)
for (i in 1:131) {
  type_c[which(net_cancer4[,2]==gene_in_cancer_p2[i])]<-"protein"
  type_c[which(net_cancer4[,1]==gene_in_cancer_p2[i])]<-"protein"
}
net_cancer4<-cbind(net_cancer4,type_c)
write.table(net_cancer4,"net_cancer_final_p.txt",sep = "\t",quote = F,row.names = F,col.names = T)

###############################protein_RNA#####################
data_inter_rr_all<-read.delim("E:\\immunohistochemistryimages\\RNAInter\\Download_data_RR.txt",header = T,sep = "\t")
#data_inter_rr_all<-NULL
data_inter2<-matrix(0,1,13)
colnames(data_inter2)<-colnames(data_inter_rr_all)
for (i in 1:4) {
  g1<-gene_2[i]
  d1<-which(data_inter_rr_all$Interactor1.Symbol==g1)
  d2<-which(data_inter_rr_all$Interactor2.Symbol==g1)
  data_inter1<-data_inter_rr_all[c(d1,d2),]
  data_inter2<-rbind(data_inter2,data_inter1)
  print(i)
}
data_inter_rr<-data_inter2[-1,]
write.table(data_inter_rr,"data_inter_rr.txt",sep = "\t",quote = F,row.names = F,col.names = T)

data_inter_rp_all<-read.delim("E:\\immunohistochemistryimages\\RNAInter\\Download_data_RP.txt",header = T,sep = "\t")
data_inter2<-matrix(0,1,13)
colnames(data_inter2)<-colnames(data_inter_rp_all)
for (i in 1:4) {
  g1<-gene_2[i]
  d1<-which(data_inter_rp_all$Interactor1.Symbol==g1)
  d2<-which(data_inter_rp_all$Interactor2.Symbol==g1)
  data_inter1<-data_inter_rp_all[c(d1,d2),]
  data_inter2<-rbind(data_inter2,data_inter1)
  print(i)
}
data_inter_rp<-data_inter2[-1,]
write.table(data_inter_rp,"data_inter_rp.txt",sep = "\t",quote = F,row.names = F,col.names = T)

data_inter_rr<-read.delim("data_inter_rr.txt",header = T,sep = "\t")
data_inter_rp<-read.delim("data_inter_rp.txt",header = T,sep = "\t")
data_inter_rh<-read.delim("E:\\immunohistochemistryimages\\RNAInter\\Download_data_RH.txt",header = T,sep = "\t")
data_inter_rc<-read.delim("E:\\immunohistochemistryimages\\RNAInter\\Download_data_RC.txt",header = T,sep = "\t")
data_inter_rd<-read.delim("E:\\immunohistochemistryimages\\RNAInter\\Download_data_RD.txt",header = T,sep = "\t")
data_inter_all<-rbind(data_inter_rr,data_inter_rp,data_inter_rh,data_inter_rc,data_inter_rd)
data_inter2<-matrix(0,1,13)
colnames(data_inter2)<-colnames(data_inter_all)
for (i in 1:4) {
  g1<-gene_2[i]
  d1<-which(data_inter_all$Interactor1.Symbol==g1)
  d2<-which(data_inter_all$Interactor2.Symbol==g1)
  data_inter1<-data_inter_all[c(d1,d2),]
  data_inter2<-rbind(data_inter2,data_inter1)
  print(i)
}
data_inter<-matrix(0,1,13)
colnames(data_inter)<-colnames(data_inter_all)
for (i in 1:4) {
  g1<-gene_2[i]
  d1<-data_inter2[which(data_inter2$Interactor2.Symbol==g1),]
  d2<-d1[which(d1$Category2=="protein"|d1$Category2=="RBP"|d1$Category2=="TF"),]
  data_inter<-rbind(data_inter,d2)
  print(i)
}
data_inter<-data_inter[-1,]
write.table(data_inter,"data_inter_protein.txt",sep = "\t",quote = F,row.names = F,col.names = T)
data_inter_shiyan<-data_inter[-which(data_inter$weak=="N/A"),]
##inter
data_inter<-read.delim("data_inter_RNA.txt",header = T,sep = "\t")
g11<-data_inter$Interactor1.Symbol
g11<-gsub("hsa-miR-","MIR",g11)
g11<-gsub("hsa-let-","MIRLET",g11)
g11<-gsub("hsa-mir-","MIR",g11)
g11[which(data_inter$Category1=="miRNA")]<-toupper(g11[which(data_inter$Category1=="miRNA")])
g11<-gsub("-5P","",g11)
g11<-gsub("-3P","",g11)
data_inter$Interactor1.Symbol<-g11

gene_in_cancer_final<-list()
for (i in 1:4) {
  g1<-gene_2[i]
  gene1<-gene_in_cancer[[i]]
  gene2<-unique(data_inter$Interactor1.Symbol[which(data_inter$Interactor2.Symbol==g1)])
  gene_in_cancer_final[[i]]<-intersect(gene1,gene2)
  print(i)
}
gene_in_cancer_final1<-list()
for (i in 1:4) {
  g1<-gene_2[i]
  gene1<-gene_in_cancer[[i]]
  gene2<-unique(data_inter_shiyan$Interactor1.Symbol[which(data_inter_shiyan$Interactor2.Symbol==g1)])
  gene_in_cancer_final1[[i]]<-intersect(gene1,gene2)
  print(i)
}

net<-NULL
for (i in 1:4) {
  g1<-gene_2[i]
  gene_in1<-gene_in_cancer_final[[i]]
  net<-rbind(net,cbind(rep(g1,length(unique(gene_in1))),unique(gene_in1)))
  print(i)
}
net_cancer<-net

type_c<-rep("mRNA",nrow(net_cancer))
for (i in 1:nrow(net_cancer)) {
  g1<-net_cancer[i,2]
  l1<-data_inter[which(data_inter$Interactor1.Symbol==g1),3]
  type_c[i]<-unique(l1)
}
net_cancer<-cbind(net_cancer,type_c)
colnames(net_cancer)<-c("p1","p2","type")

#
gene_cancer<-unique(c(gene_in_cancer_final[[1]],gene_in_cancer_final[[2]],
                      gene_in_cancer_final[[3]],gene_in_cancer_final[[4]]))
gene_net_cancer_mi<-intersect(gene_cancer,gene_mi1)
gene_net_cancer_lnc<-intersect(gene_cancer,gene_lnc1)
gene_net_cancer<-c(gene_2,gene_cancer)
data_net_tumor<-data_tumor[match(gene_cancer,rownames(data_tumor)),]

p1 <- matrix(NA,nrow(data_net_tumor),nrow(data_net_tumor))
R<-matrix(NA,nrow(data_net_tumor),nrow(data_net_tumor))
for (j in 1:nrow(data_net_tumor)) {
  x<-matrix(NA,nrow(data_net_tumor),1)
  y<- matrix(NA,nrow(data_net_tumor),1)
  for(i in 1:nrow(data_net_tumor)){
    p<-cor.test(as.numeric(data_net_tumor[i,]),as.numeric(data_net_tumor[j,]),method = ("pearson"))
    x[i,1] <- p$p.value
    y[i,1] <- p$estimate
    
  }
  p1[,j] <- x
  R[,j]  <- y
  print(j)
}
p_c<-as.data.frame(p1)
colnames(p_c)<-rownames(data_net_tumor)
rownames(p_c)<-rownames(data_net_tumor)
R_c<-as.data.frame(R)
colnames(R_c)<-rownames(data_net_tumor)
rownames(R_c)<-rownames(data_net_tumor)

gene_c<-list()
for (i in 1:61) {
  g1<-which(p_c[,i]<0.05)
  gene_c[[i]]<-cbind(rownames(data_net_tumor)[g1],R_c[g1,i])
  gene_c[[i]]<-gene_c[[i]][which(gene_c[[i]][,2]>0.3),]
}
names(gene_c)<-rownames(data_net_tumor)

data_inter_rr_all<-read.delim("E:\\immunohistochemistryimages\\RNAInter\\Download_data_RR.txt",header = T,sep = "\t")
#data_inter_rr_all<-NULL
data_inter2<-matrix(0,1,13)
colnames(data_inter2)<-colnames(data_inter_rr_all)
for (i in 1:61) {
  g1<-gene_cancer[i]
  d1<-which(data_inter_rr_all$Interactor1.Symbol==g1)
  d2<-which(data_inter_rr_all$Interactor2.Symbol==g1)
  data_inter1<-data_inter_rr_all[c(d1,d2),]
  data_inter2<-rbind(data_inter2,data_inter1)
  print(i)
}
data_inter_rr_cancer<-data_inter2[-1,]

data_inter_rp_all<-read.delim("E:\\immunohistochemistryimages\\RNAInter\\Download_data_RP.txt",header = T,sep = "\t")
data_inter2<-matrix(0,1,13)
colnames(data_inter2)<-colnames(data_inter_rp_all)
for (i in 1:61) {
  g1<-gene_cancer[i]
  d1<-which(data_inter_rp_all$Interactor1.Symbol==g1)
  d2<-which(data_inter_rp_all$Interactor2.Symbol==g1)
  data_inter1<-data_inter_rp_all[c(d1,d2),]
  data_inter2<-rbind(data_inter2,data_inter1)
  print(i)
}
data_inter_rp_cancer<-data_inter2[-1,]

data_inter_cancer<-rbind(data_inter_rr_cancer,data_inter_rp_cancer,data_inter_rh,data_inter_rc,data_inter_rd)

data_inter2<-matrix(0,1,13)
colnames(data_inter2)<-colnames(data_inter_cancer)
for (i in 1:length(gene_cancer)) {
  g1<-gene_cancer[i]
  d1<-which(data_inter_cancer$Interactor1.Symbol==g1)
  data_inter1<-data_inter_cancer[d1,]
  n<-c()
  for (j in 1:length(gene_cancer)) {
    g2<-gene_cancer[j]
    d2<-which(data_inter1$Interactor2.Symbol==g2)
    n<-c(n,d2)
  }
  data_inter1<-data_inter1[n,]
  data_inter2<-rbind(data_inter2,data_inter1)
}
data_inter_cancer1<-data_inter2[-1,]
write.table(data_inter_cancer1,"data_inter_cancer.txt",sep = "\t",quote = F,row.names = F,col.names = T)

data_inter_cancer1<-read.delim("data_inter_cancer.txt",header = T,sep = "\t")
###########cancer_net
gene_in<-list()
data_in<-list()
for (i in 1:61) {
  g1<-rownames(data_net_tumor)[i]
  d1<-which(data_inter_cancer1$Interactor1.Symbol==g1)
  d2<-which(data_inter_cancer1$Interactor2.Symbol==g1)
  data_in[[i]]<-data_inter_cancer1[c(d1,d2),]
  gene_in[[i]]<-unique(c(data_in[[i]]$Interactor1.Symbol,data_in[[i]]$Interactor2.Symbol))
  gene_in[[i]]<-setdiff(gene_in[[i]],g1)
}
names(gene_in)<-rownames(data_net_tumor)
names(data_in)<-rownames(data_net_tumor)
gene_in_cancer<-list()
for (i in 1:61) {
  gene_in_cancer[[i]]<-intersect(gene_in[[i]],gene_c[[i]])
}

net<-NULL
for (i in 1:61) {
  g1<-rownames(data_net_tumor)[i]
  gene_in1<-gene_in_cancer[[i]]
  net<-rbind(net,cbind(rep(g1,length(unique(gene_in1))),unique(gene_in1)))
  print(i)
}
net_cancer1<-net

colnames(net_cancer1)<-colnames(net_cancer)[1:2]
net_cancer2<-rbind(net_cancer[,1:2],net_cancer1)
net_cancer2<-unique(as.data.frame(t(apply(net_cancer2,1,sort))))

type_c<-rep("mRNA",nrow(net_cancer2))
for (i in 1:3) {
  type_c[which(net_cancer2[,2]==gene_net_cancer_lnc[i])]<-"lncRNA"
  type_c[which(net_cancer2[,1]==gene_net_cancer_lnc[i])]<-"lncRNA"
}
net_cancer2<-cbind(net_cancer2,type_c)
colnames(net_cancer2)<-c("p1","p2","type")
write.table(net_cancer2,"net_cancer_final.txt",sep = "\t",quote = F,row.names = F,col.names = T)

