library(BiocStyle)
library(HPAanalyze)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
setwd("F:\\Testis_IHC")
tissue="Breast"
dir.create("img")
#获得HPA网站中该基因的xml文件
gene<-BRCA_final_limma_filt_0_0001[,1]
hpa_url<-NULL
for (i in 1:length(gene)) {
  g<-gene[i]
  hpa1<-hpaXmlGet(g)
  hpa_target_gene_fig_url<-hpaXmlTissueExpr(hpa1)
  if(length(hpa_target_gene_fig_url)!=0){
    hpa_target_gene_fig_url<-as.data.frame(hpa_target_gene_fig_url[[1]])
    if(nrow(hpa_target_gene_fig_url)==0){
      hpa_url<-hpa_url
    }
    else{
      geneid<-rep(g,nrow(hpa_target_gene_fig_url))
      url1<-cbind(geneid,hpa_target_gene_fig_url[,1:8],hpa_target_gene_fig_url$tissueDescription1,hpa_target_gene_fig_url$tissueDescription2)
      hpa_url<-rbind(hpa_url,url1)
    }
  }
  print(i)
}
write.csv(hpa_url,"IHC_result_BRCA.csv")
write.table(hpa_url,"hpa_url_brca.txt",sep = "\t",quote = F,row.names = F,col.names = T)

#IHC筛选
url_normal<-hpa_url[which(hpa_url$tissueDescription1=="Normal tissue, NOS"),]#5014
result1<-url_normal[which(url_normal$intensity=="Strong"),]
result2<-url_normal[which(url_normal$intensity=="Moderate"),]
result<-rbind(result1,result2)
result<-result[which(result$quantity==">75%"),]#3752
colnames(result)[10:14]<-c("tissueDescription1","tissueDescription2","tissueDescription3","tissueDescription4","tissueDescription5")
brca1<-result[which(result$tissueDescription2=="Breast"),]
brca2<-result[which(result$tissueDescription3=="Breast"),]
brca3<-result[which(result$tissueDescription4=="Breast"),]
brca4<-result[which(result$tissueDescription5=="Breast"),]
result_brca<-rbind(brca1,brca2,brca3,brca4)
brca_normal<-result_brca[which(result_brca$tissueDescription1=="Normal tissue, NOS"),]
brca_cancer<-result_brca[-which(result_brca$tissueDescription1=="Normal tissue, NOS"),]
gene<-intersect(unique(brca_cancer[,1]),unique(brca_normal[,1]))
n<-c()
for (i in 1:length(gene)) {
  g1<-gene[i]
  n1<-which(result_brca[,1]==g1)
  n<-c(n,n1)
}
hpa_url_brca<-result_brca[n,]
write.csv(hpa_url_brca,"IHC_brca.csv")

setwd("F:\\Testis_IHC")
hpa_url_brca<-as.data.frame(read.csv("IHC_brca.csv")[,-1])
hpa_brca_normal<-hpa_url_brca[which(hpa_url_brca$tissueDescription1=="Normal tissue, NOS"),]
hpa_brca_cancer<-hpa_url_brca[-which(hpa_url_brca$tissueDescription1=="Normal tissue, NOS"),]
gene<-unique(hpa_url_brca[,1])

#################################location#########################################
gene1 <- read_excel("E:/immunohistochemistry images/data/gene1.xlsx")
gene_all<-as.data.frame(gene1)
result<-NULL
for (i in 1:nrow(gene_all)) {
  g1<-gene_all[i,1]
  f1<-gene_all[i,2]
  l1<-gene_all[i,3]
  l<-strsplit(l1,split=",")
  a<-rbind(rep(g1,length(l[[1]])),rep(f1,length(l[[1]])),t(l[[1]]))
  result<-rbind(result,t(a))
  #print(i)
}
gene_in<-intersect(gene,unique(result[,1]))

n<-c()
for (i in 1:length(gene_in)) {
  g1<-gene_in[i]
  g2<-which(result[,1]==g1)
  n<-c(n,g2)
}
result1<-result[n,]
a<-which(result1[,3]=="Nuclear membrane"|result1[,3]=="Nucleoplasm"|result1[,3]=="Nucleoli"
         |result1[,3]== "Nuclear bodies"|result1[,3]=="Nuclear speckles"|result1[,3]=="Nucleoli rim"
         |result1[,3]== "Nucleoli fibrillar center"|result1[,3]== "Nuclear membrane")
result1[a,3]<-"Nuclear"
b<-which(result1[,3]=="Cytosol"|result1[,3]=="Cell Junctions"
         |result1[,3]=="Lipid droplets"|result1[,3]=="Actin filaments"|result1[,3]=="Intermediate filaments"
         |result1[,3]=="Focal adhesion sites"|result1[,3]=="Microtubules"|result1[,3]=="Cytokinetic bridge"
         |result1[,3]=="Midbody"|result1[,3]=="Mitotic chromosome"|result1[,3]=="Mitotic spindle"
         |result1[,3]=="Cytoplasmic bodies"|result1[,3]=="Midbody ring"|result1[,3]=="Rods & Rings"
         |result1[,3]=="Centriolar satellite"|result1[,3]=="Centrosome")
result1[b,3]<-"Cytoplasm"
d<-which(result1[,3]=="Endosomes"|result1[,3]=="Vesicles"|result1[,3]=="Peroxisomes"|result1[,3]=="Lysosomes")
result1[d,3]<-"Vesicles"
data_en<-result1[which(result1[,2]=="Enhanced"),]
data_ap<-result1[which(result1[,2]=="Approved"),]
data_su<-result1[which(result1[,2]=="Supported"),]

#################enhanced_gene IHC下载#####################################
gene_en<-unique(data_en[,1])
file<-c()
for (j in 1:296) {
  setwd("F:\\Testis_IHC\\Breast\\normal_local_enhanced")
  g1<-gene_en[j]
  dir.create(g1)
  setwd(g1)
  n1<-which(hpa_brca_normal[,1]==g1)
  file<-hpa_brca_normal[n1,]
  for (i in 1:nrow(file)) {
    file_url<-file$imageUrl[i]
    file_dir<-paste(file$geneid[i],file$patientId[i],file$tissueDescription1[i],".png",sep = "_")
    download.file(url = file_url,destfile = file_dir,mode = "wb")
  }
  print(paste(j,g1,sep = ","))
}

file<-c()
for (j in 1:296) {
  setwd("F:\\Testis_IHC\\Breast\\cancer_local_enhanced")
  g1<-gene_en[j]
  dir.create(g1)
  setwd(g1)
  n1<-which(hpa_brca_cancer[,1]==g1)
  file<-hpa_brca_cancer[n1,]
  for (i in 1:nrow(file)) {
    file_url<-file$imageUrl[i]
    file_dir<-paste(file$geneid[i],file$patientId[i],file$tissueDescription1[i],".png",sep = "_")
    download.file(url = file_url,destfile = file_dir,mode = "wb")
  }
  print(paste(j,g1,sep = ","))
}
