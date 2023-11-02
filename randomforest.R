gene_2<-c("NSUN5","UTP6","CCNT1","PRPF4")
clinic3 <- read_excel("clinic3.xlsx")
clinic_n<-clinic3[which(clinic3$`Sample Type`=="Solid Tissue Normal"),]
clinic_c<-clinic3[which(clinic3$`Sample Type`=="Primary Tumor"),]
sample_normal<-clinic_n$`Sample Submitter ID`
sample_tumor<-clinic_c$`Sample Submitter ID`
data<-read.delim('C:/Users/Administrator/Desktop/汇报-2022下/结果/3/CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv',header = T,check.names = F)
rownames(data)<-data[,1]
data<-data[,-1]
data<-data[,-c(307:312)]
data_p2<-data[,seq(0,ncol(data),2)]#Unshared Log Ratio
colnames(data_p2)<-substr(colnames(data_p2),1,30)
data_p_n2<-data_p2[,match(sample_normal,colnames(data_p2))]
n<-match(sample_tumor,colnames(data_p2))
n<-n[-which(is.na(n))]
data_p_c2<-data_p2[,n]
data_all<-cbind(data_p_c2,data_p_n2)
##############################diff############################################
#foldchange
fc<-rep(0,nrow(data_all))
for(i in 1:nrow(data_all)){
  fc[i]<-mean(as.numeric(data_p_c2[i,]))/mean(as.numeric(data_p_n2[i,]));
}
fc_log<-log(fc,base=2)
#t-test
n<-c(1:nrow(data_all))
p.t<-rep(NA,nrow(data_all))
for ( j in n){  
  p.t[j] <-t.test(data_p_c2[j,],data_p_n2[j,])$p.value;  
}
fdr.t<-p.adjust(p.t,method="fdr",length(p.t))
re<-cbind(fc_log,fdr.t,data_all)
diff<-re[which(abs(re[,1])>1&re[,2]<0.05),]#求差异基因
write.table(diff,"diff_141.txt",col.names = T,row.names = T,sep ="\t",quote=F)
diff_up<-rownames(diff)[which(diff[,1]>1)]
diff_down<-rownames(diff)[which(diff[,1]<1)]
####################################protein########################################
###留一法
#all
data4<-as.data.frame(t(data_exp1))
feature<-c(rep(1,123),rep(0,18))
pre_prob<-c()
pre_label1<-NA
for (i in 1:141) {
  train<-data4[-i,]
  test<-data4[i,]
  group<-factor(feature[-i]) 
  colnames(train) <- make.names(colnames(train))
  mydata.rf<-randomForest(group~.,data=train,importance=TRUE, proximity=TRUE)
  pre_ran <- predict(mydata.rf, newdata=test, probability = TRUE)
  pre_prob <-rbind(pre_prob,pre_ran)##attr改变数据类型
  pre_label1[i]<-as.character(pre_ran[1])
}
pred <- prediction(pre_prob[, 1], feature)
perf <- performance(pred, "tpr", "fpr")
auc.tmp <- performance(pred, "auc")
auc.value <- auc.tmp@y.values
round(auc.value[[1]], digits = 4)

#equal
auc_1_18_all<-c()
for (j in 1:10) {
  a0<-sample(1:123,18,replace=F)
  a1<-c(a0,124:141)
  data5<-data4[a1,]
  feature1<-feature[a1]
  pre_prob<-c()
  pre_label1<-NA
  for (i in 1:36) {
    train<-data5[-i,]
    test<-data5[i,]
    group<-factor(feature1[-i]) 
    colnames(train) <- make.names(colnames(train))
    mydata.rf<-randomForest(group~.,data=train,importance=TRUE, proximity=TRUE)
    pre_ran <- predict(mydata.rf, newdata=test, probability = TRUE)
    pre_prob <-rbind(pre_prob,pre_ran)##attr改变数据类型
    pre_label1[i]<-as.character(pre_ran[1])
  }
  pred <- prediction(pre_prob[, 1], feature1)
  perf <- performance(pred, "tpr", "fpr")
  auc.tmp <- performance(pred, "auc")
  auc.value <- auc.tmp@y.values
  auc_1_18_all[j]<-round(auc.value[[1]], digits = 4)
}
####################################others########################################
#####diff
library(ROCR)
library(class)
library(e1071)
library(pROC)
library(randomForest)
gene_diff<-rownames(diff)
gene_diff1<-setdiff(gene_diff,gene_2)
gene_diff1<-gsub("-",".",gene_diff1)

###留一法
#all
auc_1_22<-c()
for (j in 1:100) {
  a1<-sample(1:3383,2)
  a2<-sample(1:3831,2)
  gene1<-gene_diff1[a1]
  gene2<-gene_ndiff1[a2]
  data8<-data6[,c(gene1,gene2)]
  feature<-c(rep(1,123),rep(0,18))
  pre_prob<-c()
  pre_label1<-NA
  for (i in 1:141) {
    train<-data8[-i,]
    test<-data8[i,]
    group<-factor(feature[-i]) 
    colnames(train) <- make.names(colnames(train))
    mydata.rf<-randomForest(group~.,data=train,importance=TRUE, proximity=TRUE)
    pre_ran <- predict(mydata.rf, newdata=test, probability = TRUE)
    pre_prob <-rbind(pre_prob,pre_ran)##attr???????????
    pre_label1[i]<-as.character(pre_ran[1])
  }
  pred <- prediction(pre_prob[, 1], feature)
  perf <- performance(pred, "tpr", "fpr")
  auc.tmp <- performance(pred, "auc")
  auc.value <- auc.tmp@y.values
  auc_1_22[j]<-round(auc.value[[1]], digits = 4)
  print(j)
}

#equal
auc_1_18<-c()
for (j in 1:100) {
  a1<-sample(1:3383,2)
  a2<-sample(1:3831,2)
  gene1<-gene_diff1[a1]
  gene2<-gene_ndiff1[a2]
  data8<-data6[,c(gene1,gene2)]
  a0<-sample(1:123,18,replace=F)
  a1<-c(a0,124:141)
  data5<-data8[a1,]
  feature1<-feature[a1]
  pre_prob<-c()
  pre_label1<-NA
  for (i in 1:36) {
    train<-data5[-i,]
    test<-data5[i,]
    group<-factor(feature1[-i]) 
    colnames(train) <- make.names(colnames(train))
    mydata.rf<-randomForest(group~.,data=train,importance=TRUE, proximity=TRUE)
    pre_ran <- predict(mydata.rf, newdata=test, probability = TRUE)
    pre_prob <-rbind(pre_prob,pre_ran)##attr???????????
    pre_label1[i]<-as.character(pre_ran[1])
  }
  pred <- prediction(pre_prob[, 1], feature1)
  perf <- performance(pred, "tpr", "fpr")
  auc.tmp <- performance(pred, "auc")
  auc.value <- auc.tmp@y.values
  auc_1_18[j]<-round(auc.value[[1]], digits = 4)
  print(j)
}