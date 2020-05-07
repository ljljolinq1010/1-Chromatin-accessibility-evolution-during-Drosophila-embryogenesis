library("data.table")
library("RColorBrewer") 
library("plyr")
library("qvalue")
library("topGO")
library("ROCR")
library("tseries")
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3")) 

setwd("/Users/jliu5/Documents/PhD_thesis/final_phd/Chromatin-accessibility-evolution-during-Drosophila-embryogenesis/data")

#####** peak summary **#####
## all peaks
# dmel
dmelT1Peak<-fread("bedFiles/all_peaks/dmel24_IDR_0.05.bed")
dmelT2Peak<-fread("bedFiles/all_peaks/dmel68_IDR_0.05.bed")
dmelT3Peak<-fread("bedFiles/all_peaks/dmel1012_IDR_0.05.bed")
dmelT4Peak<-fread("bedFiles/all_peaks/dmel1416_IDR_0.05.bed")
dmelT5Peak<-fread("bedFiles/all_peaks/dmel1820_IDR_0.05.bed")
# dvir
dvirT1Peak<-fread("bedFiles/all_peaks/dvir25_IDR_0.05.bed")
dvirT2Peak<-fread("bedFiles/all_peaks/dvir710_IDR_0.05.bed")
dvirT3Peak<-fread("bedFiles/all_peaks/dvir1417_IDR_0.05.bed")
dvirT4Peak<-fread("bedFiles/all_peaks/dvir1922_IDR_0.05.bed")
dvirT5Peak<-fread("bedFiles/all_peaks/dvir2528_IDR_0.05.bed")
# plot
par(mfrow=c(1,1))
par(mar=c(4, 6, 4, 2) + 0.1)
bp<-barplot(c(nrow(dvirT1Peak),nrow(dvirT2Peak),nrow(dvirT3Peak),nrow(dvirT4Peak),nrow(dvirT5Peak)),yaxt='n',ylim=c(-25000,25000),main="Peak counts",col=c(pal[3],pal[7],pal[4],pal[5],pal[2]))
barplot(c(-nrow(dmelT1Peak),-nrow(dmelT2Peak),-nrow(dmelT3Peak),-nrow(dmelT4Peak),-nrow(dmelT5Peak)),yaxt='n',add=T,col=c(pal[3],pal[7],pal[4],pal[5],pal[2]))
axis(2, at=seq(-25000,25000,by = 5000),labels=c(seq(25000,5000,by = -5000),seq(0,25000,by = 5000)), las=2)
text(x=bp,y=-25100,cex.lab=1.2,adj = 1,labels = c("TP1", "TP2","TP3", "TP4","TP5"),xpd = TRUE)

## enhancer
# dmel
dmelT1Enh<-fread("bedFiles/enhancers/all_enhancers/dmel24_IDR_0.05_enhancer.bed")
dmelT2Enh<-fread("bedFiles/enhancers/all_enhancers/dmel68_IDR_0.05_enhancer.bed")
dmelT3Enh<-fread("bedFiles/enhancers/all_enhancers/dmel1012_IDR_0.05_enhancer.bed")
dmelT4Enh<-fread("bedFiles/enhancers/all_enhancers/dmel1416_IDR_0.05_enhancer.bed")
dmelT5Enh<-fread("bedFiles/enhancers/all_enhancers/dmel1820_IDR_0.05_enhancer.bed")
# dvir
dvirT1Enh<-fread("bedFiles/enhancers/all_enhancers/dvir25_IDR_0.05_enhancer.bed")
dvirT2Enh<-fread("bedFiles/enhancers/all_enhancers/dvir710_IDR_0.05_enhancer.bed")
dvirT3Enh<-fread("bedFiles/enhancers/all_enhancers/dvir1417_IDR_0.05_enhancer.bed")
dvirT4Enh<-fread("bedFiles/enhancers/all_enhancers/dvir1922_IDR_0.05_enhancer.bed")
dvirT5Enh<-fread("bedFiles/enhancers/all_enhancers/dvir2528_IDR_0.05_enhancer.bed")
# plot
par(mfrow=c(1,1))
par(mar=c(4, 6, 4, 2) + 0.1)
bp<-barplot(c(nrow(dvirT1Enh),nrow(dvirT2Enh),nrow(dvirT3Enh),nrow(dvirT4Enh),nrow(dvirT5Enh)),yaxt='n',ylim=c(-15000,15000),main="Enhancer counts",col=c(pal[3],pal[7],pal[4],pal[5],pal[2]))
barplot(c(-nrow(dmelT1Enh),-nrow(dmelT2Enh),-nrow(dmelT3Enh),-nrow(dmelT4Enh),-nrow(dmelT5Enh)),yaxt='n',add=T,col=c(pal[3],pal[7],pal[4],pal[5],pal[2]))
axis(2, at=seq(-15000,15000,by = 5000),labels=c(seq(15000,5000,by = -5000),seq(0,15000,by = 5000)), las=2)
text(x=bp,y=-15100,cex.lab=1.2,adj = 1,labels = c("TP1", "TP2","TP3", "TP4","TP5"),xpd = TRUE)


#####** stage specific enhancer turnover **#####
#####* alignment rate *#####
## satge specific enhancers
dmelT1EnhSpec<-fread("bedFiles/enhancers/stage_specific_enhancers/dmel24_IDR_0.05_specific_enhancer.bed")
dmelT2EnhSpec<-fread("bedFiles/enhancers/stage_specific_enhancers/dmel68_IDR_0.05_specific_enhancer.bed")
dmelT3EnhSpec<-fread("bedFiles/enhancers/stage_specific_enhancers/dmel1012_IDR_0.05_specific_enhancer.bed")
dmelT4EnhSpec<-fread("bedFiles/enhancers/stage_specific_enhancers/dmel1416_IDR_0.05_specific_enhancer.bed")
dmelT5EnhSpec<-fread("bedFiles/enhancers/stage_specific_enhancers/dmel1820_IDR_0.05_specific_enhancer.bed")
## satge specific enhancers' ortholog regions 
dmelT1SpeciEnhAlign<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_ortholog_annotation/dmel24_IDR_0.05_specific_enhancer_dm6_droVir3.bed")
dmelT2SpeciEnhAlign<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_ortholog_annotation/dmel68_IDR_0.05_specific_enhancer_dm6_droVir3.bed")
dmelT3SpeciEnhAlign<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_ortholog_annotation/dmel1012_IDR_0.05_specific_enhancer_dm6_droVir3.bed")
dmelT4SpeciEnhAlign<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_ortholog_annotation/dmel1416_IDR_0.05_specific_enhancer_dm6_droVir3.bed")
dmelT5SpeciEnhAlign<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_ortholog_annotation/dmel1820_IDR_0.05_specific_enhancer_dm6_droVir3.bed")

## satge specific enhancers
dvirT1EnhSpec<-fread("bedFiles/enhancers/stage_specific_enhancers/dvir25_IDR_0.05_specific_enhancer.bed")
dvirT2EnhSpec<-fread("bedFiles/enhancers/stage_specific_enhancers/dvir710_IDR_0.05_specific_enhancer.bed")
dvirT3EnhSpec<-fread("bedFiles/enhancers/stage_specific_enhancers/dvir1417_IDR_0.05_specific_enhancer.bed")
dvirT4EnhSpec<-fread("bedFiles/enhancers/stage_specific_enhancers/dvir1922_IDR_0.05_specific_enhancer.bed")
dvirT5EnhSpec<-fread("bedFiles/enhancers/stage_specific_enhancers/dvir2528_IDR_0.05_specific_enhancer.bed")
## satge specific enhancers' ortholog regions 
dvirT1SpeciEnhAlign<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_ortholog_annotation/dvir25_IDR_0.05_specific_enhancer_droVir3_dm6.bed")
dvirT2SpeciEnhAlign<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_ortholog_annotation/dvir710_IDR_0.05_specific_enhancer_droVir3_dm6.bed")
dvirT3SpeciEnhAlign<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_ortholog_annotation/dvir1417_IDR_0.05_specific_enhancer_droVir3_dm6.bed")
dvirT4SpeciEnhAlign<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_ortholog_annotation/dvir1922_IDR_0.05_specific_enhancer_droVir3_dm6.bed")
dvirT5SpeciEnhAlign<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_ortholog_annotation/dvir2528_IDR_0.05_specific_enhancer_droVir3_dm6.bed")

## dmel
par(mfrow=c(1,1))
par(mar=c(4, 6, 4, 2) + 0.1)
uniqAlignEnh<-c(nrow(dmelT1SpeciEnhAlign),nrow(dmelT2SpeciEnhAlign),nrow(dmelT3SpeciEnhAlign),nrow(dmelT4SpeciEnhAlign),nrow(dmelT5SpeciEnhAlign))
allEnh<-c(nrow(dmelT1EnhSpec),nrow(dmelT2EnhSpec),nrow(dmelT3EnhSpec),nrow(dmelT4EnhSpec),nrow(dmelT5EnhSpec))
dmelEnh<-data.frame(uniqAlignEnh, allEnh)   

bp<-barplot(dmelEnh$uniqAlignEnh/dmelEnh$allEnh,ylim=c(0,1),
            main="Proportion of uniquely aligned enhancers\n(D.melanogaster to D.virilis)",col=c(pal[3],pal[7],pal[4],pal[5],pal[2]))
text(x=bp,y=-0.1,cex.lab=1.2,adj = 1,labels = c("TP1", "TP2","TP3", "TP4","TP5"),xpd = TRUE)

## dvir
par(mfrow=c(1,1))
par(mar=c(4, 6, 4, 2) + 0.1)
uniqAlignEnh<-c(nrow(dvirT1SpeciEnhAlign),nrow(dvirT2SpeciEnhAlign),nrow(dvirT3SpeciEnhAlign),nrow(dvirT4SpeciEnhAlign),nrow(dvirT5SpeciEnhAlign))
allEnh<-c(nrow(dvirT1EnhSpec),nrow(dvirT2EnhSpec),nrow(dvirT3EnhSpec),nrow(dvirT4EnhSpec),nrow(dvirT5EnhSpec))
dvirEnh<-data.frame(uniqAlignEnh, allEnh)   

bp<-barplot(dvirEnh$uniqAlignEnh/dvirEnh$allEnh,
            main="Proportion of uniquely aligned enhancers\n(D.virilis to D.melanogaster)",col=c(pal[3],pal[7],pal[4],pal[5],pal[2]))
text(x=bp,y=-0.1,cex.lab=1.2,adj = 1,labels = c("TP1", "TP2","TP3", "TP4","TP5"),xpd = TRUE)

#####* proportion of conserved enhancers *#####
#####  strict conservation (at least one bp overlap between two ortholog enhancers) #####
dmelT1SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_overlap_annotation/dmel24_IDR_0.05_specific_enhancer_overlap_dvir.bed")
dmelT2SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_overlap_annotation/dmel68_IDR_0.05_specific_enhancer_overlap_dvir.bed")
dmelT3SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_overlap_annotation/dmel1012_IDR_0.05_specific_enhancer_overlap_dvir.bed")
dmelT4SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_overlap_annotation/dmel1416_IDR_0.05_specific_enhancer_overlap_dvir.bed")
dmelT5SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_overlap_annotation/dmel1820_IDR_0.05_specific_enhancer_overlap_dvir.bed")

## find the number of overlapped peaks
n=0
consEnh<-c(nrow(dmelT1SpeciEnhAnno[dmelT1SpeciEnhAnno$V4>n,]),nrow(dmelT2SpeciEnhAnno[dmelT2SpeciEnhAnno$V4>n,]),nrow(dmelT3SpeciEnhAnno[dmelT3SpeciEnhAnno$V4>n,]),
           nrow(dmelT4SpeciEnhAnno[dmelT4SpeciEnhAnno$V4>n,]),nrow(dmelT5SpeciEnhAnno[dmelT5SpeciEnhAnno$V4>n,]))
## calculate jaccard index
allEnh<-c(nrow(dmelT1SpeciEnhAnno)+nrow(dvirT1SpeciEnhAlign),nrow(dmelT2SpeciEnhAnno)+nrow(dvirT2SpeciEnhAlign),nrow(dmelT3SpeciEnhAnno)+nrow(dvirT3SpeciEnhAlign),
          nrow(dmelT4SpeciEnhAnno)+nrow(dvirT4SpeciEnhAlign),nrow(dmelT5SpeciEnhAnno)+nrow(dvirT5SpeciEnhAlign))
allEnh<-allEnh-consEnh
dmelEnh<-data.frame(consEnh, allEnh)   

## plot
par(mfrow=c(1,1))
par(mar=c(5, 5, 4, 2) + 0.1)
bplot<-barplot(dmelEnh$consEnh/dmelEnh$allEnh, ylim=c(0,0.1),main="Strict conservation",ylab = "Proportion of conserved enhancers",col=c(pal[3],pal[7],pal[4],pal[5],pal[2]),cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
text(x=bplot,y=-0.01,cex=1.5,adj = 1,labels = c("TP1", "TP2","TP3","TP4","TP5"),xpd = TRUE)

##### relax conservation (the distance between two peaks smaller than 1kb) #####
dmelT1SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_distance_annotation/dmel24_IDR_0.05_specific_enhancer_distance_dvir.bed")
dmelT1SpeciEnhAnno$V4 <- ifelse(dmelT1SpeciEnhAnno$V7<1000 , 1, 0)
dmelT1SpeciEnhAnno<-dmelT1SpeciEnhAnno[,c(1:4)]
dmelT1SpeciEnhAnno<-unique(dmelT1SpeciEnhAnno)

dmelT2SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_distance_annotation/dmel68_IDR_0.05_specific_enhancer_distance_dvir.bed")
dmelT2SpeciEnhAnno$V4 <- ifelse(dmelT2SpeciEnhAnno$V7<1000 , 1, 0)
dmelT2SpeciEnhAnno<-dmelT2SpeciEnhAnno[,c(1:4)]
dmelT2SpeciEnhAnno<-unique(dmelT2SpeciEnhAnno)

dmelT3SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_distance_annotation/dmel1012_IDR_0.05_specific_enhancer_distance_dvir.bed")
dmelT3SpeciEnhAnno$V4 <- ifelse(dmelT3SpeciEnhAnno$V7<1000 , 1, 0)
dmelT3SpeciEnhAnno<-dmelT3SpeciEnhAnno[,c(1:4)]
dmelT3SpeciEnhAnno<-unique(dmelT3SpeciEnhAnno)

dmelT4SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_distance_annotation/dmel1416_IDR_0.05_specific_enhancer_distance_dvir.bed")
dmelT4SpeciEnhAnno$V4 <- ifelse(dmelT4SpeciEnhAnno$V7<1000 , 1, 0)
dmelT4SpeciEnhAnno<-dmelT4SpeciEnhAnno[,c(1:4)]
dmelT4SpeciEnhAnno<-unique(dmelT4SpeciEnhAnno)

dmelT5SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_distance_annotation/dmel1820_IDR_0.05_specific_enhancer_distance_dvir.bed")
dmelT5SpeciEnhAnno$V4 <- ifelse(dmelT5SpeciEnhAnno$V7<1000 , 1, 0)
dmelT5SpeciEnhAnno<-dmelT5SpeciEnhAnno[,c(1:4)]
dmelT5SpeciEnhAnno<-unique(dmelT5SpeciEnhAnno)

## find the number of overlapped peaks
consEnh<-c(nrow(dmelT1SpeciEnhAnno[dmelT1SpeciEnhAnno$V4==1,]),nrow(dmelT2SpeciEnhAnno[dmelT2SpeciEnhAnno$V4==1,]),nrow(dmelT3SpeciEnhAnno[dmelT3SpeciEnhAnno$V4==1,]),
           nrow(dmelT4SpeciEnhAnno[dmelT4SpeciEnhAnno$V4==1,]),nrow(dmelT5SpeciEnhAnno[dmelT5SpeciEnhAnno$V4==1,]))

## calculate jaccard index
allEnh<-c(nrow(dmelT1SpeciEnhAnno)+nrow(dvirT1SpeciEnhAlign),nrow(dmelT2SpeciEnhAnno)+nrow(dvirT2SpeciEnhAlign),nrow(dmelT3SpeciEnhAnno)+nrow(dvirT3SpeciEnhAlign),
          nrow(dmelT4SpeciEnhAnno)+nrow(dvirT4SpeciEnhAlign),nrow(dmelT5SpeciEnhAnno)+nrow(dvirT5SpeciEnhAlign))
allEnh<-allEnh-consEnh
dmelEnh<-data.frame(consEnh, allEnh)

## plot
par(mfrow=c(1,1))
par(mar=c(5, 5, 4, 2) + 0.1)
bplot<-barplot(dmelEnh$consEnh/dmelEnh$allEnh, ylim=c(0,0.1),main="Relaxed conservation",ylab = "Proportion of conserved enhancers",col=c(pal[3],pal[7],pal[4],pal[5],pal[2]),cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
text(x=bplot,y=-0.01,cex=1.5,adj = 1,labels = c("TP1", "TP2","TP3","TP4","TP5"),xpd = TRUE)


#####** positive selection **#####
#####* ROC analysis for SVM models *#####
##### function used for ROC analysis #####
rocPrArea<-function(positive, negative) {

  positive$state<-rep(1,tiems=nrow(positive))
  negative$state<-rep(0,tiems=nrow(negative))
  aucResult<-c()
  pr_plotValueX<-list()
  pr_plotValueY<-list()

  for (i in 2:5) {
    pos<-cbind(positive[[i]],positive[,6])
    neg<-cbind(negative[[i]],negative[,6])
    tempData<-rbind(pos,neg)
    pred <- prediction(tempData$V1,tempData$state) 
    perf <- performance( pred, "tpr", "fpr" )
    
    auc_temp <-performance( pred, measure = "auc")
    aucResult[i-1]<-round(unlist(slot(auc_temp, "y.values")),digits = 3)

    pr_temp<-data.frame(unlist(perf@x.values),unlist(perf@y.values))
    names(pr_temp)<-c("x","y")
    pr_temp<-na.omit(pr_temp)
    pr_plotValueX[[i-1]]<-pr_temp$x
    pr_plotValueY[[i-1]]<-pr_temp$y
  } 
  return(list(pr_plotValueX,pr_plotValueY,aucResult))
}

##### take T1 enhancers as an example #####
## T1 model cross validation 
cv<-fread("svmModels/dmel24_specific_enhancer_gkmtrain.cvpred.txt")
crossV<-cv[,c(2,3)]
## T1 enhancers prediction based on T2 model 
pos2<-fread("svmModels/dmel68_model_gkmpredict_dmel24_specific_enhancer_posSet.txt")
neg2<-fread("svmModels/dmel68_model_gkmpredict_dmel24_specific_enhancer_negSet.txt")
## T1 enhancers prediction based on T3 model 
pos3<-fread("svmModels/dmel1012_model_gkmpredict_dmel24_specific_enhancer_posSet.txt")
neg3<-fread("svmModels/dmel1012_model_gkmpredict_dmel24_specific_enhancer_negSet.txt")
## T1 enhancers prediction based on T4 model 
pos4<-fread("svmModels/dmel1416_model_gkmpredict_dmel24_specific_enhancer_posSet.txt")
neg4<-fread("svmModels/dmel1416_model_gkmpredict_dmel24_specific_enhancer_negSet.txt")
## T1 enhancers prediction based on T5 model 
pos5<-fread("svmModels/dmel1820_model_gkmpredict_dmel24_specific_enhancer_posSet.txt")
neg5<-fread("svmModels/dmel1820_model_gkmpredict_dmel24_specific_enhancer_negSet.txt")

## calculate AUC
# T1 model cross validation 
pred <- prediction(crossV$V2, crossV$V3) 
perf <- performance( pred, "tpr", "fpr" )
auc_temp <-performance( pred, measure = "auc")
aucT1<-unlist(slot(auc_temp, "y.values"))
# T1 enhancers prediction based on other models
pos<-cbind(pos2,pos3$V2,pos4$V2,pos5$V2)
neg<-cbind(neg2,neg3$V2,neg4$V2,neg5$V2)
stage1<-rocPrArea(pos,neg)
aucOthers<-list(stage1[[1]],stage1[[2]],stage1[[3]])

## plot
plotColors<-c(pal[3],pal[7],pal[4],pal[5],pal[2])
plotLegend<-c("T1 model", "T2 model", "T3 model", "T4 model", "T5 model")

par(mar=c(5,5,2,2))
plot(unlist(perf@x.values),unlist(perf@y.values),lwd = 3,cex=1.5,col=plotColors[1],main="D.melanogaster T1 enhancers",xlab="False positive rate",ylab="Ture positive rate",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
legend(0.2,0.5,paste(plotLegend[1]," ","(AUC = ",round(aucT1,digits = 3),")",sep=""),lwd = 3,bty="n",cex=1.4,col = plotColors[1]) 

for (i in 1:4) {
  lines(unlist(aucOthers[[1]][i]),unlist(aucOthers[[2]][i]),col=plotColors[i+1],lwd=3)
}

legendY<-c(0.4,0.3,0.2,0.1)
for (i in 1:4) {
  legendName<-paste(plotLegend[i+1]," ","(AUC = ",round(aucOthers[[3]][i],digits = 3),")",sep="")
  legend(0.2,legendY[i],legendName,lwd = 3,bty="n",cex=1.4,col = plotColors[[i+1]]) 
  
}


#####* deltaSVM summary *#####
##### function used for deltaSVM analysis #####
dataMod<-function(deltaSVM) {
  ## change the first column into bed format
  deltaSVM<-as.data.frame(deltaSVM)
  names(deltaSVM)<-c("ID","deltaSVM","varNumb","pValue")
  splString<-strsplit(deltaSVM$ID,"_",fixed=TRUE)
  splString<-data.frame(unlist(splString))
  ID.bed<-matrix(splString$unlist.splString., ncol=5, byrow=TRUE)
  ID.bed<-ID.bed[,c(1:3)]
  ID.bed[,2]<-as.numeric(ID.bed[,2])-1
  deltaSVM<-cbind(ID.bed, deltaSVM[,c(2:4)])
  colnames(deltaSVM)[c(1:3)]<-c("chr","start","end")
  deltaSVM$start<-as.numeric(as.character(deltaSVM$start))
  deltaSVM$end<-as.numeric(as.character(deltaSVM$end))
  
  ## qvalue
  deltaSVM<-deltaSVM[order(deltaSVM$pValue),]
  deltaSVM$qValue<-qvalue(deltaSVM$pValue,pi0=1)$qvalues
  return(deltaSVM)
}

dmelT1Enh_deltaSVM<-fread("deltaSVM/dmel24_specific_enhancer_dm3_deltaSVM_highertailTest.txt")
dmelT1Enh_deltaSVM<-dataMod(dmelT1Enh_deltaSVM)

dmelT2Enh_deltaSVM<-fread("deltaSVM/dmel68_specific_enhancer_dm3_deltaSVM_highertailTest.txt")
dmelT2Enh_deltaSVM<-dataMod(dmelT2Enh_deltaSVM)

dmelT3Enh_deltaSVM<-fread("deltaSVM/dmel1012_specific_enhancer_dm3_deltaSVM_highertailTest.txt")
dmelT3Enh_deltaSVM<-dataMod(dmelT3Enh_deltaSVM)

dmelT4Enh_deltaSVM<-fread("deltaSVM/dmel1416_specific_enhancer_dm3_deltaSVM_highertailTest.txt")
dmelT4Enh_deltaSVM<-dataMod(dmelT4Enh_deltaSVM)

dmelT5Enh_deltaSVM<-fread("deltaSVM/dmel1820_specific_enhancer_dm3_deltaSVM_highertailTest.txt")
dmelT5Enh_deltaSVM<-dataMod(dmelT5Enh_deltaSVM)

## plot deltaSVM and pvalue
par(mfrow=c(1,2))
par(mar=c(7,5,4,2))
hist(dmelT1Enh_deltaSVM$deltaSVM,breaks = 50,main="D.melanogaster T1",xlab="deltaSVM",xlim=c(-50,50),
     cex.lab=2,cex.axis=2,cex.main=2,col=pal[3])
hist(dmelT1Enh_deltaSVM$pValue,breaks = 50,main="D.melanogaster T1",xlab="Pvalue",xlim=c(0,1),
     cex.lab=2,cex.axis=2,cex.main=2,col=pal[3])

hist(dmelT2Enh_deltaSVM$deltaSVM,breaks = 50,main="D.melanogaster T2",xlab="deltaSVM",xlim=c(-50,50),
     cex.lab=2,cex.axis=2,cex.main=2,col=pal[7])
hist(dmelT2Enh_deltaSVM$pValue,breaks = 50,main="D.melanogaster T2",xlab="Pvalue",xlim=c(0,1),
     cex.lab=2,cex.axis=2,cex.main=2,col=pal[7])

hist(dmelT3Enh_deltaSVM$deltaSVM,breaks = 50,main="D.melanogaster T3",xlab="deltaSVM",xlim=c(-50,50),
     cex.lab=2,cex.axis=2,cex.main=2,col=pal[4])
hist(dmelT3Enh_deltaSVM$pValue,breaks = 50,main="D.melanogaster T3",xlab="Pvalue",xlim=c(0,1),
     cex.lab=2,cex.axis=2,cex.main=2,col=pal[4])

hist(dmelT4Enh_deltaSVM$deltaSVM,breaks = 50,main="D.melanogaster T4",xlab="deltaSVM",xlim=c(-50,50),
     cex.lab=2,cex.axis=2,cex.main=2,col=pal[5])
hist(dmelT4Enh_deltaSVM$pValue,breaks = 50,main="D.melanogaster T4",xlab="Pvalue",xlim=c(0,1),
     cex.lab=2,cex.axis=2,cex.main=2,col=pal[5])

hist(dmelT5Enh_deltaSVM$deltaSVM,breaks = 50,main="D.melanogaster T5",xlab="deltaSVM",xlim=c(-50,50),
     cex.lab=2,cex.axis=2,cex.main=2,col=pal[2])
hist(dmelT5Enh_deltaSVM$pValue,breaks = 50,main="D.melanogaster T5",xlab="Pvalue",xlim=c(0,1),
     cex.lab=2,cex.axis=2,cex.main=2,col=pal[2])

#####* mk test validation *#####
##### take T1 enhancers as an example #####
dmelT1Enh_polys<-fread("mkTest/dmel24_IDR_0.05_specific_enhancer_polyNub_dm3.bed")
dmelT1Enh_polys$V1<-paste0(rep("chr",nrow(dmelT1Enh_polys)),dmelT1Enh_polys$V1)
dmelT1Enh_snvs<-fread("mkTest/dmel24_specific_enhancer_dm3_SNVs.txt")
colnames(dmelT1Enh_polys)<-c("chr","start","end","polyNumb")
colnames(dmelT1Enh_snvs)<-c("chr","start","end","snvNumb")
dmelT1Enh_polys_snvs<-merge(dmelT1Enh_polys,dmelT1Enh_snvs,by=c("chr","start","end"))

dmelT1Enh_polys_snvs$polyNumb<-as.numeric(as.character(dmelT1Enh_polys_snvs$polyNumb))
dmelT1Enh_polys_snvs<-na.omit(dmelT1Enh_polys_snvs)
dmelT1Enh_polys_snvs_deltaSVM<-merge(dmelT1Enh_polys_snvs,dmelT1Enh_deltaSVM,by=c("chr","start","end"))

## plot
data1<-subset(dmelT1Enh_polys_snvs_deltaSVM,dmelT1Enh_polys_snvs_deltaSVM$qValue<0.05)
data2<-subset(dmelT1Enh_polys_snvs_deltaSVM,dmelT1Enh_polys_snvs_deltaSVM$qValue>=0.05)

snvNumb<-c(sum(data1$snvNumb),sum(data2$snvNumb))
snpNumb<-c(sum(data1$polyNumb),sum(data2$polyNumb))

par(mfrow=c(1,1))
par(mar=c(9,6,4,2))
bp<-barplot(snvNumb/snpNumb,ylim=c(0,5),main="D.melanogaster T1",cex.lab=1.5,cex.main=1.5,ylab="# substitution / # polymorphism",col=pal[3],cex.axis=1.5)
text(x=bp,y=0-5/15,cex=1.5,srt = 45,adj = 1,labels = c("Posiitve sites", "Non-positive sites"),xpd = TRUE)
# fisher exact test
fTest1<-fisher.test(matrix(c(snvNumb[1],snpNumb[1],snvNumb[2],snpNumb[2]),ncol = 2,nrow = 2))
legend("topleft",cex=1.5,legend=paste("p=",signif(fTest1$p.value,3)),bty = 'n')

#####* proportion of positive selection in different stages *#####
dmelT1Enh_deltaSVM<-fread("deltaSVM/dmelT1Enh_deltaSVM_dm6.txt")
colnames(dmelT1Enh_deltaSVM)<-c("chr","start","end","deltaSVM","varNumb","qValue")
dmelT2Enh_deltaSVM<-fread("deltaSVM/dmelT2Enh_deltaSVM_dm6.txt")
colnames(dmelT2Enh_deltaSVM)<-c("chr","start","end","deltaSVM","varNumb","qValue")
dmelT3Enh_deltaSVM<-fread("deltaSVM/dmelT3Enh_deltaSVM_dm6.txt")
colnames(dmelT3Enh_deltaSVM)<-c("chr","start","end","deltaSVM","varNumb","qValue")
dmelT4Enh_deltaSVM<-fread("deltaSVM/dmelT4Enh_deltaSVM_dm6.txt")
colnames(dmelT4Enh_deltaSVM)<-c("chr","start","end","deltaSVM","varNumb","qValue")
dmelT5Enh_deltaSVM<-fread("deltaSVM/dmelT5Enh_deltaSVM_dm6.txt")
colnames(dmelT5Enh_deltaSVM)<-c("chr","start","end","deltaSVM","varNumb","qValue")

##### for strictly conserved enhancers (T1 specific enhancers in dmel are also T1 specific enhancers in dvir) #####
dmelT1SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_overlap_annotation/dmel24_IDR_0.05_specific_enhancer_overlap_dvir.bed")
dmelT2SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_overlap_annotation/dmel68_IDR_0.05_specific_enhancer_overlap_dvir.bed")
dmelT3SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_overlap_annotation/dmel1012_IDR_0.05_specific_enhancer_overlap_dvir.bed")
dmelT4SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_overlap_annotation/dmel1416_IDR_0.05_specific_enhancer_overlap_dvir.bed")
dmelT5SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_overlap_annotation/dmel1820_IDR_0.05_specific_enhancer_overlap_dvir.bed")

dmelT1SpeciEnhCons_strict<-dmelT1SpeciEnhAnno[dmelT1SpeciEnhAnno$V4>0]
colnames(dmelT1SpeciEnhCons_strict)<-c("chr","start","end","anno" )
dmelT1SpeciEnhCons_strict_deltaSVM<-merge(dmelT1SpeciEnhCons_strict,dmelT1Enh_deltaSVM,by=c("chr","start","end"))

dmelT2SpeciEnhCons_strict<-dmelT2SpeciEnhAnno[dmelT2SpeciEnhAnno$V4>0]
colnames(dmelT2SpeciEnhCons_strict)<-c("chr","start","end","anno" )
dmelT2SpeciEnhCons_strict_deltaSVM<-merge(dmelT2SpeciEnhCons_strict,dmelT2Enh_deltaSVM,by=c("chr","start","end"))

dmelT3SpeciEnhCons_strict<-dmelT3SpeciEnhAnno[dmelT3SpeciEnhAnno$V4>0]
colnames(dmelT3SpeciEnhCons_strict)<-c("chr","start","end","anno" )
dmelT3SpeciEnhCons_strict_deltaSVM<-merge(dmelT3SpeciEnhCons_strict,dmelT3Enh_deltaSVM,by=c("chr","start","end"))

dmelT4SpeciEnhCons_strict<-dmelT4SpeciEnhAnno[dmelT4SpeciEnhAnno$V4>0]
colnames(dmelT4SpeciEnhCons_strict)<-c("chr","start","end","anno" )
dmelT4SpeciEnhCons_strict_deltaSVM<-merge(dmelT4SpeciEnhCons_strict,dmelT4Enh_deltaSVM,by=c("chr","start","end"))

dmelT5SpeciEnhCons_strict<-dmelT5SpeciEnhAnno[dmelT5SpeciEnhAnno$V4>0]
colnames(dmelT5SpeciEnhCons_strict)<-c("chr","start","end","anno" )
dmelT5SpeciEnhCons_strict_deltaSVM<-merge(dmelT5SpeciEnhCons_strict,dmelT5Enh_deltaSVM,by=c("chr","start","end"))


posEnh<-c(nrow(dmelT1SpeciEnhCons_strict_deltaSVM[dmelT1SpeciEnhCons_strict_deltaSVM$qValue<0.05]),
          nrow(dmelT2SpeciEnhCons_strict_deltaSVM[dmelT2SpeciEnhCons_strict_deltaSVM$qValue<0.05]),
          nrow(dmelT3SpeciEnhCons_strict_deltaSVM[dmelT3SpeciEnhCons_strict_deltaSVM$qValue<0.05]),
          nrow(dmelT4SpeciEnhCons_strict_deltaSVM[dmelT4SpeciEnhCons_strict_deltaSVM$qValue<0.05]),
          nrow(dmelT5SpeciEnhCons_strict_deltaSVM[dmelT5SpeciEnhCons_strict_deltaSVM$qValue<0.05]))
allEnh<-c(nrow(dmelT1SpeciEnhCons_strict_deltaSVM),nrow(dmelT2SpeciEnhCons_strict_deltaSVM),nrow(dmelT3SpeciEnhCons_strict_deltaSVM),nrow(dmelT4SpeciEnhCons_strict_deltaSVM),nrow(dmelT5SpeciEnhCons_strict_deltaSVM))

par(mfrow=c(1,1))
par(mar=c(5, 5, 4, 2) + 0.1)
bplot<-barplot(posEnh/allEnh, ylim=c(0, 0.5), main="Strictly conserved enhancers",ylab = "Proportion of positive selection",col=c(pal[3],pal[7],pal[4],pal[5],pal[2]),
               cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
text(x=bplot,y=0-0.4*0.08,cex=1.5,adj = 1,labels = c("TP1", "TP2","TP3","TP4","TP5"),xpd = TRUE)
text(x=bplot,y=0+0.4*0.05,cex=1.2,labels = paste0(posEnh,"/",allEnh),xpd = TRUE)

##### for non-strictly conserved enhancers (T1 specific enhancers in dmel have overlapped enhancers in dvir) #####
dmelT1SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_overlap_annotation/dmel24_IDR_0.05_specific_enhancer_annoALlDvirSatges.bed")
dmelT2SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_overlap_annotation/dmel68_IDR_0.05_specific_enhancer_annoALlDvirSatges.bed")
dmelT3SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_overlap_annotation/dmel1012_IDR_0.05_specific_enhancer_annoALlDvirSatges.bed")
dmelT4SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_overlap_annotation/dmel1416_IDR_0.05_specific_enhancer_annoALlDvirSatges.bed")
dmelT5SpeciEnhAnno<-fread("bedFiles/enhancers/stage_specific_enhancers/enhancer_overlap_annotation/dmel1820_IDR_0.05_specific_enhancer_annoALlDvirSatges.bed")

dmelT1SpeciEnhCons<-dmelT1SpeciEnhAnno[dmelT1SpeciEnhAnno$V4>0|dmelT1SpeciEnhAnno$V5>0|dmelT1SpeciEnhAnno$V6>0|dmelT1SpeciEnhAnno$V7>0|dmelT1SpeciEnhAnno$V8>0,]
dmelT1SpeciEnhCons<-dmelT1SpeciEnhCons[,c(1:4)]
colnames(dmelT1SpeciEnhCons)<-c("chr","start","end","anno" )
dmelT1SpeciEnhCons_deltaSVM<-merge(dmelT1SpeciEnhCons,dmelT1Enh_deltaSVM,by=c("chr","start","end"))
dmelT1SpeciEnhCons_relax_deltaSVM<-dmelT1SpeciEnhCons_deltaSVM[!(dmelT1SpeciEnhCons_deltaSVM$chr%in%dmelT1SpeciEnhCons_strict_deltaSVM$chr&
                                                                   dmelT1SpeciEnhCons_deltaSVM$start%in%dmelT1SpeciEnhCons_strict_deltaSVM$start&
                                                                   dmelT1SpeciEnhCons_deltaSVM$end%in%dmelT1SpeciEnhCons_strict_deltaSVM$end),]


dmelT2SpeciEnhCons<-dmelT2SpeciEnhAnno[dmelT2SpeciEnhAnno$V4>0|dmelT2SpeciEnhAnno$V5>0|dmelT2SpeciEnhAnno$V6>0|dmelT2SpeciEnhAnno$V7>0|dmelT2SpeciEnhAnno$V8>0,]
dmelT2SpeciEnhCons<-dmelT2SpeciEnhCons[,c(1:4)]
colnames(dmelT2SpeciEnhCons)<-c("chr","start","end","anno" )
dmelT2SpeciEnhCons_deltaSVM<-merge(dmelT2SpeciEnhCons,dmelT2Enh_deltaSVM,by=c("chr","start","end"))
dmelT2SpeciEnhCons_relax_deltaSVM<-dmelT2SpeciEnhCons_deltaSVM[!(dmelT2SpeciEnhCons_deltaSVM$chr%in%dmelT2SpeciEnhCons_strict_deltaSVM$chr&
                                                                   dmelT2SpeciEnhCons_deltaSVM$start%in%dmelT2SpeciEnhCons_strict_deltaSVM$start&
                                                                   dmelT2SpeciEnhCons_deltaSVM$end%in%dmelT2SpeciEnhCons_strict_deltaSVM$end),]

dmelT3SpeciEnhCons<-dmelT3SpeciEnhAnno[dmelT3SpeciEnhAnno$V4>0|dmelT3SpeciEnhAnno$V5>0|dmelT3SpeciEnhAnno$V6>0|dmelT3SpeciEnhAnno$V7>0|dmelT3SpeciEnhAnno$V8>0,]
dmelT3SpeciEnhCons<-dmelT3SpeciEnhCons[,c(1:4)]
colnames(dmelT3SpeciEnhCons)<-c("chr","start","end","anno" )
dmelT3SpeciEnhCons_deltaSVM<-merge(dmelT3SpeciEnhCons,dmelT3Enh_deltaSVM,by=c("chr","start","end"))
dmelT3SpeciEnhCons_relax_deltaSVM<-dmelT3SpeciEnhCons_deltaSVM[!(dmelT3SpeciEnhCons_deltaSVM$chr%in%dmelT3SpeciEnhCons_strict_deltaSVM$chr&
                                                                   dmelT3SpeciEnhCons_deltaSVM$start%in%dmelT3SpeciEnhCons_strict_deltaSVM$start&
                                                                   dmelT3SpeciEnhCons_deltaSVM$end%in%dmelT3SpeciEnhCons_strict_deltaSVM$end),]

dmelT4SpeciEnhCons<-dmelT4SpeciEnhAnno[dmelT4SpeciEnhAnno$V4>0|dmelT4SpeciEnhAnno$V5>0|dmelT4SpeciEnhAnno$V6>0|dmelT4SpeciEnhAnno$V7>0|dmelT4SpeciEnhAnno$V8>0,]
dmelT4SpeciEnhCons<-dmelT4SpeciEnhCons[,c(1:4)]
colnames(dmelT4SpeciEnhCons)<-c("chr","start","end","anno" )
dmelT4SpeciEnhCons_deltaSVM<-merge(dmelT4SpeciEnhCons,dmelT4Enh_deltaSVM,by=c("chr","start","end"))
dmelT4SpeciEnhCons_relax_deltaSVM<-dmelT4SpeciEnhCons_deltaSVM[!(dmelT4SpeciEnhCons_deltaSVM$chr%in%dmelT4SpeciEnhCons_strict_deltaSVM$chr&
                                                                   dmelT4SpeciEnhCons_deltaSVM$start%in%dmelT4SpeciEnhCons_strict_deltaSVM$start&
                                                                   dmelT4SpeciEnhCons_deltaSVM$end%in%dmelT4SpeciEnhCons_strict_deltaSVM$end),]

dmelT5SpeciEnhCons<-dmelT5SpeciEnhAnno[dmelT5SpeciEnhAnno$V4>0|dmelT5SpeciEnhAnno$V5>0|dmelT5SpeciEnhAnno$V6>0|dmelT5SpeciEnhAnno$V7>0|dmelT5SpeciEnhAnno$V8>0,]
dmelT5SpeciEnhCons<-dmelT5SpeciEnhCons[,c(1:4)]
colnames(dmelT5SpeciEnhCons)<-c("chr","start","end","anno" )
dmelT5SpeciEnhCons_deltaSVM<-merge(dmelT5SpeciEnhCons,dmelT5Enh_deltaSVM,by=c("chr","start","end"))
dmelT5SpeciEnhCons_relax_deltaSVM<-dmelT5SpeciEnhCons_deltaSVM[!(dmelT5SpeciEnhCons_deltaSVM$chr%in%dmelT5SpeciEnhCons_strict_deltaSVM$chr&
                                                                   dmelT5SpeciEnhCons_deltaSVM$start%in%dmelT5SpeciEnhCons_strict_deltaSVM$start&
                                                                   dmelT5SpeciEnhCons_deltaSVM$end%in%dmelT5SpeciEnhCons_strict_deltaSVM$end),]


posEnh<-c(nrow(dmelT1SpeciEnhCons_relax_deltaSVM[dmelT1SpeciEnhCons_relax_deltaSVM$qValue<0.05]),
          nrow(dmelT2SpeciEnhCons_relax_deltaSVM[dmelT2SpeciEnhCons_relax_deltaSVM$qValue<0.05]),
          nrow(dmelT3SpeciEnhCons_relax_deltaSVM[dmelT3SpeciEnhCons_relax_deltaSVM$qValue<0.05]),
          nrow(dmelT4SpeciEnhCons_relax_deltaSVM[dmelT4SpeciEnhCons_relax_deltaSVM$qValue<0.05]),
          nrow(dmelT5SpeciEnhCons_relax_deltaSVM[dmelT5SpeciEnhCons_relax_deltaSVM$qValue<0.05]))
allEnh<-c(nrow(dmelT1SpeciEnhCons_relax_deltaSVM),nrow(dmelT2SpeciEnhCons_relax_deltaSVM),nrow(dmelT3SpeciEnhCons_relax_deltaSVM),nrow(dmelT4SpeciEnhCons_relax_deltaSVM),nrow(dmelT5SpeciEnhCons_relax_deltaSVM))

par(mfrow=c(1,1))
par(mar=c(5, 5, 4, 2) + 0.1)
bplot<-barplot(posEnh/allEnh, ylim=c(0, 0.5), main="Conserved enhancers",ylab = "Proportion of positive selection",col=c(pal[3],pal[7],pal[4],pal[5],pal[2]),
               cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
text(x=bplot,y=0-0.4*0.08,cex=1.5,adj = 1,labels = c("TP1", "TP2","TP3","TP4","TP5"),xpd = TRUE)
text(x=bplot,y=0+0.4*0.05,cex=1.2,labels = paste0(posEnh,"/",allEnh),xpd = TRUE)

##### for non-conserved enhancers (T1 specific enhancer in dmel without overlaped enhancers in dvir) #####
dmelT1SpeciEnh_nonCons_deltaSVM<-dmelT1Enh_deltaSVM[!(dmelT1Enh_deltaSVM$chr%in%dmelT1SpeciEnhCons_deltaSVM$chr&
                                                        dmelT1Enh_deltaSVM$start%in%dmelT1SpeciEnhCons_deltaSVM$start&
                                                        dmelT1Enh_deltaSVM$end%in%dmelT1SpeciEnhCons_deltaSVM$end),]
dmelT2SpeciEnh_nonCons_deltaSVM<-dmelT2Enh_deltaSVM[!(dmelT2Enh_deltaSVM$chr%in%dmelT2SpeciEnhCons_deltaSVM$chr&
                                                        dmelT2Enh_deltaSVM$start%in%dmelT2SpeciEnhCons_deltaSVM$start&
                                                        dmelT2Enh_deltaSVM$end%in%dmelT2SpeciEnhCons_deltaSVM$end),]
dmelT3SpeciEnh_nonCons_deltaSVM<-dmelT3Enh_deltaSVM[!(dmelT3Enh_deltaSVM$chr%in%dmelT3SpeciEnhCons_deltaSVM$chr&
                                                        dmelT3Enh_deltaSVM$start%in%dmelT3SpeciEnhCons_deltaSVM$start&
                                                        dmelT3Enh_deltaSVM$end%in%dmelT3SpeciEnhCons_deltaSVM$end),]
dmelT4SpeciEnh_nonCons_deltaSVM<-dmelT4Enh_deltaSVM[!(dmelT4Enh_deltaSVM$chr%in%dmelT4SpeciEnhCons_deltaSVM$chr&
                                                        dmelT4Enh_deltaSVM$start%in%dmelT4SpeciEnhCons_deltaSVM$start&
                                                        dmelT4Enh_deltaSVM$end%in%dmelT4SpeciEnhCons_deltaSVM$end),]
dmelT5SpeciEnh_nonCons_deltaSVM<-dmelT5Enh_deltaSVM[!(dmelT5Enh_deltaSVM$chr%in%dmelT5SpeciEnhCons_deltaSVM$chr&
                                                        dmelT5Enh_deltaSVM$start%in%dmelT5SpeciEnhCons_deltaSVM$start&
                                                        dmelT5Enh_deltaSVM$end%in%dmelT5SpeciEnhCons_deltaSVM$end),]
posEnh<-c(nrow(dmelT1SpeciEnh_nonCons_deltaSVM[dmelT1SpeciEnh_nonCons_deltaSVM$qValue<0.05]),
          nrow(dmelT2SpeciEnh_nonCons_deltaSVM[dmelT2SpeciEnh_nonCons_deltaSVM$qValue<0.05]),
          nrow(dmelT3SpeciEnh_nonCons_deltaSVM[dmelT3SpeciEnh_nonCons_deltaSVM$qValue<0.05]),
          nrow(dmelT4SpeciEnh_nonCons_deltaSVM[dmelT4SpeciEnh_nonCons_deltaSVM$qValue<0.05]),
          nrow(dmelT5SpeciEnh_nonCons_deltaSVM[dmelT5SpeciEnh_nonCons_deltaSVM$qValue<0.05]))
allEnh<-c(nrow(dmelT1SpeciEnh_nonCons_deltaSVM),nrow(dmelT2SpeciEnh_nonCons_deltaSVM),nrow(dmelT3SpeciEnh_nonCons_deltaSVM),nrow(dmelT4SpeciEnh_nonCons_deltaSVM),nrow(dmelT5SpeciEnh_nonCons_deltaSVM))

par(mfrow=c(1,1))
par(mar=c(5, 5, 4, 2) + 0.1)
bplot<-barplot(posEnh/allEnh, ylim=c(0, 0.5), main="Non-conserved enhancers",ylab = "Proportion of positive selection",col=c(pal[3],pal[7],pal[4],pal[5],pal[2]),
               cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
text(x=bplot,y=0-0.4*0.08,cex=1.5,adj = 1,labels = c("TP1", "TP2","TP3","TP4","TP5"),xpd = TRUE)
text(x=bplot,y=0+0.4*0.05,cex=1.2,labels = paste0(posEnh,"/",allEnh),xpd = TRUE)



