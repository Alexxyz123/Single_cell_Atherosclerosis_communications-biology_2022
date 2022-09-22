library(dplyr)
library(cluster)
library(pheatmap)
library(tibble)
library(factoextra)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(ggpubr)
library(limma)
library(GEOquery)
library(tableone)
library(sva)
#####1.get the data#####
setwd("/home/data/vip07/single-cell/data/bulk_RNA")
rm(list = ls())
GSE90074 <-read.csv("GSE90074.csv")
names(GSE90074)[1] <- "ID"
load("gset")
gset=gset[[1]]
pdata <- pData(gset)
#View(pdata)
######2.combate and visulazition####
exprSet=exprs(gset)
####class=0
batch = pdata$`batchid:ch2`
exprSet <- ComBat(dat = exprSet,batch = batch)
boxplot(exprSet,outline=FALSE, notch=T, las=2)
#######3.read the annotation file#####
GPL <- data.table::fread(file = "GPL6480.txt",skip = 17) %>% select(c("ID","GENE_SYMBOL"))
exprSet = as.data.frame(exprSet) %>% rownames_to_column(var = "ID") 
exp <- merge(GPL,exprSet,by="ID") %>% select(-ID) %>%
  .[.$GENE_SYMBOL!="",]
####max value
exprSet_symbol <- aggregate(x =exp[,2:ncol(exp)],
                            by = list(exp$GENE_SYMBOL),
                            FUN =max)
names(exprSet_symbol)[1]<-"ID"
exprSet_symbol <- exprSet_symbol %>% column_to_rownames(var = "ID")
exprSet_symbolt=normalizeBetweenArrays(exprSet_symbol)
boxplot(exprSet_symbol,outline=FALSE, notch=T, las=2)
#画主成分分析图需要加载这两个包
library(FactoMineR)
library(factoextra) 
# PCA
dat.pca <- PCA(t(exprSet_symbol), graph = FALSE)
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = pdata$`batchid:ch2`, # color by groups
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups")
######4.clustering####
TF_name <- c("JUN","ETS1","EPAS1","JUND","STAT3","RARG"  ,"KLF13","NR3C1")
library(ConsensusClusterPlus)
dataset<-exprSet_symbol[TF_name,] %>% as.matrix()
#####回归分析###
dataset <- t(dataset) %>% as.data.frame() %>% 
  rownames_to_column(var = "ID")
pdata=GSE90074
data <- merge(dataset,pdata,by="ID")
####logistic回归分析
b=rep(0,8)
data$event <- data$Obstructive_cad
for (i in 1:length(TF_name)) {
  tmp <- TF_name[i]
  a <- glm(event~data[,tmp],data = data)
  a <- summary(a)
  b[i] <- if_else(as.numeric(a$coefficients[,"Pr(>|t|)"][2])<0.05,1,0)
}
b <- as.data.frame(b)
#####begin clustering 
dataset <- dataset %>% column_to_rownames(var = "ID") %>% as.matrix() %>% t()
cluster_method=c("km","pam")
distance=c("pearson", "spearman", "euclidean","maximum","canberra",
           "manhattan","minkowski")
#####data prepare
object=list()
PAC_result=list()
Kvec=2:10
x1 = 0.1
x2 = 0.9 
for (i in 1:2) {
  for (j in 1:7) {
    result=ConsensusClusterPlus(dataset,maxK = 10,reps =50,
                                pItem = 0.8,pFeature = 0.8,
                                clusterAlg = cluster_method[i],
                                title = paste(cluster_method[i],distance[j],sep = "-"),
                                distance = distance[j],writeTable = F,plot="png")
    object[[paste("result",cluster_method[i],distance[j],sep = "_")]]=result
    PAC = rep(NA,length(Kvec)) 
    names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
    ####obtain the PAC with different number of clusters
    for(h in Kvec){
      M = result[[h]]$consensusMatrix
      Fn = ecdf(M[lower.tri(M)])
      PAC[h-1] = Fn(x2) - Fn(x1)
    }
    PAC_result[[paste("PAC",cluster_method[i],distance[j],sep = "_")]]=PAC
  }
  
}
####find the best clusters####
result=matrix(nrow = 14,ncol = 9) 
result=as.data.frame(result)
rownames(result)=names(PAC_result)
colnames(result)=names(PAC_result[[1]])
for (i in 1:14) {
  result[i,]=PAC_result[[i]]
}
####obtained the lowest PAC value and get the place
result_place=as.data.frame(matrix(nrow = 14,ncol = 4))
for (i in 1:14) {
  result_place[i,1]=which.min(result[i,])
  result_place[i,2]=colnames(result[which.min(result[i,])])
  result_place[i,3]=rownames(result[i,])
  result_place[i,4]=result[i,which.min(result[i,])]
}
####the lowest PAC
min(as.numeric(result_place[,4]))

#######km_spearman
rcc <- ConsensusClusterPlus(dataset,maxK = 10,reps =50,
                     pItem = 0.8,pFeature = 1,
                     clusterAlg ="km",
                     title ="example",
                     distance = "euclidean",writeTable = F,plot="png")
###clutser =2
P_value <-matrix(nrow = 1,ncol = 14)
n <- 1
for (i in 1:2) {
  for (j in 1:7) {
    method <- paste("result",cluster_method[i],distance[j],sep = "_")
    cluster=object[[method]][[2]]$consensusClass
    ID <- colnames(exprSet_symbol)
    group_info <- data.frame(ID=ID,cluster=cluster)
    group_info <- merge(group_info,GSE90074,by="ID") %>% as.data.frame()
    # write.csv(group_info,file = "group_info.csv")
    # group_info <- read.csv("group_info.csv")
    group_info$cluster <- as.factor(group_info$cluster)
    group_info$event <- as.factor(group_info$Obstructive_cad)
    #####
    tableone_groups <- CreateTableOne(
      vars = c("event"),#指定纳入的变量
      strata = 'cluster', #指定分组变量#若不指定则对总体分析做表#
      data = group_info, #指定数据集
    ) #指定分类变量
    tmp <- print(tableone_groups$CatTable)
    P_value[1,n] <- tmp[2,3]
    n <- n+1
  }
}
