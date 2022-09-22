####首先是读入我们的测试集-GSE159677，包括3个动脉粥样斑块组织和其配对的动脉粥样斑块组织
######change the directory####
setwd("~/single-cell/data")
####加载包
library(Seurat)
library(tidyverse)
library(patchwork)
rm(list=ls())

##===创建seurat对象列表===##
##设置文件目录与样本名称
dir <- dir("GSE159677/")
dir <- paste0("GSE159677/", dir)
#查看文件顺序
dir                         
#按文件顺序给样本命名，名称不要以数字开头，中间不能有空格 
#AC:Calcified atherosclerotic core (AC) plaques
#PA:patient-matched proximal adjacent (PA) portions of carotid artery
samples_name = c('AC_01', 'PA_01', 'AC_02', 'PA_02', 'AC_03',
                 'PA_03')

##使用循环命令批量创建seurat对象
scRNAlist <- list()
HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") # 人类血液常见红细胞基因
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  #不设置min.cells过滤基因会导致CellCycleScoring报错：
  #Insufficient data values to produce 24 bins.  
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i],
                                       min.cells=3, min.features = 200)
  #给细胞barcode加个前缀，防止合并后barcode重名
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])   
  #计算线粒体基因比例，小鼠用第二行命令
  if(T){    
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-") 
    #scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^mt-")
  }
  #计算核糖体基因比例，小鼠用第二行命令
  if(T){
    scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^RB[SL]")
    #scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^Rb[sl]")
  }
  HB_m <- match(HB.genes_total,rownames(counts@assays$RNA))
  HB.genes <- rownames(counts@assays$RNA)[HB_m]
  HB.genes <- HB.genes[!is.na(HB.genes)]
  pbmc[["percent.HB"]]<-PercentageFeatureSet(counts,features=HB.genes)
}
#给列表命名并保存数据
names(scRNAlist) <- samples_name
#save(scRNAlist, file = "scRNAlist0.Rdata")
saveRDS(scRNAlist, file = "scRNAlist0.rds")

