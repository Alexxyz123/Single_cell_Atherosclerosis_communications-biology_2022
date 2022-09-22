######这儿我们首先分别对每个样本进行质控，在进行SCT标准化+锚点整合#####
#####对GSE159677进行数据质控###
library(Seurat)
library(tidyverse)
library(patchwork)
library(future)
##数据质控
setwd("/home/data/ssy008/single-cell/data")
rm(list = ls())
options(future.globals.maxSize = 20 * 1024^3) #将全局变量上限调至20G
scRNAlist <- readRDS("scRNAlist0.rds")
setwd("./SCTtransform")
#设置可能用到的主题
theme.set1 = theme(axis.title.x=element_blank(), 
                   axis.text.x=element_blank(), 
                   axis.ticks.x=element_blank())
theme.set2 = theme(axis.title.x=element_blank())
#设置绘图元素
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb")
cell_number_before_QC=vector()
#加入分组栏
group = "orig.ident"
for (i in 1:length(scRNAlist)) {
  print("每个样本细胞数")
  print(as.character(unique(scRNAlist[[i]]@meta.data$orig.ident)))
  print(nrow(scRNAlist[[i]]@meta.data))
  cell_number_before_QC[i]=nrow(scRNAlist[[i]]@meta.data)
  ####number of genes detected per UMI: 
  #这个度量让我们对数据集的复杂性有了一个概念(每个UMI检测到的基因越多，
  #我们的数据就越复杂)
  scRNAlist[[i]]$log10GenesPerUMI <- log10(scRNAlist[[i]]$nFeature_RNA) / log10(scRNAlist[[i]]$nCount_RNA)
  ##绘制质控小提琴图
  plots = list()
  for(j in seq_along(plot.featrures)){
    plots[[j]] = VlnPlot(scRNAlist[[i]], group.by=group,
                         features = plot.featrures[j]) + theme.set2 + NoLegend()}
  violin <- wrap_plots(plots = plots, nrow=2) 
  name=paste0("before_QC_",as.character(unique(scRNAlist[[i]]@meta.data$orig.ident)))
  ggsave(paste0(name,".png"), plot = violin, width = 10, height = 8)  
  
}
######细胞筛选——设置质控标准####
minGene=c()
maxGene=c()
pctMT=c()
pctRB=c()
for (i in 1:length(scRNAlist)) {minGene[i] <- quantile(scRNAlist[[i]]$nFeature_RNA,0.1)
maxGene[i] <- quantile(scRNAlist[[i]]$nFeature_RNA,0.9)
pctMT[i]=quantile(scRNAlist[[i]]$percent.mt,0.95)
pctRB[i]=quantile(scRNAlist[[i]]$percent.rb,0.9)
}
minGenesPerUMI=rep(0.8,6)
##数据质控并绘制小提琴图
cell_number_after_QC=vector()
for (i in 1:length(scRNAlist)) {
  ####number of genes detected per UMI: 
  #这个度量让我们对数据集的复杂性有了一个概念(每个UMI检测到的基因越多，
  #我们的数据就越复杂)
  scRNAlist[[i]]$log10GenesPerUMI <- log10(scRNAlist[[i]]$nFeature_RNA) / log10(scRNAlist[[i]]$nCount_RNA)
  ###开始筛选
  scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = nFeature_RNA > minGene & nFeature_RNA < 
                    maxGene & percent.mt < pctMT & percent.rb < pctRB & log10GenesPerUMI>
                    minGenesPerUMI)
  #####筛选后的细胞数
  print("质控后每个样本细胞数")
  print(as.character(unique(scRNAlist[[i]]@meta.data$orig.ident)))
  print(nrow(scRNAlist[[i]]@meta.data))
  cell_number_after_QC[i]=nrow(scRNAlist[[i]]@meta.data)
  print("质控后的细胞数百分比")
  print(cell_number_after_QC[i]/cell_number_before_QC[i])
  ####质控后作图######
  plots = list()
  for(j in seq_along(plot.featrures)){
    plots[[j]] = VlnPlot(scRNAlist[[i]], group.by=group,
                         features = plot.featrures[j]) + theme.set2 + NoLegend()}
  violin <- wrap_plots(plots = plots, nrow=2) 
  name=paste0("after_QC_",as.character(unique(scRNAlist[[i]]@meta.data$orig.ident)))
  ggsave(paste0(name,".png"), plot = violin, width = 10, height = 8)  
}

##保存分别质控后的结果#####
saveRDS(scRNAlist, "scRNA_01_SCT.rds")
