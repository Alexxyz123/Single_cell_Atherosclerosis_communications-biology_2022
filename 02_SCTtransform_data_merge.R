######多种方法，这儿我们使用SCT标准化后，进行锚点整合,再做降维聚类
#####对GSE159677进行数据质控###
library(Seurat)
library(tidyverse)
library(patchwork)
library(future)
##初始设置
setwd("/home/data/vip07/single-cell/data/SCTtransform")
rm(list = ls())
options(future.globals.maxSize = 20 * 1024^3) #将全局变量上限调至20G
#plan(multiprocess, workers = 20) 
#####整合数据#####
scRNAlist <- readRDS("scRNA_01_SCT.rds")
##SCTransform
scRNAlist <- lapply(X=scRNAlist, FUN=function(x) SCTransform(x))
##FindAnchors
scRNA.features <- SelectIntegrationFeatures(scRNAlist, nfeatures = 6000)
scRNAlist <- PrepSCTIntegration(scRNAlist, anchor.features = scRNA.features, verbose = FALSE)
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, normalization.method="SCT", 
                                        anchor.features=scRNA.features)
##Integrate
scRNA <- IntegrateData(anchorset = scRNA.anchors, normalization.method="SCT")
#####降维#####
##PCA
scRNA <- RunPCA(scRNA, npcs = 100, verbose = FALSE)
ggsave("ElbowPlot.png", plot = ElbowPlot(scRNA, ndims=100), width=8, height=6)
##set varable
pc.num=1:50
set.seed(846)
##Cluster
scRNA <- RunTSNE(scRNA,dims=pc.num,) %>% 
  RunUMAP(dims=pc.num) %>%
  FindNeighbors(dims=pc.num) %>% 
  FindClusters(resolution =c(seq(1.6,0.2,-0.2)))
#查看每群的细胞数
table(scRNA$seurat_clusters)
#找到最佳的resolution，clusteree包
library(clustree)
cluster_tree <- clustree(scRNA@meta.data, prefix = "integrated_snn_res.")
ggsave("cluster.png", cluster_tree, width = 8, height = 8)
#####添加变量
scRNA$groups <- scRNA$orig.ident
scRNA$groups <- recode(scRNA$orig.ident, 
                       'AC_01' = "AC", 
                       'PA_01' = "PA",
                       'AC_02'= "AC",
                       'PA_02' = "PA",
                       'AC_03' = "AC",
                       'PA_03' = "PA",)
scRNA$batches <- scRNA$orig.ident
scRNA$batches <- recode(scRNA$batches, 
                       'AC_01' = "batch_1", 
                       'PA_01' = "batch_1",
                       'AC_02'= "batch_2",
                       'PA_02' = "batch_2",
                       'AC_03' = "batch_3",
                       'PA_03' = "batch_3")
##RunTSNE
scRNA <- RunTSNE(scRNA, dims = pc.num, check_duplicates = FALSE) %>% RunUMAP(dims = pc.num, check_duplicates = FALSE)
plot1 = DimPlot(scRNA, reduction="tsne", group.by='orig.ident') + NoLegend() 
ggsave("tSNE.png", plot = plot1, width = 6, height = 5)
plot2 = DimPlot(scRNA, reduction="tsne", split.by='orig.ident', group.by='orig.ident', ncol=2) + NoLegend() 
ggsave("tSNE_sample.png", plot = plot2, width = 6, height = 15)
####UMAP图
plot1 = DimPlot(scRNA, reduction="umap",label = T) 
ggsave("UMAP.png", plot = plot1, width = 6, height = 5)
plot2 = DimPlot(scRNA, reduction="umap", split.by='orig.ident', group.by='orig.ident', ncol=2) + NoLegend() 
ggsave("UMAP_sample.png", plot = plot2, width = 6, height = 15)
######查看批次效应
p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "batches")
p2 <- DimPlot(scRNA, reduction = "umap", group.by = "batches")
pc = p1 + p2 + plot_layout(guides = "collect")
ggsave("batch_overview.png", pc, width = 8, height = 4)
#####查看正常组和疾病组
p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "groups")
p2 <- DimPlot(scRNA, reduction = "umap", group.by = "groups")
pc = p1 + p2 + plot_layout(guides = "collect")
ggsave("group_overview.png", pc, width = 8, height = 4)

##计算number of genes detected per UMI: 这个度量让我们对数据集的复杂性有了一个概念(每个UMI检测到的基因越多，我们的数据就越复杂)
scRNA$log10GenesPerUMI <- log10(scRNA$nFeature_RNA) / log10(scRNA$nCount_RNA)
##save seurat object
saveRDS(scRNA, "scRNA_SCTransform.rds")

