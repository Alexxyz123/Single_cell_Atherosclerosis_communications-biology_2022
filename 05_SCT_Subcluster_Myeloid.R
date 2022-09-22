library(Seurat)
library(tidyverse)
library(patchwork)
library(future)
options(future.globals.maxSize = 20 * 1024^3) #将全局变量上限调至20G
color_manual=c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
               "#F781BF","#999999","lightgrey","lightblue","#FC8D62","#F2F2F2")
##===提取细胞子集===##
setwd("/home/data/vip07/single-cell/data")
source("Resource/function.R")
rm(list=ls())
scRNA <- get(load("CellType/scRNA_celltype.Rdata"))
dir.create("Subset_Myeloid")
#####细胞亚群#######
Idents(scRNA) <- "seurat_clusters"
scRNAsub <- subset(scRNA, idents = c(5,7,8))
######净化数据+SCT标准化#######
DefaultAssay(scRNAsub) <- "RNA"
scRNAsub@assays$SCT <- NULL
colnames(scRNAsub@meta.data)
scRNAsub@meta.data <- scRNAsub@meta.data[,c("orig.ident", "nCount_RNA",
                                            "nFeature_RNA", "percent.mt", "percent.rb",
                                            "groups","celltype","batches")]
scRNAsub$main.celltype <- scRNAsub$celltype
scRNAsub$celltype <- NULL

##提取子集之后需要从头分析吗？需要重头分析
# var.all <- VariableFeatures(scRNA, assay = "RNA")
# var.sub <- FindVariableFeatures(scRNAsub, assay = "RNA") %>% 
#   VariableFeatures(assay = "RNA")
# var.share <- intersect(var.all, var.sub)
# length(var.share)

##数据标准化
#log标准化
# scRNAsub <- NormalizeData(scRNAsub) %>% FindVariableFeatures() %>%
#   ScaleData(features = rownames(scRNAsub))
###将scRNAsub重新分割为scRNAlist结构
Idents(scRNAsub) <- "orig.ident"
sample=unique(scRNAsub@meta.data$orig.ident)
scRNAlist=list()
for (i in 1:6) {scRNAlist[[i]] <- subset(scRNAsub, idents = sample[i])}
#SCT标准化
scRNAlist <- lapply(X=scRNAlist, FUN=function(x) SCTransform(x))
#保存SCT标准化之后的结果
save(scRNAlist, file = "Subset_Myeloid/scRNAlist_Myeloid_SCT.Rdata")

#######使用CCA整合样本#######
##FindAnchors
scRNA.features <- SelectIntegrationFeatures(scRNAlist, nfeatures = 3000)
scRNAlist <- PrepSCTIntegration(scRNAlist, anchor.features = scRNA.features, verbose = FALSE)
##样本较小，需要调整K.filter
k.num <- min(sapply(scRNAlist, ncol))
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, 
                                        normalization.method="SCT", 
                                        anchor.features=scRNA.features,
                                        k.filter =100 )
##Integrate
scRNAsub <- IntegrateData(anchorset = scRNA.anchors, normalization.method="SCT")
##PCA
scRNAsub <- RunPCA(scRNAsub, npcs = 100, verbose = FALSE)
ElbowPlot(scRNAsub, ndims = 100)
pc.num=1:40

##降维聚类
scRNAsub <-RunTSNE(scRNAsub,dims=pc.num) %>% 
  RunUMAP(dims=pc.num) %>%
  FindNeighbors(dims=pc.num) %>% 
  FindClusters(resolution =c(seq(0.7,0.1,-0.1)))
#查看每群的细胞数
table(scRNAsub$seurat_clusters)
#找到最佳的resolution，clusteree包
library(clustree)
clustree(scRNAsub@meta.data, prefix = "integrated_snn_res.")
#####分辨率为0.1
p1 <- DimPlot(scRNAsub, reduction = "tsne", group.by = "integrated_snn_res.0.1", 
              label = T,seed = 520) + ggtitle("integrated_snn_res.0.1")
p2 <- DimPlot(scRNAsub, reduction = "tsne", group.by = "integrated_snn_res.0.2", 
              label = T,seed = 520) + ggtitle("integrated_snn_res.0.2")

p3 <- DimPlot(scRNAsub, reduction = "tsne", group.by = "integrated_snn_res.0.3", 
              label = T,seed = 520) + ggtitle("integrated_snn_res.0.3")
p4 <- DimPlot(scRNAsub, reduction = "tsne", group.by = "integrated_snn_res.0.4", 
              label = T,seed = 520) + ggtitle("integrated_snn_res.0.4")
p=p1+p2+p3+p4
ggsave("Subset_Myeloid/Resolution_test.png", p, width = 12, height = 6)
ggsave("Subset_Myeloid/RES_0.3.png", p1, width = 12, height = 6)
#查看结果后采用0.1的分辨率
scRNAsub$seurat_clusters <- scRNAsub$integrated_snn_res.0.1
Idents(scRNAsub) <- "integrated_snn_res.0.1"

##查看CCA的整合效果
p1 <- DimPlot(scRNAsub, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(scRNAsub, reduction = "tsne", group.by = "orig.ident")
pc = p1 + p2 + plot_layout(guides = "collect")
ggsave("Subset_Myeloid/CCA_integr.png", pc, width = 8, height = 4)
#使用分面图查看效果
p <- DimPlot(scRNAsub, reduction = "tsne", group.by = "groups", 
             split.by = "orig.ident", ncol = 3) + NoLegend()
ggsave("Subset_Myeloid/CCA_facet.png", p, width = 9, height = 12)
####disease Vs. normal
p1 <- DimPlot(scRNAsub, reduction = "umap", group.by = "groups")
p2 <- DimPlot(scRNAsub, reduction = "tsne", group.by = "groups")
pc = p1 + p2 + plot_layout(guides = "collect")
ggsave("Subset_Myeloid/group_integr.png", pc, width = 8, height = 4)
##保存CCA之后的结果
save(scRNAsub, file = "Subset_Myeloid/scRNAsub_Myeloid_CCA.Rdata") 
scRNAsub <- get(load("Subset_Myeloid/scRNAsub_Myeloid_CCA.Rdata"))
#####鉴定数据亚群的细胞类型######
##提取各个Cluster的marker genes
ClusterMarker <- FindAllMarkers(scRNAsub, assay = "integrated", slot = "data", only.pos = T,
                                logfc.threshold = 0.25, min.pct = 0.25)
ClusterMarker <- ClusterMarker[,c(7,1:6)]
write.csv(ClusterMarker,'Subset_Myeloid/ClusterMarker_Myeloid.csv', row.names=F)
#ClusterMarker <- read.csv('Subset_Myeloid/ClusterMarker.csv')
#提取差异显著并且没有核糖体的Markers
ClusterMarker_noRibo <- ClusterMarker[!grepl("^RP[SL]", 
                                             ClusterMarker$gene, ignore.case = F),]
top = 15   #可根据需要调整
TopMarkers1 <- ClusterMarker_noRibo %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_logFC)
TopMarkers2 <- ClusterMarker_noRibo %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_logFC)
TopMarkers_noRibo <- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% 
  arrange(cluster)
write.csv(TopMarkers_noRibo,'Subset_Myeloid/TopMarkers_Myeloid.csv', row.names=F)

##导出一个方便对比的表格(每行一个cluster的marker gene)
TopMarkers2Lines(data = TopMarkers_noRibo, output = "Subset_Myeloid")
source("Resource/MarkerGeneList.R")
TopMarkers_noRibo <- read.csv('Subset_Myeloid/TopMarkers_Myeloid.csv')
MatchMarkerGene(data=TopMarkers2, output = "Subset_Myeloid")
##导入人工鉴定结果
cell.type <- c("TREM2-high",
               "Resident-like","Inflammatory")
Idents(scRNAsub) <- "seurat_clusters"
names(cell.type) <- levels(scRNAsub)
scRNAsub <- RenameIdents(scRNAsub, cell.type)
scRNAsub$celltype <- Idents(scRNAsub)
Idents(scRNAsub) <- "seurat_clusters"

##鉴定结果可视化
p1 <- DimPlot(scRNAsub, reduction = "umap", label = T,
              pt.size = 0.01,cols =color_manual) + NoLegend()
p2 <- DimPlot(scRNAsub, reduction = "umap", label = F, label.size = 2.5,
              group.by = "celltype",,pt.size = 0.01,cols =color_manual)
pc = p1|p2
ggsave("Subset_Myeloid/CellType_Custom.png", pc, width = 10, height = 4)

DimPlot(scRNAsub, reduction = 'tsne', group.by = 'celltype', 
        split.by = 'orig.ident', ncol = 3)
##保存结果####
save(scRNAsub, file = "Subset_Myeloid/scRNAsub_classify.Rdata")
load("~/single-cell/data/Subset_Myeloid/scRNAsub_classify.Rdata")
##===结果可视化===##
##鉴定结果可视化
pdf(file = "Subset_Myeloid/01_Myeloid_group.pdf",5,4)
DimPlot(scRNAsub, reduction = "umap",seed = 520,group.by = "groups",pt.size = 0.00001)+ 
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.1,"cm"))
dev.off()

pdf(file="Subset_Myeloid/02_Myeloid_celltype.pdf",6,4)
DimPlot(scRNAsub, reduction = "umap", label = T, label.size = 3,
        group.by = "celltype",pt.size = 0.00001)+ 
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.1,"cm")) 
dev.off()

pdf(file="Subset_Myeloid/03_Myeloid_cluster.pdf",5,4)
DimPlot(scRNAsub, reduction = "umap", label = T, label.size = 3,pt.size = 0.00001)+ 
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.1,"cm")) 
dev.off()
p1 <- DimPlot(scRNAsub, reduction = "umap", label = T) + NoLegend()

p2 <- DimPlot(scRNAsub, reduction = "umap", label = T, label.size = 2.5,
              group.by = "celltype") + NoLegend()
pc = p1|p2
ggsave("Subset_Myeloid/CellType_Custom.png", pc, width = 10, height = 4)

DimPlot(scRNAsub, reduction = 'umap', group.by = 'celltype', 
        split.by = 'orig.ident', ncol = 3)

##保存结果####
save(scRNAsub, file = "Subset_Myeloid/scRNAsub_classify.Rdata")

load("Subset_Myeloid/scRNAsub_classify.Rdata")
##===结果可视化===##
##Stackbar细胞丰度柱状图
tmp <-dplyr::select(scRNAsub@meta.data,c("groups", "celltype"))
df <- data.frame()
for(i in unique(tmp$groups)){
  df_i <- subset(tmp, tmp$groups==i) %>% pull(celltype) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("sample", "celltype", "value")
  df <- rbind(df, df_i)
}
source("Resource/function.R")
#按样本统计细胞类型
pdf("Subset_Myeloid/04_stackbar_myeloid.pdf",5,5)
ggplot(df, aes(x=sample, y=value, fill=celltype)) + 
  geom_bar(stat= "identity", position = "fill",width=0.4) + 
  scale_fill_manual(values = color_manual) +scale_y_continuous(expand=c(0,0))+
  labs(x = 'Group', y = 'Relative Abundance', title = 'Group composition') +
  theme_classic()+theme(legend.key.size=unit(0.4,'cm'))+
  guides(fill = guide_legend(title = 'Cell Type'))+#修改图例大小
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
dev.off()
ggsave('Subset_Myeloid/Stackbar_celltype.png', p, width = 8, height = 4.5)

##Marker基因可视化
Idents(scRNAsub) <- "celltype"
markers <-unique(c( "Trem2", "Cd9", "Lgals3", 
                    "Ctsb", "ABCA1",
                    "Cx3cr1","Folr2","MRC1","F13a1", "GAS6","Nlrp3", 
                   "Nfkbia", "IL1b", 
                   "Cxcl2", "S100A9","S100A12","S100A8","S100A4"))
markers <- CaseMatch(markers, rownames(scRNAsub))

# ## Featureplot标注marker基因
pdf(file ="Subset_Myeloid/Markers_features.pdf",21,9 )
FeaturePlot(scRNAsub, reduction = "umap", features = markers, ncol = 6,
                 cols = c("grey","#EF3B2C","firebrick3"))
dev.off()
###AS组有更多衰老的细胞
#Heatmap图
p <- DoHeatmap(scRNAsub, features = markers, size = 4,assay = "SCT",slot = "scale.data")
ggsave("Subset_Myeloid/Markers_heatmap.png", width = 12, height = 7)
Idents(scRNAsub) <- "seurat_clusters"
pdf(file = "Subset_Myeloid/05_complexheatmap.pdf",10,10)
DoHeatmap(scRNAsub, features=markers,slot ="scale.data",assay = "integrated",
          angle = 0,size=4)+
  NoLegend()
dev.off()
#Marker基因气泡图
mark_gene <- c("TREM2","ABCA1","CX3CR1","GAS6","NFKBIA",
               "IL1B")
pdf(file = "Subset_Myeloid/06_dotplot.pdf",12,6)
DotPlot(scRNAsub, features = mark_gene,cols = color_manual,group.by ="celltype",
        assay = "integrated" )
dev.off()


##Marker基因小提琴图
pdf(file = "Subset_Myeloid/07_Myeloid_vlnplot.pdf",8,8)
VlnPlot(scRNAsub, features = mark_gene, pt.size = 0.1,group.by = "celltype",ncol = 3,
        assay = "integrated",slot = "data")+theme(axis.title.x = element_blank())
dev.off()