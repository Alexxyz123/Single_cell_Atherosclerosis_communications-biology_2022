#设置目录
setwd("~/single-cell/data")
#加载数据
rm(list=ls())
load("scRNA_celltype.Rdata")
#加载包
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
#####1.提取T细胞子集做亚群鉴定####
##提取细胞子集####
dir.create("Subset_Tcell")
Idents(scRNA) <- "seurat_clusters"
table(scRNA@meta.data$SCT_snn_res.0.3)
cluster <- scRNA@meta.data %>% filter(celltype=="T_cells")
table(cluster$SCT_snn_res.0.3)
scRNAsub <- subset(scRNA, idents = c(0, 2, 13))
###方案二：也可以选使用celltype提取
if(F){
  names(table(scRNA$celltype))
  Idents(scRNA) <- "celltype"
  scRNAsub <- subset(scRNA, idents = c("T_naive", "CD4 T_em", "T_reg", 
                                       "CD8_Trm", "CD8_T_naive", "T cells"))
}
##净化数据
DefaultAssay(scRNAsub) <- "RNA"
scRNAsub@assays$SCT <- NULL
colnames(scRNAsub@meta.data)
scRNAsub@meta.data <- scRNAsub@meta.data[,c("orig.ident", "nCount_RNA","batches",
                                            "nFeature_RNA", "percent.mt",  "percent.rb","proj",
                                            "S.Score", "G2M.Score", "Phase","celltype")]
scRNAsub$main.celltype <- scRNAsub$celltype
scRNAsub$celltype <- NULL

##提取子集之后需要从头分析吗？需要从头分析
var.all <- VariableFeatures(scRNA, assay = "RNA")
var.sub <- FindVariableFeatures(scRNAsub, assay = "RNA") %>% 
  VariableFeatures(assay = "RNA")
#大的scRNA数据集和子集的高变异基因差别大，所以需要重新分析
var.share <- intersect(var.all, var.sub)
length(var.share)

####数据标准化####
#log标准化
scRNAsub <- NormalizeData(scRNAsub) %>% FindVariableFeatures() %>%
  ScaleData(features = rownames(scRNAsub))
#SCT标准化
scRNAsub <- SCTransform(scRNAsub, return.only.var.genes = F)
save(scRNAsub, file = "Subset_Tcell/scRNAsub_SCT.Rdata")

##使用harmony整合数据
library(harmony)
scRNAsub <- RunPCA(scRNAsub, npcs=50, verbose=FALSE)
#sample区别
DimPlot(scRNAsub, reduction = "pca", group.by = "orig.ident")
#正常和对照组
DimPlot(scRNAsub, reduction = "pca", group.by = "batches")
##校正batch或者sample
scRNAsub <- RunHarmony(scRNAsub, group.by.vars="batches", assay.use="SCT",
                       max.iter.harmony = 15)
ElbowPlot(scRNAsub, ndims = 50)
pc.num=1:35

##降维聚类
scRNAsub <- RunTSNE(scRNAsub, reduction="harmony", dims=pc.num) %>% 
  RunUMAP(reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num) %>% 
  FindClusters(resolution=0.3) %>% FindClusters(resolution=0.5) %>%
  FindClusters(resolution=0.8) %>% FindClusters(resolution=1.2)
p1 <- DimPlot(scRNAsub, reduction = "umap", group.by = "SCT_snn_res.0.3", 
              label = T) + ggtitle("SCT_snn_res.0.3")
p2 <- DimPlot(scRNAsub, reduction = "umap", group.by = "SCT_snn_res.0.5", 
              label = T) + ggtitle("SCT_snn_res.0.5")
p3 <- DimPlot(scRNAsub, reduction = "umap", group.by = "SCT_snn_res.0.8", 
              label = T) + ggtitle("SCT_snn_res.0.8")
p4 <- DimPlot(scRNAsub, reduction = "umap", group.by = "SCT_snn_res.1.2", 
              label = T) + ggtitle("SCT_snn_res.1.2")
p <- (p1|p2)/(p3|p4)
ggsave("Subset_Tcell/Resolution_test.png", p, width = 10, height = 8)
#查看结果后采用0.3的分辨率
scRNAsub$seurat_clusters <- scRNAsub$SCT_snn_res.0.3
Idents(scRNAsub) <- "SCT_snn_res.0.3"

##查看harmony的整合效果
p1 <- DimPlot(scRNAsub, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(scRNAsub, reduction = "umap", group.by = "batches")
pc = p1 + p2 + plot_layout(guides = "collect")
ggsave("Subset_Tcell/Harmony_integr.png", pc, width = 8, height = 4)
#使用分面图查看效果
p <- DimPlot(scRNAsub, reduction = "umap", group.by = "orig.ident", 
             split.by = "orig.ident", ncol = 3) + NoLegend()
ggsave("Subset_Tcell/Harmony_facet.png", p, width = 9, height = 12)

##保存结果
save(scRNAsub, file = "Subset_Tcell/scRNAsub_harmony.Rdata")     


##===鉴定数据子集的细胞类型===##
##提取各个Cluster的marker genes
ClusterMarker <- FindAllMarkers(scRNAsub, assay = "SCT", slot = "data", only.pos = T,
                                logfc.threshold = 0.25, min.pct = 0.1)
ClusterMarker <- ClusterMarker[,c(7,1:6)]
write.csv(ClusterMarker,'Subset_Tcell/ClusterMarker.csv', row.names=F)
#ClusterMarker <- read.csv('Subset/ClusterMarker.csv')
#提取差异显著的marker genes
top = 15   #可根据需要调整
TopMarkers1 <- ClusterMarker %>% filter(p_val_adj == 0) %>% group_by(cluster) %>% 
  top_n(n = top, wt = avg_logFC)
TopMarkers2 <- ClusterMarker %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>%
  top_n(n = top, wt = avg_logFC)
TopMarkers <- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(TopMarkers,'Subset_Tcell/TopMarkers.csv', row.names=F)

##提取没有核糖体的Markers
ClusterMarker_noRibo <- ClusterMarker[!grepl("^RP[SL]", 
                                             ClusterMarker$gene, ignore.case = F),]
top = 15   #可根据需要调整
TopMarkers1 <- ClusterMarker_noRibo %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_logFC)
TopMarkers2 <- ClusterMarker_noRibo %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_logFC)
TopMarkers_noRibo <- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% 
  arrange(cluster)
write.csv(TopMarkers_noRibo,'Subset_Tcell/TopMarkers_noRibo.csv', row.names=F)

##辅助查找marker基因
source("Resource/MarkerGeneList.R")
source("Resource/function.R")
##导出一个方便对比的表格(每行一个cluster的marker gene)
TopMarkers2Lines(data = TopMarkers_noRibo, output = "Subset_Tcell")

#TopMarkers_noRibo <- read.csv('Subset/TopMarkers_noRibo.csv')
MatchMarkerGene(data=TopMarkers_noRibo, output = "Subset_Tcell")

##使用main.celltype作为参考
p1 <- DimPlot(scRNAsub, reduction = "umap", label = T) + NoLegend()
p2 <- DimPlot(scRNAsub, reduction = "umap", label = T, label.size = 2.5,
              group.by = "main.celltype") + NoLegend()
pc = p1|p2
ggsave("Subset_Tcell/CellType_Ref.png", pc, width = 10, height = 4)

##导入人工鉴定结果
#######cell.type <- c(+++++++)
Idents(scRNAsub) <- "seurat_clusters"
names(cell.type) <- levels(scRNAsub)
scRNAsub <- RenameIdents(scRNAsub, cell.type)
scRNAsub$celltype <- Idents(scRNAsub)
Idents(scRNAsub) <- "seurat_clusters"

##鉴定结果可视化
p1 <- DimPlot(scRNAsub, reduction = "umap", label = T) + NoLegend()
p2 <- DimPlot(scRNAsub, reduction = "umap", label = T, label.size = 2.5,
              group.by = "celltype") + NoLegend()
pc = p1|p2
ggsave("Subset/CellType_Custom.png", pc, width = 10, height = 4)

DimPlot(scRNAsub, reduction = 'umap', group.by = 'celltype', 
        split.by = 'orig.ident', ncol = 4)

##保存结果
save(scRNAsub, file = "Subset/scRNAsub_classify.Rdata")


##===结果可视化===##
##Stackbar细胞丰度柱状图
tmp <- select(scRNAsub@meta.data, c("orig.ident", "celltype"))
df <- data.frame()
for(i in unique(tmp$orig.ident)){
  df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(celltype) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("sample", "celltype", "value")
  df <- rbind(df, df_i)
}
source("Resource/function.R")
#按样本统计细胞类型
p <- ggplot(df, aes(x=sample, y=value, fill=celltype)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = col21) +
  labs(x = 'Samples', y = 'Relative Abundance', title = 'Samples composition') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45))
ggsave('Subset/Stackbar_celltype.png', p, width = 8, height = 4.5)
#按细胞类型统计样本
p <- ggplot(df, aes(x=celltype, y=value, fill=sample)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = col21) +
  labs(x = 'Cell type', y = 'Relative Abundance', title = 'Source of cell type') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45))
ggsave('Subset/Stackbar_sample.png', p, width = 8, height = 4.5)

##Marker基因可视化
Idents(scRNAsub) <- "celltype"
markers <- c("CCR7", "FOXP3", "CD8A", "GZMK", "NKG7", "GNLY", "TIMP1", "MT-ND3",
             "KLRB1", "HSPA1A", "ICA1", "PLCG2", "ISG15", "STMN1")
markers <- CaseMatch(markers, rownames(scRNAsub))

## Featureplot标注marker基因
p <- FeaturePlot(scRNAsub, reduction = "umap", features = markers, ncol = 4)
ggsave("Subset/Markers_features.png", p, width = 12, height = 10)

#Heatmap图
p <- DoHeatmap(scRNAsub, features = markers, size = 4)
ggsave("Subset/Markers_heatmap.png", width = 12, height = 7)

#Marker基因气泡图
p <- DotPlot(scRNAsub, features = markers) + 
  theme(axis.text.x = element_text(angle=45, vjust = 0.6))
ggsave("Subset/Markers_dotplot.png", width = 12, height = 6)

##Marker基因小提琴图
p <- VlnPlot(scRNAsub, features = markers, pt.size = 0, stack = T)
ggsave("Subset/Markers_vlnplot.png", width = 10, height = 6)
