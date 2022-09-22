##细胞鉴定和可视化####
##加载包
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
library(scCATCH)
library(RColorBrewer)
#加载数据
setwd("/home/data/ssy008/single-cell/data")
dir.create("CellType")
rm(list=ls())
source("Resource/function.R")
scRNA <- readRDS("SCTtransform/scRNA_SCTransform.rds")
#####SingleR预测####
scRNA$seurat_clusters=scRNA$integrated_snn_res.0.4
Idents(object = scRNA) <- "integrated_snn_res.0.4"
refdata <- get(load("Resource/SingleR_ref/ref_Human_all.RData"))
scRNA <- cell_identify(scRNA, refdata, output = "CellType")
p1 <- DimPlot(scRNA, reduction = "umap", label = T,group.by = "seurat_clusters") +
  NoLegend()
p2 <- DimPlot(scRNA, reduction = "umap", group.by = "SingleR", label = T,
              repel = T) + ggtitle("Predicted by SingleR") 
p <- p1 + p2
ggsave("CellType/CellType_SingleR.png", p, width = 10, height = 4)

####ssCATCH鉴定#####
# Step 1:寻找每一群得markers基因
clu_markers <- findmarkergenes(object = scRNA,
                               species = 'Human',
                               cluster = 'All')
##对应组织的细胞种类要全面
clu_markers <- findmarkergenes(object = scRNA,
                               species = 'Human',
                               cluster = 'All',
                               match_CellMatch = TRUE,
                               cancer = NULL,
                               tissue = c('Blood vessel',"Adventitia","Artery"),
                               cell_min_pct = 0.25,
                               logfc = 0.25,
                               pvalue = 0.05)
# Step 2:通过marker评分注释每一群得结果
clu_ann <- scCATCH(object = clu_markers$clu_markers,
                   species = 'Human',
                   cancer = NULL,
                   tissue = c('Blood vessel',"Adventitia","Artery"))
write.csv(clu_ann,file = "CellType/cluster_anno.csv")
###Marker基因鉴定####
##提取各个Cluster的marker genes
scRNA <- FindClusters(scRNA,resolution = 0.4)
ClusterMarker <- FindAllMarkers(scRNA, assay = "SCT", slot = "scale.data", 
                                only.pos = T,
                                logfc.threshold = 0.25, min.pct = 0.25)
ClusterMarker <- ClusterMarker[,c(7,1:6)]
write.csv(ClusterMarker,'CellType/ClusterMarker.csv', row.names=F)
##提取没有核糖体的Markers
ClusterMarker_noRibo <- ClusterMarker[!grepl("^RP[SL]", 
                                             ClusterMarker$gene, ignore.case = F),]
top = 15 #可根据需要调整
#P值为0.00
TopMarkers1 <- ClusterMarker_noRibo %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_logFC)
#P值<0.01
TopMarkers2 <- ClusterMarker_noRibo %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_logFC)
TopMarkers_noRibo <- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% 
  arrange(cluster)
write.csv(TopMarkers_noRibo,'CellType/TopMarkers_noRibo.csv', row.names=F)

##导出一个方便对比的表格(每行一个cluster的marker gene)
TopMarkers2Lines(data = TopMarkers2, output = "CellType")

##辅助查找marker基因
source("Resource/MarkerGeneList.R")
source("Resource/function.R")
TopMarkers_noRibo <- read.csv('CellType/TopMarkers_noRibo.csv')
MatchMarkerGene(data=TopMarkers2, output = "CellType")

####导入人工鉴定,manual gating结果#####
scRNA$seurat_clusters <- scRNA$integrated_snn_res.0.4
cell.type <- c("T Cells","T Cells","Endothelial Cells",
               "Smooth Muscle Cells","Smooth Muscle Cells",
               "Myeloid Cells","T Cells","Myeloid Cells","Myeloid Cells",
               "B Cells","T Cells","NK Cells","Endothelial Cells",
               "Dendritic Cells","Stem Cells","Smooth Muscle Cells",
               "Endothelial Cells","T Cells","Dendritic Cells")
Idents(scRNA) <- "seurat_clusters"
names(cell.type) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, cell.type)
scRNA$celltype <- Idents(scRNA)
Idents(scRNA) <- "seurat_clusters"
setwd("/home/data/vip07/single-cell/data/CellType")
####结果用t-SNE图展示#####
scRNA$groups <- scRNA$orig.ident
scRNA$groups <- recode(scRNA$groups, 
                      'AC_01' = "AC", 
                      'PA_01' = "PA",
                      'AC_02'= "AC",
                      'PA_02' = "PA",
                      'AC_03' = "AC",
                      'PA_03' = "PA")
p1 <- DimPlot(scRNA, reduction = "umap", label = T,pt.size = 0.00001)+
  labs(fill="Cluster")+ theme(legend.title = element_text(size = 8),
                               legend.text = element_text(size = 6),
                               legend.key.size = unit(0.2, "cm"),
                               legend.key.width = unit(0.1,"cm"))
pdf(file = "01_CellType_UMAP.pdf",5,4)
p1
dev.off()
ggsave("01_CellType_UMAP.png", p1, width = 5, height = 4)

pdf(file ="02_CellType_groups.pdf",5,4 )
DimPlot(scRNA, reduction = "umap", label = F, label.size = 2.5,
              group.by = "groups",pt.size = 0.00001)+ 
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.1,"cm"))
dev.off()
ggsave("02_CellType_groups.png", p2, width = 5, height = 4)

pdf(file ="03_CellType_type.pdf",5,4 )
DimPlot(scRNA, reduction = "umap", label = F, label.size = 2.5,
        group.by = "celltype",pt.size = 0.00001)+ 
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.1,"cm"))
dev.off()
ggsave("02_CellType_groups.png", p2, width = 5, height = 4)

##保存结果
save(scRNA, file = "scRNA_classify.Rdata")

#####Stackbar细胞丰度柱状图可视化####
tmp <- select(scRNA@meta.data, c("groups", "celltype"))
df <- data.frame()
for(i in unique(tmp$groups)){
  df_i <- subset(tmp, tmp$group==i) %>% pull(celltype) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("sample", "celltype", "value")
  df <- rbind(df, df_i)
}
source("/home/data/vip07/single-cell/data/Resource/function.R")
#按样本统计细胞类型
pdf(file = "04_Stackbar_celltype.pdf",4,4)
ggplot(df, aes(x=sample, y=value, fill=celltype)) + 
  geom_bar(stat= "identity", position = "fill",width=0.4) + 
  scale_fill_manual(values = color_manual) +scale_y_continuous(expand=c(0,0))+
  labs(x = 'Group', y = 'Relative Abundance') +
  theme_classic()+theme(legend.title = element_text(size = 10),
                        legend.text = element_text(size = 8),
                        legend.key.size = unit(0.4, "cm"),
                        legend.key.width = unit(0.4,"cm"))+
  guides(fill = guide_legend(title = 'Cell Type'))#修改图例大小
dev.off()
ggsave('Stackbar_celltype.png', p, width = 4, height = 4)
#按细胞类型统计样本
pdf()
p <- ggplot(df, aes(x=celltype, y=value, fill=sample)) +
  geom_bar(stat= "identity", position = "fill",width=0.4) + 
  scale_fill_manual(values = color_manual) +scale_y_continuous(expand=c(0,0))+
  labs(x = 'Group', y = 'Relative Abundance') +
  theme_classic()+theme(legend.key.size=unit(0.4,'cm'))+
  guides(fill = guide_legend(title = 'Cell Type'))#修改图例大小
ggsave('Stackbar_sample.png', p, width = 8, height = 4.5)

####Marker基因Featureplot可视化####
setwd("./CellType")
Idents(scRNA) <- "celltype"
markers <- c("TAGLN","ACTA2","MYH11","PDGFRB","CD34","PECAM1",
             "PTPRC","CD3E","CD8A","CD4","CD79A","CD79B",
             "CD68","CD14","LYZ","CX3CR1",
            "KLRB1","NKG7","FCGR3A","IL3RA","NRP1","KIT","CD44")
markers <- CaseMatch(markers, rownames(scRNA))
## Featureplot标注marker基因
p <- FeaturePlot(scRNA, reduction = "tsne", features = markers, ncol = 3,
                 cols = c("grey","#EF3B2C","firebrick3"))
ggsave("Markers_features.png", p, width = 7, height = 16)

###热图
class(Idents(scRNA))
head(Idents(scRNA))
Idents(scRNA) <- factor(Idents(scRNA),levels=c("Smooth Muscle Cells",
                                              "Endothelial Cells","T Cells",
                                              "B Cells","Myeloid Cells",
                                              "Dendritic Cells","NK Cells",
                                              "Stem Cells"))

pdf(file = "05_complexheatmap.pdf",18,10)
DoHeatmap(scRNA, features=markers,slot ="scale.data",assay = "SCT",angle = 45,
          size=4)
dev.off()

#Marker基因气泡图
pdf(file = "06_dotplot.pdf",18,10)
DotPlot(scRNA, features = markers,cols = col21)
dev.off()
##Marker基因小提琴图
marker_gene <- c("TAGLN","CD34","PECAM1",
                 "CD3E","CD79A",
                 "CD68","CD1C",
                 "KLRB1","KIT")
pdf(file = "04_Cellmarker_type.pdf",10,10)
VlnPlot(scRNA, features = marker_gene,pt.size = 0,
        sort="increasing")
dev.off()

ggsave("Markers_vlnplot.png", width = 15, height = 15)
####保存结果###
save(scRNA,file = "scRNA_celltype.Rdata")

