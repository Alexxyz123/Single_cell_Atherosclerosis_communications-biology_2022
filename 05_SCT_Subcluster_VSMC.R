library(Seurat)
library(tidyverse)
library(patchwork)
library(future)
rm(list = ls())
###设置手动颜色
color_manual=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
               "#F781BF","#999999","lightgrey","lightblue","#FC8D62","#F2F2F2")
options(future.globals.maxSize = 20 * 1024^3) #将全局变量上限调至20G
##===提取细胞子集===##
setwd("/home/data/vip07/single-cell/data")
scRNA <- get(load("CellType/scRNA_celltype.Rdata"))
dir.create("Subset_VSMC")
#####细胞亚群#######
Idents(scRNA) <- "seurat_clusters"
scRNAsub <- subset(scRNA, idents = c(3,4,15))
######净化数据+SCT标准化#######
DefaultAssay(scRNAsub) <- "RNA"
scRNAsub@assays$SCT <- NULL
colnames(scRNAsub@meta.data)
scRNAsub@meta.data <- scRNAsub@meta.data[,c("orig.ident", "nCount_RNA",
                                            "nFeature_RNA", "percent.mt", 
                                            "percent.rb","celltype","groups")]
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
save(scRNAlist, file = "Subset_VSMC/scRNAlist_VSMC_SCT.Rdata")
#获得scRNAlist数据
scRNAlist <- get(load("Subset_VSMC/scRNAlist_VSMC_SCT.Rdata"))

#######使用CCA整合样本#######
##FindAnchors
scRNA.features <- SelectIntegrationFeatures(scRNAlist, nfeatures = 5000)
scRNAlist <- PrepSCTIntegration(scRNAlist, anchor.features = scRNA.features, verbose = FALSE)
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, normalization.method="SCT", 
                                        anchor.features=scRNA.features)
##Integrate
scRNAsub <- IntegrateData(anchorset = scRNA.anchors, normalization.method="SCT")
##PCA
scRNAsub <- RunPCA(scRNAsub, npcs = 100, verbose = FALSE)
ElbowPlot(scRNAsub, ndims = 100)
pc.num=1:50

##降维聚类
scRNAsub <-RunTSNE(scRNAsub,dims=pc.num) %>% 
  RunUMAP(dims=pc.num) %>%
  FindNeighbors(dims=pc.num) %>% 
  FindClusters(resolution =c(seq(0.7,0.1,-0.1)))
#查看每群的细胞数
table(scRNAsub$seurat_clusters)
#找到最佳的resolution，clusteree包
library(clustree)
p1=clustree(scRNAsub@meta.data, prefix = "integrated_snn_res.")
ggsave("Subset_VSMC/cluster_tree.6000.png", p1, width = 12, height = 6)
save(scRNAsub,file ="Subset_VSMC/scRNAsub.fea=6000.Rdata" )
#save(scRNAsub,file ="Subset_VSMC/scRNAsub.fea=5000.Rdata" )
#save(scRNAsub,file ="Subset_VSMC/scRNAsub.fea=3000.Rdata" )
get(load("Subset_VSMC/scRNAsub.fea=5000.Rdata"))
####比较正常组和对照组的UMAP和t-SNE#####
scRNAsub$groups <- scRNAsub$orig.ident
scRNAsub$groups <- recode(scRNAsub$group, 
                      'AC_01' = "AC", 
                      'PA_01' = "PA",
                      'AC_02'= "AC",
                      'PA_02' = "PA",
                      'AC_03' = "AC",
                      'PA_03' = "PA")
p1 <- DimPlot(scRNAsub, reduction = "umap", group.by = "integrated_snn_res.0.1", 
              label = T,seed = 520) + ggtitle("integrated_snn_res.0.1")
p2 <- DimPlot(scRNAsub, reduction = "umap", group.by = "integrated_snn_res.0.2", 
              label = T,seed = 520) + ggtitle("integrated_snn_res.0.2")

p3 <- DimPlot(scRNAsub, reduction = "umap", group.by = "integrated_snn_res.0.3", 
              label = T,seed = 520) + ggtitle("integrated_snn_res.0.3")
p4 <- DimPlot(scRNAsub, reduction = "umap", group.by = "integrated_snn_res.0.4", 
              label = T,seed = 520) + ggtitle("integrated_snn_res.0.4")
p=p1+p2+p3+p4
ggsave("Subset_VSMC/Resolution_test.png", p, width = 12, height = 6)
ggsave("Subset_VSMC/RES_0.2.png", p2, width = 12, height = 6)
#查看结果后采用0.2的分辨率
scRNAsub$seurat_clusters <- scRNAsub$integrated_snn_res.0.2
Idents(scRNAsub) <- "integrated_snn_res.0.2"

##查看CCA的整合效果
p1 <- DimPlot(scRNAsub, reduction = "umap", group.by = "groups",seed = 520)
p2 <- DimPlot(scRNAsub, reduction = "tsne", group.by = "groups",seed = 520)
pc = p1 + p2 + plot_layout(guides = "collect")
ggsave("Subset_VSMC/CCA_integr.png", pc, width = 8, height = 4)
#使用分面图查看效果
p <- DimPlot(scRNAsub, reduction = "tsne", group.by = "groups", 
             split.by = "orig.ident", ncol = 3,seed = 520) + NoLegend()
ggsave("Subset_VSMC/CCA_facet.png", p, width = 9, height = 12)

##保存CCA之后的结果
save(scRNAsub, file = "Subset_VSMC/scRNAsub_VSMC_CCA.Rdata")     
#scRNAsub <- get(load("Subset_VSMC/scRNAsub_VSMC_CCA.Rdata"))
#####鉴定数据亚群的细胞类型######
##提取各个Cluster的marker genes
ClusterMarker <- FindAllMarkers(scRNAsub, assay = "integrated", slot = "data", only.pos = T,
                                logfc.threshold = 0.1, min.pct = 0.25)
ClusterMarker <- ClusterMarker[,c(7,1:6)]
write.csv(ClusterMarker,'Subset_VSMC/ClusterMarker_VSMC.csv', row.names=F)
#ClusterMarker <- read.csv('Subset_VSMC/ClusterMarker.csv')
#提取差异显著并且没有核糖体的Markers
ClusterMarker_noRibo <- ClusterMarker[!grepl("^RP[SL]", 
                                             ClusterMarker$gene, ignore.case = F),]
top = 20   #可根据需要调整
TopMarkers1 <- ClusterMarker_noRibo %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_logFC)
TopMarkers2 <- ClusterMarker_noRibo %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_logFC)
TopMarkers_noRibo <- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% 
  arrange(cluster)
write.csv(TopMarkers_noRibo,'Subset_VSMC/TopMarkers_VSMC.csv', row.names=F)

##导出一个方便对比的表格(每行一个cluster的marker gene)
TopMarkers2Lines(data = TopMarkers_noRibo, output = "Subset_VSMC")

##辅助查找marker基因
source("Resource/MarkerGeneList.R")
source("Resource/function.R")
#TopMarkers_noRibo <- read.csv('Subset_VSMC/TopMarkers_noRibo.csv')
MatchMarkerGene(data=TopMarkers_noRibo, output = "Subset_VSMC")

##使用main.celltype作为参考
p1 <- DimPlot(scRNAsub, reduction = "umap", label = T) + NoLegend()
p2 <- DimPlot(scRNAsub, reduction = "umap", label = T, label.size = 2.5,
              group.by = "groups") + NoLegend()
pc = p1|p2
ggsave("Subset_VSMC/CellType_Ref.png", pc, width = 10, height = 4)

##导入人工鉴定结果
cell.type <- c("Contractile", "Macrophage-like",
               "Fibroblast-like","Chondrocyte-like","Macrophage-like")
Idents(scRNAsub) <- "seurat_clusters"
names(cell.type) <- levels(scRNAsub)
scRNAsub <- RenameIdents(scRNAsub, cell.type)
scRNAsub$celltype <- Idents(scRNAsub)
Idents(scRNAsub) <- "seurat_clusters"

##鉴定结果可视化
load("~/single-cell/data/Subset_VSMC/scRNAsub_classify.Rdata")
pdf(file = "Subset_VSMC/01_VSMC_group.pdf",5,4)
DimPlot(scRNAsub, reduction = "umap",seed = 520,group.by = "groups",pt.size = 0.00001)+ 
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.1,"cm"))
dev.off()

pdf(file="Subset_VSMC/02_VSMC_celltype.pdf",6,4)
DimPlot(scRNAsub, reduction = "umap", label = T, label.size = 2.5,
        group.by = "celltype",pt.size = 0.00001)+ 
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.1,"cm")) 
dev.off()

pdf(file="Subset_VSMC/03_VSMC_cluster.pdf",5,4)
DimPlot(scRNAsub, reduction = "umap", label = T, label.size = 2.5,pt.size = 0.00001)+ 
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm"),
        legend.key.width = unit(0.1,"cm")) 
dev.off()
p1 <- DimPlot(scRNAsub, reduction = "umap", label = T) + NoLegend()

p2 <- DimPlot(scRNAsub, reduction = "umap", label = T, label.size = 2.5,
              group.by = "celltype") + NoLegend()
pc = p1|p2
ggsave("Subset_VSMC/CellType_Custom.png", pc, width = 10, height = 4)

DimPlot(scRNAsub, reduction = 'umap', group.by = 'celltype', 
        split.by = 'orig.ident', ncol = 3)

##保存结果####
save(scRNAsub, file = "Subset_VSMC/scRNAsub_classify.Rdata")

load("~/single-cell/data/Subset_VSMC/scRNAsub_classify.Rdata")
##===结果可视化===##
##Stackbar细胞丰度柱状图
tmp <-dplyr::select(scRNAsub@meta.data,c("groups", "celltype"))
df <- data.frame()
for(i in unique(tmp$groups)){
  df_i <- subset(tmp, tmp$groups==i) %>% pull(celltype) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("groups", "celltype", "value")
  df <- rbind(df, df_i)
}
source("Resource/function.R")
#按样本统计细胞类型
pdf("Subset_VSMC/Stackbar_celltype.pdf",5,10)
ggplot(df, aes(x=groups, y=value, fill=celltype)) + 
  geom_bar(stat= "identity", position = "fill",width=0.4) + 
  scale_fill_manual(values = color_manual) +scale_y_continuous(expand=c(0,0))+
  labs(x = 'Group', y = 'Relative Abundance', title = 'Group composition') +
  theme_classic()+theme(legend.key.size=unit(0.4,'cm'))+
  guides(fill = guide_legend(title = 'Cell Type'))#修改图例大小
dev.off()
##Marker基因可视化
Idents(scRNAsub) <- "celltype"
markers <-unique(c("MYH11","ATCA2","TAGLN","CD68","LGALS3","CXCR4","PTPRC",
                   "PDGFRB","MYH10", "ELN", 
                   "FN1","RUNX2","BMP2"))
markers <- CaseMatch(markers, rownames(scRNAsub))

# ## Featureplot标注marker基因
pdf("Subset_VSMC/featureplot.pdf",20,6)
FeaturePlot(scRNAsub, reduction = "umap", features = markers, ncol = 6,
                 cols = c("grey","#EF3B2C","firebrick3"))
dev.off()
###Featureplot展示细胞衰老和增殖情况
pdf("Subset_VSMC/age_prolideration.pdf",10,10)
FeaturePlot(scRNAsub, reduction = "umap", features = c("IGF1","CDKN1A","CDNK2B"), ncol = 3,
            cols = c("grey","#EF3B2C","firebrick3"),split.by = "groups")
dev.off()
###AS组有更多衰老的细胞
FeaturePlot(scRNAsub, reduction = "umap", features = c("Ki-67","PNCA"), ncol = 3,
            cols = c("grey","#EF3B2C","firebrick3"),split.by = "groups")
markers <- CaseMatch(c("MKI67","PCNA","TP53"), rownames(scRNAsub))
markers <- c("IGF1","CDKN1A")

Idents(scRNAsub)="groups"
scRNAsub$groups <- recode(scRNAsub$groups,
                          AS="AC",
                          normal="PA")

library(reshape2)
vln.df=as.data.frame(scRNAsub@assays$integrated@scale.data[markers,])
vln.df$gene <- markers
vln.df=reshape2::melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")

anno=scRNAsub@meta.data[,c("groups")] %>% as.data.frame() %>% 
  rename(Group=".")

anno$CB <- colnames(scRNAsub)
vln.df=inner_join(vln.df,anno,by="CB")

#为了控制画图的基因顺序
vln.df$gene=factor(vln.df$gene,levels = markers)
my_comparisons <- list(c("PA", "AC"))
vln.df%>%ggplot(aes(Group,exp))+geom_violin(aes(fill=gene),
                                            scale = "width",
                                            width=0.4)+
  facet_grid(vln.df$gene~.,scales = "free_y")+
  scale_fill_npg()+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_classic()+
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "ns")),
                     method.args = list(alternative = "greater"),
                     label = "p.signif"
  )+
  theme(
    strip.text.y = element_text(size = 8), # 设置分面的字字体大小、颜色、背景、边框，
    axis.text.y = element_text(size = 8),
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1,size = 8),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )
ggsave("output/05_vln_overall.pdf",width = 16,height =12,units = "cm")




#Heatmap图
p <- DoHeatmap(scRNAsub, features = markers, size = 4,assay = "SCT",slot = "scale.data")
ggsave("Subset_VSMC/Markers_heatmap.png", width = 12, height = 7)
Idents(scRNAsub) <- "seurat_clusters"
pdf(file = "Subset_VSMC/04_complexheatmap.pdf",10,10)
DoHeatmap(scRNAsub, features=markers,slot ="scale.data",assay = "SCT",angle = 0,size=4)+
  NoLegend()
dev.off()
#Marker基因气泡图
mark_gene <- c("TAGLN","ACTA2","CD68","CD74","FN1",
               "BMP2")
pdf(file = "Subset_VSMC/05_dotplot.pdf",12,6)
DotPlot(scRNAsub, features = mark_gene,cols = color_manual,group.by ="celltype" )
dev.off()

##Marker基因小提琴图
pdf(file = "Subset_VSMC/06_VSMC_vlnplot.pdf",6,10)
VlnPlot(scRNAsub, features = mark_gene, pt.size = 0,group.by = "celltype",ncol = 2,
        assay = "integrated",slot = "data")
dev.off()
p <- VlnPlot(scRNAsub, features = mark_gene, pt.size = 0,group.by = "celltype")
ggsave("Subset_VSMC/Markers_vlnplot.png", width = 10, height = 12)

