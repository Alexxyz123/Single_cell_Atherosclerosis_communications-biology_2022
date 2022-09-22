###vignette(https://github.com/saeyslab/nichenetr/blob/master/vignettes/circos.md)
####加载包####
library(nichenetr)
library(Seurat)
library(tidyverse)
library(circlize)
rm(list=ls())
setwd("/home/data/vip07/single-cell/data")
####读取总的单细胞数据#####
scRNA = get(load("CellType/scRNA_celltype.Rdata"))
#####加载既往已有的ligand-target模型####
ligand_target_matrix = readRDS("Nichenetr/ligand_target.rds")
#存在ligand-receptor数据库
lr_network = readRDS("Nichenetr/lr_network.rds")
# weighted_networks列表包含两个数据框，
#lr_sig是配体-受体权重信号网络，gr是配体-靶基因权重调控网络
weighted_networks = readRDS("Nichenetr/weighted_networks.rds")
#直接进行分析
####更改进一步的细胞分类
scRNA_VSMC <- get(load("Subset_VSMC/dorothea/scRNAsub_doro.Rdata"))
scRNA_Myeloid <- get(load("Subset_Myeloid/scRNAsub_classify.Rdata"))
scRNA_Myeloid$groups <- scRNA_Myeloid$orig.ident
scRNA_Myeloid$groups <- recode(scRNA_Myeloid$groups, 
                          'AC_01' = "AC", 
                          'PA_01' = "PA",
                          'AC_02'= "AC",
                          'PA_02' = "PA",
                          'AC_03' = "AC",
                          'PA_03' = "PA")
#####overall 分析#########
# cell_total <- scRNA@meta.data %>% rownames_to_column(var = "ID")
# cell_VSMC <- scRNA_VSMC@meta.data %>% subset(.,select=c("celltype","groups"))%>% rownames_to_column(var = "ID")
# cell_Myeloid <- scRNA_Myeloid@meta.data %>% subset(.,select=c("celltype","groups")) %>% rownames_to_column(var = "ID")
# cells_name <- setdiff(cell_total$ID,c(cell_Myeloid$ID,cell_VSMC$ID))
# cell_other <- cell_total[cell_total$ID%in%cells_name,] %>% subset(.,select=c("celltype","groups","ID"))
# total <- rbind(cell_other,cell_VSMC,cell_Myeloid)
# total$subtype <- total$subtype %>% as.character()
# names(total)[1] <-"subtype" 
# rownames(total) <- NULL
# total <- as.data.frame(total) %>% column_to_rownames(.,var = "ID")
# scRNA <-AddMetaData(scRNA,total,col.name = c('subtype',"groups"))
# rm(cells_name,cell_other)
# Idents(scRNA) <- "subtype"
# nichenet_output_total = nichenet_seuratobj_aggregate(seurat_obj = scRNA, 
#                                                      top_n_ligands = 20,
#                                                      receiver = "Macrophage-like",
#                                                      sender =unique(scRNA$subtype),
#                                                      condition_colname = "groups", 
#                                                      condition_oi = "AC", 
#                                                      condition_reference = "PA", 
#                                                      ligand_target_matrix = ligand_target_matrix, 
#                                                      lr_network = lr_network, 
#                                                      weighted_networks = weighted_networks, 
#                                                      organism = "human",assay_oi = "RNA")
DefaultAssay(scRNA) <- "RNA"
####分裂成PA和AC两个小类，分别进行分析
Idents(scRNA) <- "celltype"
nichenet_output_total = nichenet_seuratobj_aggregate(seurat_obj = scRNA, 
                                                     top_n_ligands = 20,
                                                     receiver = "Smooth Muscle Cells",
                                                     sender =unique(scRNA$celltype),
                                                     condition_colname = "groups", 
                                                     condition_oi = "AC", 
                                                     condition_reference = "PA", 
                                                     ligand_target_matrix = ligand_target_matrix, 
                                                     lr_network = lr_network, 
                                                     weighted_networks = weighted_networks, 
                                                     organism = "human",assay_oi = "RNA")

nichenet_output_total$top_ligands
# 查看top20 ligands在各个细胞亚群中表达情况
pdf("Nichenetr/01_overall_ligand.pdf",10,8)
DotPlot(scRNA,cols = "RdYlBu",features = nichenet_output_total$top_ligands,split.by ="groups" )+
  theme(legend.key.size=unit(0.4,'cm'),legend.key.height = unit(0.4,'cm'),
        legend.text = element_text(size = 8),legend.title =element_text(size = 12) )+
  RotatedAxis()
dev.off()
####SMC的受体
pdf("Nichenetr/02_overall_receptor.pdf",8,5)
nichenet_output_total$ligand_receptor_heatmap+  
  scale_fill_gradient2(low = "skyblue",  high = "firebrick3", breaks = c(0,1,1.5)) + 
  xlab("Receptors Expressed in VSMCs")+
  ylab("Prioritized Ligands Modulating VSMCs Response")
dev.off()
####Target基因的dotplot
pdf("Nichenetr/03_overall_target.pdf",30,5)
DotPlot(scRNA %>% subset(idents = c("Smooth Muscle Cells")), 
        features = nichenet_output_total$top_targets, 
        split.by = "groups",cols = "RdYlBu") + RotatedAxis()+
  theme(axis.text.x = element_text(size=6) ,legend.key.size=unit(0.4,'cm'),legend.key.height = unit(0.4,'cm'),
        legend.text = element_text(size = 8),legend.title =element_text(size = 12) )+
  RotatedAxis()
dev.off()

#########VSMC和endothelial以及monocyte分群分析########

####获得cells的ID号
cell_VSMC <- scRNA_VSMC@meta.data %>% subset(.,select=c("celltype","groups"))
cell_Myeloid <- scRNA_Myeloid@meta.data %>% subset(.,select=c("celltype","groups")) 
#####获取细胞子集，并进行亚群分类
ID <- c(rownames(cell_VSMC),rownames(cell_Myeloid))
scRNAsub <- scRNA[,ID]
######净化数据+SCT标准化#######
DefaultAssay(scRNAsub) <- "RNA"
scRNAsub@assays$SCT <- NULL
colnames(scRNAsub@meta.data)
scRNAsub@meta.data <- scRNAsub@meta.data[,c("orig.ident", "nCount_RNA",
                                            "nFeature_RNA", "percent.mt", 
                                            "percent.rb","celltype","groups")]
scRNAsub$main.celltype <- scRNAsub$celltype
scRNAsub$celltype <- NULL
###将scRNAsub重新分割为scRNAlist结构
Idents(scRNAsub) <- "orig.ident"
sample=unique(scRNAsub@meta.data$orig.ident)
scRNAlist=list()
for (i in 1:6) {scRNAlist[[i]] <- subset(scRNAsub, idents = sample[i])}
#SCT标准化
scRNAlist <- lapply(X=scRNAlist, FUN=function(x) SCTransform(x))
##FindAnchors
scRNA.features <- SelectIntegrationFeatures(scRNAlist, nfeatures = 3000)
library(future)
options(future.globals.maxSize = 20 * 1024^3) #将全局变量上限调至20G
scRNAlist <- PrepSCTIntegration(scRNAlist, anchor.features = scRNA.features, verbose = FALSE)
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, normalization.method="SCT", 
                                        anchor.features=scRNA.features)
##Integrate
scRNAsub <- IntegrateData(anchorset = scRNA.anchors, normalization.method="SCT")
####添加细胞亚群分类信息
p_data <- rbind(cell_Myeloid,cell_VSMC)
scRNAsub <- AddMetaData(scRNAsub,p_data,col.name = c('celltype',"groups"))
###进行细胞亚群的cell-talk
table(scRNAsub@meta.data$celltype)
save(scRNAsub,file="Nichenetr/scRNAsub_VSMCandMyeloid.Rdata")
scRNAsub <- get(load("Nichenetr/scRNAsub_VSMCandMyeloid.Rdata"))
Idents(scRNAsub) <- "celltype"
DefaultAssay(scRNAsub) <- "integrated"
nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = scRNAsub, 
                                               top_n_ligands = 40,
                                               receiver = "Macrophage-like",
                                               sender = c("TREM2-high","Resident-like","Inflammatory","Contractile",
                                                          "Chondrocyte-like","Fibroblast-like"),
                                               condition_colname = "groups", 
                                               condition_oi = "AC", 
                                               condition_reference = "PA",
                                               ligand_target_matrix = ligand_target_matrix, 
                                               lr_network = lr_network, 
                                               weighted_networks = weighted_networks, 
                                               organism = "human",assay_oi = "integrated")
# top_n_ligands参数指定用于后续分析的高活性配体的数量                  
## 查看配体活性分析结果
# 主要参考pearson指标，bona_fide_ligand=True代表有文献报道的配体-受体，
# bona_fide_ligand=False代表PPI预测未经实验证实的配体-受体。
x <- nichenet_output$ligand_activities
write.csv(x, "Nichenetr/macrophage_ligand_activities.csv", row.names = F)
# 查看top20 ligands
nichenet_output$top_ligands
# 查看top20 ligands在各个细胞亚群中表达情况
pdf("Nichenetr/04_macrophage_ligand_split.pdf",10,8)
DotPlot(scRNAsub, features = nichenet_output$top_ligands, cols = "RdBu",split.by = "groups") +
  theme(legend.key.size=unit(0.4,'cm'),legend.key.height = unit(0.4,'cm'),
        legend.text = element_text(size = 8),legend.title =element_text(size = 12) )+
  RotatedAxis()
dev.off()
pdf("Nichenetr/05_macrophage_ligand.pdf",10,4)
DotPlot(scRNAsub, features = nichenet_output$top_ligands, cols = "RdBu") +
  theme(legend.key.size=unit(0.4,'cm'),legend.key.height = unit(0.4,'cm'),
        legend.text = element_text(size = 8),legend.title =element_text(size = 12) )+
  RotatedAxis()
dev.off()
####Macrophage-like VSMCs receptor
pdf("Nichenetr/06_macrophage_receptor.pdf",10,7)
nichenet_output$ligand_receptor_heatmap + 
  xlab("Receptors Expressed in Macrophage-like VSMCs")+
  ylab("Prioritized ligands modulated Macrophage-like VSMCs response")+
  theme(legend.key.size=unit(0.8,'cm'),legend.key.height = unit(0.4,'cm'),
        legend.text = element_text(size = 8),legend.title =element_text(size = 12) )
dev.off()
######Macrophage-like VSMCs Target
pdf("Nichenetr/07_macrophage_Target.pdf",21,7)
nichenet_output$ligand_target_heatmap + 
  xlab("Target Genes Expressed in Macrophage-like VSMCs")+
  ylab("Prioritized Ligands")+
  theme(legend.key.size=unit(0.8,'cm'),legend.key.height = unit(0.4,'cm'),
        legend.text = element_text(size = 8),legend.title =element_text(size = 12) )
dev.off()

pdf("Nichenetr/08_dotplot_macrophage_Target.pdf",21,7)
DotPlot(scRNAsub %>% subset(idents=c("Macrophage-like")), 
        features = nichenet_output$top_targets[1:50],split.by = "groups",cols = "RdBu",assay = "SCT") +
  theme(axis.text.x =element_text(size = 8)  ,legend.key.size=unit(0.4,'cm'),legend.key.height = unit(0.4,'cm'),
        legend.text = element_text(size = 8),legend.title =element_text(size = 12) )+
  RotatedAxis()+ 
  xlab("Target 50 Genes Expressed in Macrophage-like VSMCs")
dev.off()

###SCENIC重要的靶基因
deg.scenic <-  read_csv("Subset_VSMC/SCENIC/DEG-Reg_VSMC.csv") %>% 
  mutate(TF=str_match(.$Regulon,pattern = "[A-Za-z0-9]*"))
SCENIC_gene <- read_tsv("Subset_VSMC/SCENIC/output/Step2_regulonTargetsInfo.tsv") %>% 
  subset(.,TF==as.character(deg.scenic$TF))
###monocle和SCENIC和交集
a=intersect(beam$gene_short_name,SCENIC_gene$gene)##n=25
b=intersect(SCENIC_gene$gene,nichenet_output$top_targets)
intersect(a,b)##得到7个基因
####dorothea分析
###转录因子和nichenet交集
a=intersect(Target_and_TF,nichenet_output$top_targets)
###转录因子和monocle
b=intersect(Target_and_TF,beam$gene_short_name)
intersect(a,b)##得到12个基因
p=FeaturePlot(scRNAsub, features="VCAM1", label=T
            ,reduction = 'umap')
ggsave("Nichenetr/feature.png", p, width = 8, height = 4)

grep("CD34",a)





