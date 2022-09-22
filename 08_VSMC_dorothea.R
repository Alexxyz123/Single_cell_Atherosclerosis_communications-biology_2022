## We load the required packages
library(dorothea)
library(dplyr)
library(Seurat)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)
library(patchwork)
library(ggplot2)
setwd("/home/data/ssy008/single-cell/data/Subset_VSMC")
rm(list = ls())
load("scRNAsub_classify.Rdata")
#####加入State相关信息####
##添加monocle的state信息
get(load("Monocle/mycds.Rdata"))
###获得表型信息
p_data <- subset(phenoData(mycds)@data,select='State')
scRNAsub <- AddMetaData(scRNAsub,p_data,col.name = 'State')
scRNAsub$groups <- scRNAsub$orig.ident
scRNAsub$groups <- recode(scRNAsub$groups, 
                          'AC_01' = "AC", 
                          'PA_01' = "PA",
                          'AC_02'= "AC",
                          'PA_02' = "PA",
                          'AC_03' = "AC",
                          'PA_03' = "PA")
Idents(scRNAsub) <- "orig.ident"
## We read Dorothea Regulons for Human:
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
## We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C","D"))
## We compute Viper Scores 
scRNAsub <- run_viper(scRNAsub, regulon,
                  options = list(method = "scale", minsize = 4, 
                                 eset.filter = FALSE, cores = 32, 
                                 verbose = FALSE),assay_key = "SCT")
## We compute the Nearest Neighbours to perform cluster
DefaultAssay(object = scRNAsub) <- "dorothea"
scRNAsub <- ScaleData(scRNAsub)
scRNAsub <- RunPCA(scRNAsub, features = rownames(scRNAsub), verbose = FALSE)
p <- ElbowPlot(scRNAsub,ndims = 50)
ggsave(filename = "dorothea/elbowplot.png",height = 10,width =10)
scRNAsub <- FindNeighbors(scRNAsub, dims = 1:20, verbose = FALSE)
scRNAsub <- FindClusters(scRNAsub, resolution = c(seq(0.5,0.1,-0.1)), verbose = FALSE)
scRNAsub <- RunUMAP(scRNAsub, dims = 1:30, umap.method = "uwot", 
                    metric = "cosine")
DimPlot(scRNAsub)
######PA vs. AC组的的差异分析
##PA vs. AC组的的差异分析
deg_doro <- FindMarkers(scRNAsub, ident.1 ="AC", ident.2 ="PA",
                       assay = "dorothea", 
                       group.by = "groups",
                       logfc.threshold = 0.25, min.pct = 0.25)
deg_doro <- data.frame(TF = rownames(deg_doro), deg_doro)
deg_doro <- filter(deg_doro, p_val<0.05) %>% arrange(desc(avg_logFC))
deg_doro$change <- NULL
deg_doro$change <- if_else(deg_doro$avg_logFC>0,"AC","PA")
write.csv(deg.reg, "dorothea/DEG-Reg_VSMC.csv", row.names = F)
####AS vs. Normal的macrophage-like差异分析#####
cells_AC <- subset(scRNAsub_doro@meta.data, celltype=="Macrophage-like"&groups=="AC") %>% rownames()
cells_PA <- subset(scRNAsub_doro@meta.data, celltype=="Macrophage-like"&groups=="PA") %>% rownames()
deg_TF_AC_PA <- FindMarkers(scRNAsub_doro, ident.1 = cells_AC, ident.2 = cells_PA, 
                   logfc.threshold = 0.25, min.pct = 0.25)
deg_TF_AC_PA <- data.frame(TF = rownames(deg_TF_AC_PA), deg_TF_AC_PA)
deg_TF_AC_PA <- filter(deg_TF_AC_PA, p_val<0.05) %>% arrange(desc(avg_logFC))
deg_TF_AC_PA$change <- NULL
deg_TF_AC_PA$change <- if_else(deg_TF_AC_PA $avg_logFC>0,"AC","PA")
TF_deg <- intersect(deg$TF,intersect(deg_doro$TF,deg_doro_State$TF))

####不同state的差异分析
deg_doro_State <- FindMarkers(scRNAsub, ident.1 =3, ident.2 =4, assay = "dorothea", 
                             group.by = "State",
                             logfc.threshold = 0.25, min.pct = 0.25)
deg_doro_State <- data.frame(TF = rownames(deg_doro_State), deg_doro_State)
deg_doro_State <- filter(deg_doro_State, p_val<0.05) %>% arrange(desc(avg_logFC))
deg_doro_State$change <- NULL
deg_doro_State$change <- if_else(deg_doro_State$avg_logFC>0,"State3","State4")
write.csv(deg_doro_State, "dorothea/DEG-Reg_State.csv", row.names = F)
####macrophage-like和其他细胞类型的比较
cells1 <- subset(scRNAsub_doro@meta.data, celltype=="Macrophage-like") %>% rownames()
cells2 <- subset(scRNAsub_doro@meta.data, celltype!="Macrophage-like") %>% rownames()
deg <- FindMarkers(scRNAsub_doro, ident.1 = cells1, ident.2 = cells2, 
                   logfc.threshold = 0.25, min.pct = 0.25)
deg <- data.frame(TF = rownames(deg), deg)
deg <- filter(deg, p_val<0.05) %>% arrange(desc(avg_log2FC))
deg$change <- NULL
deg$change <- if_else(deg$avg_log2FC>0,"Macro-up","Macro-down")
Macorphage_deg <- deg$TF
## Assigning cell type identity to clusters
# new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", 
#                      "FCGR3A+ Mono", "NK", "DC", "Platelet")
# names(new.cluster.ids) <- levels(scRNAsub)
# scRNAsub <- RenameIdents(scRNAsub, new.cluster.ids)
## We transform Viper scores, scaled by seurat, into a data frame to better 
## handling the results
####保存结果
save(scRNAsub,file = "dorothea/scRNAsub_doro.Rdata")
#获得结果
scRNAsub_doro <- get(load("dorothea/scRNAsub_doro.Rdata"))
##得到筛选出来的TF，进行绘图
# pheno_data <- scRNAsub@meta.data[,c("groups","celltype")] %>%
#   filter(celltype=="Macrophage-like") %>% 
#   arrange(groups) %>% arrange(celltype)
# heat_data <- viper_scores_df[,highly_variable_tfs$tf] %>% t() %>% 
#   .[,rownames(pheno_data)]

#####macrophage-like细胞的热图######
viper_scores_df <- GetAssayData(scRNAsub_doro, slot = "scale.data", 
                                assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t()
## We create a data frame containing the cells and their clusters
Idents(scRNAsub_doro) <- "celltype"
CellsClusters <- data.frame(cell = names(Idents(scRNAsub_doro)), 
                            cell_type = as.character(Idents(scRNAsub_doro)),
                            check.names = F)
## We create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)
## We summarize the Viper scores by cellpopulation
summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))
## We select the 20 most variable TFs. (20*9 populations = 180)
highly_variable_tfs <- data.frame(tf=Macorphage_deg)
## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 
write_csv(summarized_viper_scores_df %>% t() %>% 
            as.data.frame() %>% rownames_to_column(var = "ID"),
          file = "Figure4A_summarized_viper_scores_df.csv")
summarized_viper_scores_df <- read.csv("Figure4A_summarized_viper_scores_df.csv") %>% 
  column_to_rownames(var = "ID")
#pheatmap
palette_length = 100
my_color = colorRampPalette(c("darkblue", "white","red"))(palette_length)
my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))
colnames(summarized_viper_scores_df) <- c("Chondrocyte-like",
                             "Contractile","Fibroblast-like",
                             "Macrophage-like"
                             )

names <- c("Macrophage-like","Contractile","Fibroblast-like","Chondrocyte-like")

pdf("dorothea/01_Macrophage_TF_heatmap.pdf",6,8)
pheatmap(summarized_viper_scores_df[,names],fontsize=10, 
                       fontsize_row = 5, 
                       color=my_color, breaks = my_breaks, 
                       angle_col = 45,
                       treeheight_col = 0,  border_color = NA) 

dev.off()
#######macrophage-like和其他细胞类型进行取交集，韦恩图#####
a=deg %>% filter(change=="Macro-up") %>% pull(TF)
b <- deg_TF_AC_PA %>% filter(change=="AC") %>% pull(TF)
c=deg %>% filter(change=="Macro-down") %>% pull(TF)
d <- deg_TF_AC_PA %>% filter(change=="PA") %>% pull(TF)
intersect(a,d)
intersect(c,b)
write.csv(a,file = "dorothea/Up-regulated in Macrophage-like.csv")
write.csv(c,file = "dorothea/down-regulated in Macrophage-like.csv")
write.csv(b,file = "dorothea/Up-regulated in AC.csv")
write.csv(d,file = "dorothea/down-regulated in AC.csv")
library(VennDiagram)
p <- venn.diagram(list(a, b, c, d),
                  fill = c("#00468B", "#ED0000", "#42B540", "#0099B4"),
                  cat.col = "darkblue",
                  category.names=c("Up-regulated in Macrophage-like", 
                                   "Up-regulated in AC", 
                                   "down-regulated in Macrophage-like", 
                                   "down-regulated in AC"), 
                  lwd = 3,
                  col = "black",
                  alpha = 0.4, 
                  rotation.degree = 0,
                  cat.cex = 1.2, 
                  sub.col = "black",main.col = "black",
                  cat.dist = c(0.20, 0.20, 0.1, 0.1), 
                  cex =1.5, 
                  filename = NULL,
                  margin=0.2)    
pdf(file="dorothea/02_Venn_deg_TF.pdf",11, 11)
grid.draw(p)
dev.off()
#######提取macrophage,ACvs.PA########
#获得分别在macrophage-like大类中上调以及在PA/AC中同向调控的TF
TF_macro <- c(intersect(a,b),intersect(a,d),intersect(c,b),intersect(c,d))
write.csv(TF_macro,file = "Subset_VSMC/dorothea/TF_macro.csv")
####转录因子评分矩阵
viper_scores_df <- GetAssayData(scRNAsub_doro, slot = "scale.data", 
                                assay = "dorothea") %>%.[TF_macro,] %>% 
  data.frame(check.names = F) 
## We create a data frame containing the cells and their clusters
Idents(scRNAsub_doro) <- "celltype"
CellsClusters <- scRNAsub_doro@meta.data[c(cells_AC,cells_PA),c("groups","State")]%>%
  arrange(State) %>%arrange(groups) 
viper_scores_df <- viper_scores_df[,rownames(CellsClusters)]
palette_length = 100
my_color =color <- colorRampPalette(c("skyblue","white","red"))(100)
my_breaks <- c(seq(min(viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(viper_scores_df)/palette_length, 
                   max(viper_scores_df), 
                   length.out=floor(palette_length/2)))
pdf("dorothea/03_Macro_group_TF.pdf",6,8)
pheatmap(viper_scores_df,fontsize=8, color = my_color,
                       fontsize_row = 6, show_colnames = F,show_rownames = T,
                       breaks = my_breaks, annotation_col = CellsClusters,
                       treeheight_col = 0,  border_color = NA,cluster_cols = F,
         cluster_rows = F) 
dev.off()
###绘制Featureplot
#默认assay
DefaultAssay(scRNAsub_doro) <- "dorothea"
###整合坐标系
##从seurat导入整合过的umap坐标
scRNA.embed <- Embeddings(scRNAsub_doro, reduction = "umap")
scRNAsub <- get(load("scRNAsub_classify.Rdata"))
int.embed <- Embeddings(scRNAsub, reduction = "umap")
scRNAsub_doro@reductions$umap@cell.embeddings <- int.embed
# #FeaturePlot
# p1 = FeaturePlot(scRNAsub_doro, features='TBX21', label=T, 
#                  cols = c("grey","firebrick3")
#                  ,reduction = 'umap')
# p2 = DimPlot(scRNAsub_doro, reduction = 'umap', group.by = "celltype", label=T)
# plotc = p1|p2
# ggsave('dorothea/KLF4.png', plotc, width=18 ,height=4)
#######RidgePlot&VlnPlot
data <- scRNAsub_doro@assays$dorothea@scale.data %>% t() %>% as.data.frame() %>% 
  rownames_to_column(.,var = "ID")
tmp <- scRNAsub_doro@meta.data %>% as.data.frame()%>% rownames_to_column(.,var = "ID")
data <- merge(data,tmp,by="ID")
data$celltype <-  relevel(data$celltype,ref = "Macrophage-like") 

p1 <- ggplot(data = data, aes(x=celltype, y=JUN,fill=groups))+
  geom_boxplot(alpha = 1,width=0.4,outlier.size = 0.4) + #半透明
  theme_classic() + #或theme_bw()
  scale_fill_brewer(palette = "Set1") + #按类填充颜色
  scale_color_brewer(palette = "Set1") + #按类给边框着色
  theme(axis.text.x = element_text(colour="black", size = 11,
                                   #癌症名太挤，旋转45度
                                   angle = 45, hjust = .5, vjust = .5)) +
  ylab( paste("JUN", "expression"))+ggtitle("JUN")+
  theme(plot.title = element_text(hjust = 0.5))  +  #也就加上这一行
  theme(axis.title.x=element_blank(),axis.text.x = element_blank())+NoLegend()

p2 = ggplot(data = data, aes(x=celltype, y=SMAD3,fill=groups))+
  geom_boxplot(alpha = 1,width=0.4,outlier.size = 0.4) + #半透明
  theme_classic() + #或theme_bw()
  scale_fill_brewer(palette = "Set1") + #按类填充颜色
  scale_color_brewer(palette = "Set1") + #按类给边框着色
  theme(axis.text.x = element_text(colour="black", size = 11,
                                   #癌症名太挤，旋转45度
                                   angle = 45, hjust = .5, vjust = .5)) +
  ylab( paste("SMAD3", "expression"))+ggtitle("SMAD3")+
  theme(plot.title = element_text(hjust = 0.5))  +  #也就加上这一行
  theme(axis.title.x=element_blank(),axis.text.x = element_blank())

p3 = ggplot(data = data, aes(x=celltype, y=CBX2,fill=groups))+
  geom_boxplot(alpha = 1,width=0.4,outlier.size = 0.4) + #半透明
  theme_classic() + #或theme_bw()
  scale_fill_brewer(palette = "Set1") + #按类填充颜色
  scale_color_brewer(palette = "Set1") + #按类给边框着色
  
  theme(axis.text.x = element_text(colour="black", size = 11,
                                   #癌症名太挤，旋转45度
                                   angle = 45, hjust = .5, vjust = .5)) +
  ylab( paste("CBX2", "expression"))+ggtitle("CBX2")+
  theme(plot.title = element_text(hjust = 0.5))  +  #也就加上这一行
  theme(axis.title.x=element_blank())+NoLegend()
p4 = ggplot(data = data, aes(x=celltype, y=ZHX2,fill=groups))+
  geom_boxplot(alpha = 1,width=0.4,outlier.size = 0.4) + #半透明
  theme_classic() + #或theme_bw()
  scale_fill_brewer(palette = "Set1") + #按类填充颜色
  scale_color_brewer(palette = "Set1") + #按类给边框着色
  theme(axis.text.x = element_text(colour="black", size = 11,
                                   #癌症名太挤，旋转45度
                                   angle = 45, hjust = .5, vjust = .5)) +
  ylab( paste("ZHX2", "expression"))+ggtitle("ZEB2")+
  theme(plot.title = element_text(hjust = 0.5))+  #也就加上这一行
  theme(axis.title.x=element_blank())
  

plotc = (p1|p2)/(p3|p4)

pdf("dorothea/04_boxplot_TF.pdf",8,5)
plotc
dev.off()
######转录因子以及其靶基因,构建网络#######
Target_and_TF <- subset(regulon,tf%in%TF_macro)
write.csv(Target_and_TF,file = "dorothea/Target_and_TF_network.csv")
TF <- data.frame(gene=TF_macro,class='TF')
Target <- data.frame(gene=unique(Target_and_TF$target),class="gene")
gene <- rbind(TF,Target)
write.csv(gene,file = "dorothea/gene.csv")
TF_upanddown <- subset(deg_TF_AC_PA,TF%in%TF_macro)
write.csv(TF_upanddown,file = "dorothea/TF_upanddowne.csv")
