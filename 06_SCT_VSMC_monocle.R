library(monocle)
library(tidyverse)
library(patchwork)
setwd("~/single-cell/data/Subset_VSMC")
dir.create("Monocle")
##===monocle===##
# monocle不能基于整合的数据进行分析，个人不建议用于多样本分析
rm(list = ls())
######创建monocle分析对象####
load("scRNAsub_classify.Rdata")
scRNAsub$groups <- scRNAsub$orig.ident
scRNAsub$groups <- recode(scRNAsub$groups, 
                          'AC_01' = "AC", 
                          'PA_01' = "PA",
                          'AC_02'= "AC",
                          'PA_02' = "PA",
                          'AC_03' = "AC",
                          'PA_03' = "PA")
Idents(scRNAsub) <- "orig.ident"
dim(scRNAsub)
# monocle不推荐使用slot="data"
data <- GetAssayData(scRNAsub, assay = "RNA", slot = "counts")
#data <- as(data, 'sparseMatrix')   #将普通矩阵转为稀疏矩阵
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
table(scRNAsub$seurat_clusters)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())
######数据预处理#######
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=16)
######选择排序基因
#确定需要排除与包含的基因
ex1 <- c(cc.genes$s.genes, cc.genes$g2m.genes)
ex2 <- grep('^MT-', rownames(scRNAsub), value = T)
ex3 <- grep('^RP[LS]', rownames(scRNAsub), value = T)
ex <- unique(c(ex1, ex2, ex3))
##seurat确定的高变基因
order.genes <- SCTransform(scRNAsub) %>% VariableFeatures()
#order.genes <- order.genes[!order.genes %in% ex]
mycds <- setOrderingFilter(mycds, order.genes)
plot_ordering_genes(mycds)
#或者由monocle选择高变基因
if(F){
  disp_table <- dispersionTable(mycds)
  order.genes1 <- subset(disp_table, mean_expression >= 0.01 &
                          dispersion_empirical >= 1 * dispersion_fit) %>%
    pull(gene_id) %>% as.character()
  order.genes1 <- order.genes1[!order.genes1 %in% ex]
  order.genes1 <- intersect(order.genes,order.genes1)
  mycds <- setOrderingFilter(mycds, order.genes1)
  plot_ordering_genes(mycds)
}
######降维####
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#######排序
mycds <- orderCells(mycds)
###修改pseudotime的顺序
mycds <- orderCells(mycds,reverse = T)
####结果可视化####
#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State",cell_size = 0.5,
                              palette=color_manual)
ggsave("Monocle/Trajectory_State.png", plot = plot1, width = 6, height = 5)
pdf("Monocle/01_Trajectory_State.pdf",4,4)
plot_cell_trajectory(mycds, color_by = "State",cell_size = 0.5)
dev.off()
#Celltype轨迹分布图
plot2 <- plot_cell_trajectory(mycds, color_by = "celltype",cell_size = 0.5)
ggsave("Monocle/Trajectory_Celltype.png", plot = plot2, width = 6, height = 5)
pdf("Monocle/02_Trajectory_Celltype.pdf",4,4)
plot_cell_trajectory(mycds, color_by = "celltype",cell_size = 0.001)+
  facet_wrap(~celltype, nrow = 2)+NoLegend()+
  theme(legend.key.size=unit(0.5,'cm'),legend.key.height = unit(0.5,'cm'))
dev.off()
#Pseudotime轨迹图
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size = 0.5)
ggsave("Monocle/Trajectory_Pseudotime.png", plot = plot3, width = 6, height = 5)
pdf("Monocle/03_Trajectory_pseudotime.pdf",4,4)
plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size = 0.5)+
  theme(legend.key.size=unit(0.5,'cm'),legend.key.height = unit(0.4,'cm'))
dev.off()
#group轨迹分布图
plot4 <- plot_cell_trajectory(mycds, color_by = "State",cell_size = 0.5)+
  facet_wrap(~groups)
ggsave("Monocle/Trajectory_groups.png", plot = plot4, width = 6, height = 5)
pdf("Monocle/04_Trajectory_group.pdf",8,4)
plot_cell_trajectory(mycds,color_by = "State",cell_size = 0.5)+
  facet_wrap(~groups)+
  theme(legend.key.size=unit(0.5,'cm'),legend.key.height = unit(0.5,'cm'))
dev.off()

#合并作图
plotc <- plot1|plot2|plot3|plot4
ggsave("Monocle/Trajectory_Combination.png", plot = plotc, width = 10, height = 3.5)
##保存结果
save(mycds,order.genes, file = "Monocle/mycds.Rdata")
####获得之前的结果####
get(load("Monocle/mycds.Rdata"))

##指定基因的可视化
# s.genes <- c("MYH11","ATCA2","TAGLN","ACTA2",
#              "MYOCD","CD68","LGALS3","CXCR4",
#              "CD74","PTPRC",
#              "MYH10","ELN","FN1",
#              "PDGFRB","RUNX2","BMP2")
# s.genes <- CaseMatch(s.genes, rownames(scRNAsub))
# p1 <- plot_genes_jitter(mycds[s.genes,], grouping = "celltype", color_by = "celltype")
# p2 <- monocle::plot_genes_violin(mycds[s.genes,], grouping = "celltype", color_by = "celltype")
# p3 <- monocle::plot_genes_in_pseudotime(mycds[s.genes,], color_by = "celltype")
# plotc <- p1|p2|p3
# ggsave("Monocle/Genes_Jitterplot.png", plot = plotc, width = 9, height = 18)
####拟时轨迹分析-热图######
##寻找拟时差异基因。不要使用多核运算，可能警告：关闭不再使用的链结
#拟时差异基因分析，花一点时间不长
Time_diff <- differentialGeneTest(mycds[order.genes,], cores = 8,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)")
Time_diff <- Time_diff[,c(5,2,3,4,1,6)]
write.csv(Time_diff, "Monocle/Time_diff_all.csv", row.names = F)
#显著差异基因的作图
Time_genes <- top_n(Time_diff, n = -50, qval) %>%
  pull(gene_short_name) %>% as.character()
# p = plot_pseudotime_heatmap(mycds[Time_genes,], num_clusters=3,
#                             show_rownames=T, return_heatmap=T)
# ggsave("Monocle/Time_heatmap.png", p, width = 5, height = 10)
# #显著差异基因按热图结果排序并保存
hp.genes <- p$tree_row$labels[p$tree_row$order]
Time_diff_sig <- Time_diff[hp.genes, c("gene_short_name", "pval", "qval")]
write.csv(Time_diff_sig, "Monocle/Time_diff_sig.csv", row.names = F)
## BEAM分析(找到与分叉点相关的基因),我们这儿是第2个Branch
#order.genes为高变基因
beam_res <- BEAM(mycds[order.genes,], branch_point = 2, cores = 16)
beam_res  <- beam_res %>% filter(status!="FAIL")
beam_res <- beam_res[order(beam_res$qval),]
#导出与brach_2明显相关基因
write.csv(beam_res, "Monocle/BEAM_all_2.csv", row.names = F)
BEAM_genes <- top_n(beam_res, n = -100, qval) %>%
  pull(gene_short_name) %>% as.character()
###绘制不同状态热图
##筛选q value< 1e-4(monocle标准流程)
p_total=plot_genes_branched_heatmap(mycds[row.names(subset(beam_res, qval < 1e-4)),],
                            branch_states = c(3,4),
                            branch_labels = c("State 3", "State 4"),
                            num_clusters = 3, show_rownames = F, 
                            return_heatmap = T)
ggsave("Monocle/Branch_heatmap_total.png", p_total$ph_res, width = 6.5, height = 8)
pdf("Monocle/05_Trajectory_heatmapall.pdf",6,12)
p_total$ph_res
dev.off()

p1 <- plot_genes_branched_heatmap(mycds[BEAM_genes,],
                                 branch_states = c(3,4),
                                 branch_labels = c("State 3", "State 4"),
                                 num_clusters = 3, show_rownames = T, 
                                 return_heatmap = T)+scale_fill_npg()
ggsave("Monocle/Branch_heatmap_100.png", p1$ph_res, width = 6.5, height = 8)
####monocle热图
p <- plot_genes_branched_heatmap(mycds[BEAM_genes,],
                            branch_states = c(3,4),
                            branch_labels = c("State 3", "State 4"),
                            num_clusters = 3, show_rownames = T, 
                            cluster_rows = F,
                            return_heatmap = T)
pdf("Monocle/05_Trajectory_heatmap100.pdf",8,12)
p$ph_res
dev.off()
####绘制高表达的BEAM_gene在不同State和celltype之间的差别
##state3_AC,State4_PA
State3 <- c("APOE","SPP1","CCL4")
State4 <- c("APOD","KLF4","DCN")
pdf("Monocle/06_Branch_state3.pdf",6,2)
plot_genes_branched_pseudotime(mycds[State3,],
                               branch_point=2,
                               color_by="State",
                               ncol=3)
dev.off()
pdf("Monocle/07_Branch_state4.pdf",6,2)
plot_genes_branched_pseudotime(mycds[State4,],
                               branch_point=2,
                               color_by="State",
                               ncol=1)+
  theme(legend.key.size=unit(0.5,'cm'),legend.key.height = unit(0.5,'cm'))
dev.off()
pdf("Monocle/07_Branch_state.pdf",7,3)
plot_genes_branched_pseudotime(mycds[c(State4,State3),],
                               branch_point=2,
                               color_by="State",
                               ncol=3)+
  theme(legend.key.size=unit(0.5,'cm'),legend.key.height = unit(0.5,'cm'))
dev.off()


p=plot_genes_branched_pseudotime(mycds[BEAM_genes[1:10],],
                               branch_point=2,
                               color_by="celltype",
                               ncol=1)
ggsave("Monocle/Branch_celltype.png", p, width = 6.5, height = 8)
#按照不同的状态进行绘制
p=plot_genes_branched_pseudotime(mycds[BEAM_genes[1:10],],
                                 branch_point=2,
                                 color_by="State",
                                 ncol=1)
ggsave("Monocle/Branch_state.png", p, width = 6.5, height = 8)
#按照分组进行绘制
p=plot_genes_branched_pseudotime(mycds[BEAM_genes[1:10],],
                                 branch_point=2,
                                 color_by="groups",
                                 ncol=1)
ggsave("Monocle/Branch_groups.png", p, width = 6.5, height = 8)
#显著差异基因按热图结果排序(100个)并保存
hp.genes_total <- p_total$ph_res$tree_row$labels[p_total$ph_res$tree_row$order]
BEAM_sig_total <- beam_res[hp.genes_total, c("gene_short_name", "pval", "qval")]
write.csv(BEAM_sig_total, "Monocle/BEAM_sig.csv", row.names = F)

hp.genes <- p1$ph_res$tree_row$labels[p1$ph_res$tree_row$order]

BEAM_sig <- beam_res[hp.genes, c("gene_short_name", "pval", "qval")]
rownames(BEAM_sig) <- NULL
BEAM_sig$change <- 1
BEAM_sig$change[1:60] <-"PA"
BEAM_sig$change[61:73] <- "Stable"
BEAM_sig$change[74:100] <- "AC"
write.csv(BEAM_sig, "Monocle/BEAM_sig.csv", row.names = F)
#####提取每个cluster的基因,进行通路分析KEGG#####
row_cluster=cutree(p_total$ph_res$tree_row,k=3)
gene_cluster1 <-names(row_cluster[which(row_cluster==1)])
gene_cluster2 <-names(row_cluster[which(row_cluster==2)])
gene_cluster3 <-names(row_cluster[which(row_cluster==3)])
#####gene_cluster1######
gene_cluster1.k <- bitr(gene_cluster1, fromType="SYMBOL", toType="ENTREZID", 
                        OrgDb='org.Hs.eg.db') %>% pull(ENTREZID)
ekegg <- enrichKEGG(gene = gene_cluster1.k, organism = 'hsa')
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
ekegg@result$Description<- substring(ekegg@result$Description,1,70)
df.kegg <- data.frame(ekegg)
write.csv(df.kegg, 'Monocle/KEGG_gene_cluster1.csv', row.names = F)
p1 <- dotplot(ekegg, showCategory=20,font.size = 12)+ ggtitle("Gene Cluste1")+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
#####gene cluster2#####
gene_cluster2.k <- bitr(gene_cluster2, fromType="SYMBOL", toType="ENTREZID", 
                        OrgDb='org.Hs.eg.db') %>% pull(ENTREZID)
ekegg <- enrichKEGG(gene = gene_cluster2.k, organism = 'hsa')
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
ekegg@result$Description<- substring(ekegg@result$Description,1,70)
df.kegg <- data.frame(ekegg)
write.csv(df.kegg, 'Monocle/KEGG_gene_cluster2.csv', row.names = F)

p2 <- dotplot(ekegg, showCategory=20,font.size = 12)+ ggtitle("Gene Cluster 2")+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
ggsave('Monocle/KEGG_gene_cluster2.png', plot = p2, width = 8, height = 6)
#####gene_cluster3#####
gene_cluster3.k <- bitr(gene_cluster3, fromType="SYMBOL", toType="ENTREZID", 
                        OrgDb='org.Hs.eg.db') %>% pull(ENTREZID)
ekegg <- enrichKEGG(gene = gene_cluster3.k, organism = 'hsa')
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
ekegg@result$Description<- substring(ekegg@result$Description,1,70)
df.kegg <- data.frame(ekegg)
write.csv(df.kegg, 'Monocle/KEGG_gene_cluster3.csv', row.names = F)

p3 <- dotplot(ekegg, showCategory=20,font.size = 12)+ ggtitle("Gene Cluster 3")+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
ggsave('Monocle/KEGG_gene_cluster3.png', plot = p3, width = 8, height = 6)
####拼图
p=(p1/p2)
ggsave('Monocle/KEGG_gene_cluster4.png', plot = p, width =10, height = 12)

pdf("Monocle/08_State_KEGG.pdf",10,12)
p1/p2
dev.off()
write.csv(gene_cluster1,file = "Monocle/gene_cluster1.csv")
write.csv(gene_cluster2,file = "Monocle/gene_cluster2.csv")
write.csv(gene_cluster3,file = "Monocle/gene_cluster3.csv")
########堆栈图，不同State之间进行分析####
p_data <- subset(phenoData(mycds)@data,select='State')
scRNAsub <- AddMetaData(scRNAsub,p_data,col.name = 'State')
tmp <-dplyr::select(scRNAsub@meta.data,c("groups", "State"))
df <- data.frame()
for(i in unique(tmp$groups)){
  df_i <- subset(tmp, tmp$groups==i) %>% pull(State) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("sample", "State", "value")
  df <- rbind(df, df_i)
}
source("/home/data/vip07/single-cell/data/Resource/function.R")
library(ggsci)
library(viridis)
library(scales)
#按样本统计细胞类型
color_state=c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3")
pdf(file = "Monocle/09_Stackbar_State.pdf",5,5)
ggplot(df, aes(x=sample, y=value, fill=State)) + 
  geom_bar(stat= "identity", position = "fill",width=0.4)+scale_y_continuous(expand=c(0,0))+
  labs(x = 'Group', y = 'Relative Abundance', title = 'Group composition') +
  theme_classic()+theme(legend.key.size=unit(0.4,'cm'))+ 
  scale_fill_manual(values = color_state)+
  guides(fill = guide_legend(title = 'State'))+#修改图例大小
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行

dev.off()
