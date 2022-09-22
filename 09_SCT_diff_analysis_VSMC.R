##差异与富集分析##
library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
setwd("~/single-cell/data/Subset_VSMC")
dir.create("DEA_FUN")
rm(list=ls())
##一、单细胞水平差异分析######
#加载数据
load("scRNAsub_classify.Rdata")
#VSMC monocle数据
get(load("Monocle/mycds.Rdata"))

DimPlot(scRNAsub, reduction = 'umap', group.by = 'seurat_clusters', label = T)
table(scRNAsub$orig.ident)

##1. AC vs.PA平滑肌细胞组间差异 #####
cells_PA <- subset(scRNAsub@meta.data,
                 orig.ident %in% c("PA_01", "PA_02", "PA_03")) %>% rownames()
cells_AC <- subset(scRNAsub@meta.data,
                 orig.ident %in% c("AC_01", "AC_02", "AC_03")) %>% rownames()
deg <- FindMarkers(scRNAsub, ident.1 = cells_AC, ident.2 = cells_PA,
                   assay = 'SCT',logfc.threshold =0)
allDiff <- deg %>% filter(abs(avg_log2FC)>0.25) %>% 
  filter(p_val_adj<0.05)
#获得基因列表
gene <- rownames(allDiff)

## geneList 三部曲
## 1.获取基因logFC
geneList <- allDiff$avg_log2FC
## 2.命名
names(geneList) =rownames(allDiff)
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)

head(geneList)

library(dplyr)

### 选人，hallmarks，两列
library(msigdbr)
dd <- msigdbr(species = "Homo sapiens")
h_kegg <- dd %>% 
  filter(gs_cat == "H") %>% 
  select(gs_name,gene_symbol)
### 主程序GSEA
y <- GSEA(geneList,TERM2GENE =h_kegg)
yd <- as.data.frame(y)

### 看整体分布
library(ggplot2)
pdf("revision_SMC_diff_hallmarkers.pdf",7,8)
dotplot(y,showCategory=80,split=".sign")
dev.off()

#####绘图######
library(ggpubr)
logFC_cutoff <- 0.25 
deg$change <- if_else(deg$p_val>0.01,'stable',
                      if_else(deg$avg_logFC >logFC_cutoff,'AC',
                              if_else(deg$avg_logFC < -logFC_cutoff,'PA','stable')))
deg$log10PValue <- -log10(deg$p_val)
deg <- deg %>% arrange(desc(log10PValue))
deg$log10PValue[1:14] <-seq(300,314,1)
deg$symbol <- rownames(deg)
top  <-  deg %>% arrange(desc(avg_logFC)) %>% top_n(10,avg_logFC) %>% pull(symbol) %>% as.character()
low  <-  deg %>% arrange(desc(avg_logFC)) %>% top_n(-20,avg_logFC) %>% pull(symbol) %>% as.character()

pdf(file="07_DEG_volcano.pdf", width=5, height=5)
ggscatter(deg, x = "avg_logFC", y = "log10PValue",  
          color = "change",size = 0.5,
          label = "symbol", repel = T, #展示差异倍数较大的基因
          label.rectangle = F,label.select = c(top,low),
          palette = c("#E41A1C","#377EB8" ,"#999999"),
          ylab = "-log10p.value")+ylim(0,350)+xlim(-4,1.5)
dev.off()
######得到差异基因####
deg <- data.frame(gene = rownames(deg), deg)
deg <- filter(deg, p_val_adj<0.05) %>%filter(.,abs(avg_logFC)>0.25) %>% arrange(desc(avg_logFC))
write.csv(deg, "DEA_FUN/DEG_AS_VSMC.csv", row.names = F)

######2.不同细胞类型中，macrophage-like的特征基因#####
Idents(scRNAsub) <- "celltype"
deg <- FindMarkers(scRNAsub, ident.1 = "Macrophage-like",
                   assay = 'SCT',logfc.threshold =0.25)

deg <- filter(deg, p_val_adj<0.05) %>% rownames_to_column(var = "gene")

genelist <- pull(deg, gene) %>% as.character()

## GO&KEGG富集分析
#GO
ego_ALL <- enrichGO(gene          = genelist,
                    universe      = rownames(scRNAsub@assays$SCT@data),
                    OrgDb         = 'org.Hs.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego_all <- data.frame(ego_ALL)
write.csv(ego_all,'DEA_FUN/GO_AS_VSMC.csv', row.names = F)           
ego_CC <- enrichGO(gene          = genelist,
                   universe      = rownames(scRNAsub@assays$SCT@data),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_MF <- enrichGO(gene          = genelist,
                   universe      = rownames(scRNAsub@assays$SCT@data),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_BP <- enrichGO(gene          = genelist,
                   universe      = rownames(scRNAsub@assays$SCT@data),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)           
ego_CC@result$Description <- substring(ego_CC@result$Description,1,70)
ego_MF@result$Description <- substring(ego_MF@result$Description,1,70)
ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)
ego_CC@result <- ego_CC@result %>% arrange(pvalue)
ego_MF@result <- ego_MF@result %>% arrange(pvalue)
ego_BP@result <- ego_BP@result %>% arrange(pvalue)


p_BP <- barplot(ego_BP, showCategory = 10) + ggtitle("Biological process")+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
p_CC <- barplot(ego_CC, showCategory = 10) + ggtitle("Cellular component")+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
p_MF <- barplot(ego_MF, showCategory = 10) + ggtitle("Molecular function")+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
plotc <- p_BP/p_CC/p_MF
ggsave('DEA_FUN/GO_AS_VSMC.png', plotc, width = 12,height = 10)
p_All <- barplot(ego_ALL, showCategory = 10) + ggtitle("Molecular function")
ggsave('DEA_FUN/GO_Macrophage_VSMC.png', p_All, width = 12,height = 10)
####pdf Figure
p_BP <- barplot(ego_BP, showCategory = 10,font.size = 8) + ggtitle("Biological process")+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
p_CC <-barplot(ego_CC, showCategory = 10,font.size = 8) + ggtitle("Cellular component")+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
p_MF <- barplot(ego_MF, showCategory = 10,font.size = 8) + ggtitle("Molecular function")+
  theme(plot.title = element_text(hjust = 0.5))  #也就加上这一行
plotc <- p_BP/p_CC/p_MF
pdf("08_Macrophage_GO.pdf",8,10)
plotc
dev.off()

#KEGG，先把基因名转为ENTREZID
genelist.k <- bitr(genelist, fromType="SYMBOL", toType="ENTREZID", 
                   OrgDb='org.Hs.eg.db') %>% pull(ENTREZID)
ekegg <- enrichKEGG(gene = genelist.k, organism = 'hsa')
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
df.kegg <- data.frame(ekegg)
write.csv(df.kegg, 'DEA_FUN/KEGG_AS_VSMC.csv', row.names = F)
p1 <- barplot(ekegg, showCategory=20)
p2 <- dotplot(ekegg, showCategory=20)
plotc = p1/p2
ggsave('DEA_FUN/KEGG_AS_VSMC.png', plot = plotc, width = 12, height = 10)
pdf("09_Macrophage_KEGG.pdf",8,6)
dotplot(ekegg, showCategory=20,font.size = 8)
dev.off()



## GSEA富集分析
##可以选择其他的gmt
msig.set <- read.gmt("/home/data/vip07/single-cell/data/Resource/MSigDB-Entrez/c7.all.v7.2.entrez.gmt")
gene.use <- bitr(genelist, fromType="SYMBOL", toType="ENTREZID", OrgDb='org.Hs.eg.db')
names(gene.use) <- c("gene", "ENTREZID")
deg.use <- left_join(gene.use, deg)
gesa.genelist <- structure(deg.use$avg_logFC, names=deg.use$ENTREZID)
gsea <- GSEA(gesa.genelist, TERM2GENE = msig.set, minGSSize = 10,
             maxGSSize = 500, pvalueCutoff = 0.1, verbose=FALSE)
gsea <- setReadable(gsea, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
df.gsea <- data.frame(gsea)
write.csv(df.gsea, 'DEA_FUN/GSEA_AS_VSMC_kegg.csv', row.names = F)

# 批量作图
dir.create('DEA_FUN/Gsea_Plot')
for(i in seq_along(gsea@result$ID)){
  p <- gseaplot(gsea, geneSetID = i, title = gsea@result$ID[i], by = "runningScore")
  filename <- paste0('DEA_FUN/Gsea_Plot/', gsea@result$ID[i], '.png')
  ggsave(filename = filename, p, width = 8, height = 4)
}
# 美化版
if(F){
  library(enrichplot)
  for(i in seq_along(gsea@result$ID)){
    p <- gseaplot2(gsea, geneSetID = i, title = gsea@result$ID[i])
    filename <- paste0('DEA_FUN/Gsea_Plot/pretty_', gsea@result$ID[i], '.png')
    ggsave(filename = filename, p, width = 8, height = 5)
  }
}

#####monocle，不同state的差异分析#####
library(monocle)
library(Seurat)
#VSMC monocle数据
get(load("Monocle/mycds.Rdata"))
###获得表型信息
p_data <- subset(phenoData(mycds)@data,select='State')
scRNAsub <- AddMetaData(scRNAsub,p_data,col.name = 'State')
####查看PA组中的macrophage-like细胞在state4中一共有多少个







###读取PA vs. AC之间的差异基因
DEG <- read_csv("DEA_FUN/DEG_AS_VSMC.csv")
intersection=intersect(DEG$gene,rownames(sig_dge.State))
write.csv(file = "monocle_diff_gene.csv",intersection)
DimPlot(scRNAsub, reduction = 'umap', group.by = 'State', label = T)

DimPlot(scRNAsub, reduction = 'umap', group.by = 'celltype', label = T,
        split.by = "groups")




