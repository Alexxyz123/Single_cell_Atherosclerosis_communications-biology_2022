library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
dir.create("DEA_FUN")
rm(list=ls())

##===单细胞水平差异分析===##
load("scRNA_harmony.Rdata")
DimPlot(scRNA, reduction = 'umap', group.by = 'SCT_snn_res.0.4', label = T)
table(scRNA$orig.ident)

##提取需要分析的细胞（正常组vs.疾病组）
#提取某一细胞类型
table(scRNA@meta.data$SCT_snn_res.0.3)
cluster <- scRNA@meta.data %>% filter(SingleR=="T_cells")
table(cluster$SCT_snn_res.0.3)
scRNAsub <- subset(scRNA, idents = c(0, 2, 13))
cells1 <- subset(scRNA@meta.data, seurat_clusters %in% c(0, 2, 13) &
                   orig.ident %in% c("PA_01", "PA_02", "PA_03")) %>% rownames()
cells2 <- subset(scRNA@meta.data, seurat_clusters %in% c(0, 2, 13) &
                   orig.ident %in% c("AC_01", "AC_02", "AC_03")) %>% rownames()
##找出差异基因####
deg <- FindMarkers(scRNA, ident.1 = cells1, ident.2 = cells2, assay = 'SCT')
deg <- data.frame(gene = rownames(deg), deg)
deg <- filter(deg, p_val_adj<0.05) %>% arrange(desc(avg_logFC))
write.csv(deg, "DEA_FUN/DEG_pbmc_mono.csv", row.names = F)

##参数测试
if(F){
  deg1 <- FindMarkers(scRNA, ident.1 = cells1, ident.2 = cells2, assay = 'SCT')
  deg1 <- data.frame(gene = rownames(deg1), deg1)
  
  deg2 <- FindMarkers(scRNA, ident.1 = cells1, ident.2 = cells2, assay = 'SCT',
                      slot = 'counts', test.use = "DESeq2")
  deg2 <- data.frame(gene = rownames(deg2), deg2)
  
  deg3 <- FindMarkers(scRNA, ident.1 = cells1, ident.2 = cells2, assay = 'RNA',
                      slot = 'counts', test.use = "DESeq2")
  deg3 <- data.frame(gene = rownames(deg3), deg3)
  
  deg4 <- FindMarkers(scRNA, ident.1 = cells1, ident.2 = cells2, assay = 'SCT',
                      slot = 'data', test.use = "DESeq2")
  deg4 <- data.frame(gene = rownames(deg4), deg4)
  
  #qvalue过滤
  deg1 <- filter(deg1, p_val_adj<0.05)
  deg2 <- filter(deg2, p_val_adj<0.05)
  deg3 <- filter(deg3, p_val_adj<0.05)
  deg4 <- filter(deg4, p_val_adj<0.05)
  
  #venn图对比
  library(VennDiagram)
  p <- venn.diagram(list(deg1=deg1$gene, deg2=deg2$gene, deg3=deg3$gene, deg4=deg4$gene),
                    col = "transparent",
                    fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
                    cat.col = "darkblue",
                    category.names=c("deg1", "deg2", "deg3", "deg4"), 
                    alpha = 0.5, 
                    rotation.degree = 0,
                    cat.cex = 2.5, 
                    cat.dist = c(0.20, 0.20, 0.1, 0.1), 
                    cex =1.5, 
                    filename = NULL,
                    margin=0.2)    
  png(file="DEA_FUN/Venn_deg_test4.png", width=1000, height=1000)
  grid.draw(p)
  dev.off()
  
  p <- venn.diagram(list(deg2=deg2$gene, deg3=deg3$gene, deg4=deg4$gene),
                    col = "transparent",
                    fill = c("cornflowerblue", "green", "darkorchid1"),
                    cat.col = "darkblue",
                    category.names=c("deg2", "deg3", "deg4"), 
                    alpha = 0.5, 
                    rotation.degree = 0,
                    cat.cex = 2.5, 
                    cat.dist = c(0.08, 0.03, 0.05), 
                    cex =1.5, 
                    filename = NULL,
                    margin=0.2)    
  png(file="DEA_FUN/Venn_deg_test3.png", width=1000, height=1000)
  grid.draw(p)
  dev.off()
}
genelist <- pull(deg, gene) %>% as.character()

## GO&KEGG富集分析#####
#GO
ego_ALL <- enrichGO(gene          = genelist,
                    universe      = rownames(scRNA@assays$SCT@data),
                    OrgDb         = 'org.Hs.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego_all <- data.frame(ego_ALL)
write.csv(ego_all,'DEA_FUN/GO_pbmc_mono.csv', row.names = F)           
ego_CC <- enrichGO(gene          = genelist,
                   universe      = rownames(scRNA@assays$SCT@data),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_MF <- enrichGO(gene          = genelist,
                   universe      = rownames(scRNA@assays$SCT@data),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_BP <- enrichGO(gene          = genelist,
                   universe      = rownames(scRNA@assays$SCT@data),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
##截取字符串1:70，不要太长
ego_CC@result$Description <- substring(ego_CC@result$Description,1,70)
ego_MF@result$Description <- substring(ego_MF@result$Description,1,70)
ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)
#画图
p_BP <- barplot(ego_BP, showCategory = 10) + ggtitle("Biological process")
p_CC <- barplot(ego_CC, showCategory = 10) + ggtitle("Cellular component")
p_MF <- barplot(ego_MF, showCategory = 10) + ggtitle("Molecular function")
plotc <- p_BP/p_CC/p_MF
ggsave('DEA_FUN/GO_pbmc_mono.png', plotc, width = 12,height = 10)

#KEGG，先把基因名转为ENTREZID#####
genelist.k <- bitr(genelist, fromType="SYMBOL", toType="ENTREZID", 
                   OrgDb='org.Hs.eg.db') %>% pull(ENTREZID)
ekegg <- enrichKEGG(gene = genelist.k, organism = 'hsa')
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
df.kegg <- data.frame(ekegg)
write.csv(df.kegg, 'DEA_FUN/KEGG_pbmc_mono.csv', row.names = F)
p1 <- barplot(ekegg, showCategory=20)
p2 <- dotplot(ekegg, showCategory=20)
plotc = p1/p2
ggsave('DEA_FUN/KEGG_pbmc_mono.png', plot = plotc, width = 12, height = 10)

## GSEA富集分析######
msig.set <- read.gmt("Resource/MSigDB-Entrez/c2.cp.kegg.v7.2.entrez.gmt")
gene.use <- bitr(genelist, fromType="SYMBOL", toType="ENTREZID", OrgDb='org.Hs.eg.db')
names(gene.use) <- c("gene", "ENTREZID")
deg.use <- left_join(gene.use, deg)
gesa.genelist <- structure(deg.use$avg_logFC, names=deg.use$ENTREZID)
gsea <- GSEA(gesa.genelist, TERM2GENE = msig.set, minGSSize = 10,
             maxGSSize = 500, pvalueCutoff = 0.05, verbose=FALSE)
gsea <- setReadable(gsea, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
df.gsea <- data.frame(gsea)
write.csv(df.gsea, 'DEA_FUN/GSEA_pbmc_mono_kegg.csv', row.names = F)

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


##===pseudobulk RNA差异分析===##
library(Seurat)
library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(patchwork)
dir.create("bulkRNA")
rm(list=ls())
source("Resource/function.R")

load("scRNA_classify.Rdata")
DimPlot(scRNA, reduction = 'umap', group.by = 'seurat_clusters', label = T)
table(scRNA$orig.ident)
colnames(scRNA@meta.data)

##提取数据
Idents(scRNA) <- "seurat_clusters"
s.samples <- c("HNC01TIL", "HNC10TIL", "HNC20TIL", "HNC01PBMC", "HNC10PBMC", "HNC20PBMC")
s.idents <- c(3,9,15,17)
bulk.counts <- pseudobulk(scRNA, idents = s.idents, samples = s.samples)
write.csv(data.frame(Symbol=rownames(bulk.counts), bulk.counts), 
          file = "bulkRNA/pseudobulk_counts.csv", row.names = F)
head(bulk.counts)

##差异分析
data <- data.frame(bulk.counts, stringsAsFactors=F, check.names=F)
coldata <- data.frame(group = factor(c('TIL', 'TIL', 'TIL', 'PBMC', 'PBMC', 'PBMC')))
#构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~group)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- DESeq(dds)
res <- results(dds, contrast = c('group', 'TIL', 'PBMC'), pAdjustMethod = 'fdr')
res <- as.data.frame(res)
DEG <- merge(data, res, by=0)
names(DEG)[1]='gene'
DEG <- arrange(DEG, pvalue, desc(log2FoldChange))
write.csv(DEG, "bulkRNA/DEG.csv", row.names = F)
DEG.sig <- filter(DEG, pvalue < 0.05 & log2FoldChange > 2)
write.csv(DEG.sig, "bulkRNA/DEG-sig.csv", row.names = F)

##准备差异基因列表
genelist <- arrange(DEG.sig, desc(log2FoldChange)) %>% pull(gene) %>% as.character()

## GO&KEGG富集分析
#GO
ego_ALL <- enrichGO(gene          = genelist,
                    universe      = rownames(bulk.counts),
                    OrgDb         = 'org.Hs.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego_all <- data.frame(ego_ALL)
write.csv(ego_all,'bulkRNA/GO_pbmc_mono.csv', row.names = F)           
ego_CC <- enrichGO(gene          = genelist,
                   universe      = rownames(bulk.counts),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_MF <- enrichGO(gene          = genelist,
                   universe      = rownames(bulk.counts),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_BP <- enrichGO(gene          = genelist,
                   universe      = rownames(bulk.counts),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)           
ego_CC@result$Description <- substring(ego_CC@result$Description,1,70)
ego_MF@result$Description <- substring(ego_MF@result$Description,1,70)
ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)
p_BP <- barplot(ego_BP, showCategory = 10) + ggtitle("Biological process")
p_CC <- barplot(ego_CC, showCategory = 10) + ggtitle("Cellular component")
p_MF <- barplot(ego_MF, showCategory = 10) + ggtitle("Molecular function")
plotc <- p_BP/p_CC/p_MF
ggsave('bulkRNA/GO_pbmc_mono.png', plotc, width = 12,height = 10)

#KEGG，先把基因名转为ENTREZID
genelist.k <- bitr(genelist, fromType="SYMBOL", toType="ENTREZID", 
                   OrgDb='org.Hs.eg.db') %>% pull(ENTREZID)
ekegg <- enrichKEGG(gene = genelist.k, organism = 'hsa')
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
df.kegg <- data.frame(ekegg)
write.csv(df.kegg, 'bulkRNA/KEGG_pbmc_mono.csv', row.names = F)
p1 <- barplot(ekegg, showCategory=20)
p2 <- dotplot(ekegg, showCategory=20)
plotc = p1/p2
ggsave('bulkRNA/KEGG_pbmc_mono.png', plot = plotc, width = 12, height = 10)

## GSEA富集分析
msig.set <- read.gmt("Resource/MSigDB-Entrez/c2.cp.kegg.v7.2.entrez.gmt")
gene.use <- bitr(genelist, fromType="SYMBOL", toType="ENTREZID", OrgDb='org.Hs.eg.db')
names(gene.use) <- c("gene", "ENTREZID")
deg.use <- left_join(gene.use, DEG.sig) %>% arrange(desc(log2FoldChange))
gesa.genelist <- structure(deg.use$log2FoldChange, names=deg.use$ENTREZID)
gsea <- GSEA(gesa.genelist, TERM2GENE = msig.set, minGSSize = 10,
             maxGSSize = 500, pvalueCutoff = 0.05, verbose=FALSE)
gsea <- setReadable(gsea, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
df.gsea <- data.frame(gsea)
write.csv(df.gsea, 'bulkRNA/GSEA_pbmc_mono_kegg.csv', row.names = F)

# 批量作图
dir.create('bulkRNA/Gsea_Plot')
for(i in seq_along(gsea@result$ID)){
  p <- gseaplot(gsea, geneSetID = i, title = gsea@result$ID[i], by = "runningScore")
  filename <- paste0('bulkRNA/Gsea_Plot/', gsea@result$ID[i], '.png')
  ggsave(filename = filename, p, width = 8, height = 4)
}
# 美化版
if(F){
  library(enrichplot)
  for(i in seq_along(gsea@result$ID)){
    p <- gseaplot2(gsea, geneSetID = i, title = gsea@result$ID[i])
    filename <- paste0('bulkRNA/Gsea_Plot/pretty_', gsea@result$ID[i], '.png')
    ggsave(filename = filename, p, width = 8, height = 5)
  }
}
