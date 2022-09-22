library(tidyverse)
######一、大类细胞的singleR注释#######
#读取数据
scRNA$seurat_clusters=scRNA$integrated_snn_res.0.4

cluster_marker <- read_csv("CellType/ClusterMarker.csv")
deg_4 <- cluster_marker %>% filter(cluster==4) %>% 
  pull(gene)
#参考数据集
cluster_wirka <- read_csv("revision/TCF21_human.csv")
table(cluster_wirka$Cluster)
deg_SMC <- cluster_wirka %>% filter(Cluster=="SMC") %>% 
  pull(gene)
deg_Fibromyocyte <- cluster_wirka %>% filter(Cluster=="Fibromyocyte") %>% 
  pull(gene)
deg_Fibroblast1 <- cluster_wirka %>% filter(Cluster=="Fibroblast") %>% 
  pull(gene)
deg_Fibroblast2 <- cluster_wirka %>% filter(Cluster=="Fibroblast 2") %>% 
  pull(gene)
#取交集
intersect(deg_4,deg_SMC)

intersect(deg_4,deg_Fibromyocyte)
intersect(deg_4,deg_Fibroblast1)
intersect(deg_4,deg_Fibroblast2)
#####2.对VSMC不同亚型进行差异marker分析
Idents(scRNAsub)="celltype"
markers <- FindAllMarkers(scRNAsub, logfc.threshold = 0.25, min.pct = 0.1, 
                          only.pos = TRUE)
markers <- markers[!grepl("MT-", 
                  markers$gene, ignore.case = F),]
markers <- markers[!grepl("^RP[SL]", 
                          markers$gene, ignore.case = F),]
write.csv(markers,file = "revision/celltypeMarker_VSMC.csv")
#####3.singleR
library(SingleR)
library(celldex)
scRNA <- readRDS("SCTtransform/scRNA_SCTransform.rds")
scRNA$seurat_clusters=scRNA$integrated_snn_res.0.4
Idents(object = scRNA) <- "integrated_snn_res.0.4"
blueprint_encode <- BlueprintEncodeData()
hpca <- HumanPrimaryCellAtlasData()
#保存数据
save(blueprint_encode,hpca,file = "revision/single_R_refdata.Rdata")
#进行预测
library(BiocParallel)
#总的
pred.scRNA <- SingleR(test = scRNA@assays$RNA@data, 
                      ref = list(BP=blueprint_encode, HPCA=hpca), 
                      labels = list(blueprint_encode$label.main, 
                                    hpca$label.main) ,
                      clusters = scRNA@active.ident,
                      fine.tune = TRUE, 
                      BPPARAM = MulticoreParam(30))
#查看结果
pred.scRNA$pruned.labels

#查看注释准确性 
pdf(file="revision/singleR_total.pdf",width=10,height=15)   
plotScoreHeatmap(pred.scRNA, clusters=pred.scRNA@rownames, 
                 fontsize.row = 9,show_colnames = T,show.labels=F)
dev.off()
#
pred.scRNA <- SingleR(test = scRNA@assays$RNA@data, 
                      ref = hpca, 
                      labels = hpca$label.main,
                      clusters = scRNA@active.ident,
                      fine.tune = TRUE, 
                      BPPARAM = MulticoreParam(30))
#查看结果
pred.scRNA$pruned.labels

#查看注释准确性 
pdf(file="revision/singleR_HPCA.pdf",width=10,height=10)   
plotScoreHeatmap(pred.scRNA, clusters=pred.scRNA@rownames,
                 order.by = "labels",
                 fontsize.row = 9,show_colnames = T,
                 show.labels=F)
dev.off()
#######二、PA和AC组之间的细胞周期衰老#######
#加载数据
load("~/single-cell/data/Subset_VSMC/scRNAsub_classify.Rdata")

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
p1 <- vln.df%>% dplyr::filter(gene=="IGF1") %>% 
  ggplot(aes(Group,exp))+geom_violin(aes(fill=Group),
                                     scale = "width",
                                     width=0.4,outlier.shape = NA)+
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
  )+geom_text(aes(x=Group,y=3,label=""))+
theme(
  strip.text.y = element_text(size = 8), # 设置分面的字字体大小、颜色、背景、边框，
  axis.text.y = element_text(size = 8),
  axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1,size = 8),
  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
  legend.position = "none"
)
ggsave("revision/03_vln_cell_cycle.pdf",width = 10,height =16,units = "cm")
a <- vln.df %>% dplyr::filter(gene=="IGF1")

wilcox.test(exp~Group,a)
t.test(exp~Group,a)

b <- vln.df %>% dplyr::filter(gene=="CDKN1A")
wilcox.test(exp~Group,b)

t.test(exp~Group,b)


c <- vln.df %>% dplyr::filter(gene=="CDKN2B")
wilcox.test(exp~Group,c)

t.test(exp~Group,c)

#######三、KLF-4的表达差异#######
#加载数据
load("~/single-cell/data/Subset_VSMC/scRNAsub_classify.Rdata")

markers <- c("KLF4")
Idents(scRNAsub)="celltype"
scRNAsub$celltype <- recode(scRNAsub$celltype,
                            AS="AC",
                            normal="PA")

library(reshape2)
vln.df=as.data.frame(scRNAsub@assays$integrated@scale.data[markers,])
colnames(vln.df)[1] <- "KLF4"
vln.df$ID <- rownames(vln.df)

anno=scRNAsub@meta.data[,c("celltype")] %>% as.data.frame() %>% 
  rename(Group=".")

anno$ID <- colnames(scRNAsub)
vln.df=inner_join(vln.df,anno,by="ID")
vln.df$Group <- factor(vln.df$Group,
                       levels = c("Contractile","Chondrocyte-like",
                                  "Fibroblast-like","Macrophage-like"))
#为了控制画图的基因顺序
my_comparisons <- list(c("Contractile","Chondrocyte-like"),
                       c("Contractile","Fibroblast-like"),
                       c("Macrophage-like", "Contractile"),
                       c("Fibroblast-like","Macrophage-like"))
p1 <- vln.df%>% 
  ggplot(aes(Group,KLF4))+geom_violin(aes(fill=Group),
                                      scale = "width",
                                      width=0.4,outlier.shape = NA)+
  scale_fill_npg() +
  scale_x_discrete("")+scale_y_continuous("")+
  theme_minimal()+
  stat_compare_means(comparisons = my_comparisons,
                     method = "t.test",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "***")),
                     label = "p.signif"
  )+geom_text(aes(x=Group,y=3,label=""))+
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1,size = 8),
    legend.position = "none"
  )
p1
ggsave("revision/04_vln_cell_KLF4.pdf",width = 12,height =12,units = "cm")
a <- vln.df %>% filter(Group=="Contractile"|Group=="Chondrocyte-like")

wilcox.test(KLF4~Group,a)
t.test(KLF4~Group,a)

b <- vln.df %>% filter(Group=="Contractile"|Group=="Fibroblast-like")

wilcox.test(KLF4~Group,b)

t.test(KLF4~Group,b)


c <- vln.df %>% filter(Group=="Contractile"|Group=="Macrophage-like")
wilcox.test(KLF4~Group,c)

t.test(KLF4~Group,c)
#####分亚组，根据细胞类型比较AC和PA之间的差异
anno=scRNAsub@meta.data[,c("groups")] %>% as.data.frame() %>% 
  rename(Group_PA_AC=".")
anno$ID <- colnames(scRNAsub)
vln.df=inner_join(vln.df,anno,by="ID")

my_comparisons <- list(c("AC","PA"))
p <- list()
celltype <- names(table(vln.df$Group))
for (i in 1:4) {
  p[[i]] <- vln.df%>% filter(Group==celltype[i]) %>% 
    ggplot(aes(Group_PA_AC,KLF4))+geom_violin(aes(fill=Group_PA_AC),
                                              scale = "width",
                                              width=0.4,outlier.shape = NA)+
    scale_fill_npg() +
    scale_x_discrete("")+scale_y_continuous("")+
    theme_minimal()+
    stat_compare_means(comparisons = my_comparisons,
                       method = "t.test",
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                          symbols = c("****", "***", "**", "*", "NS")),
                       label = "p.signif"
    )+labs(title = paste0(celltype[i],"VSMC")) +
    theme(
      plot.title = element_text(hjust = 0.5,size = 8),
      axis.text.y = element_text(size = 8),
      axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1,size = 8),
      legend.position = "none"
    )
}
plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],labels = LETTERS[1:4])

ggsave("revision/05_vln_cell_KLF4_group.pdf",width = 16,height =16,units = "cm")
a <- vln.df %>% filter(Group=="Contractile")

wilcox.test(KLF4~Group_PA_AC,a)
t.test(KLF4~Group_PA_AC,a)

b <- vln.df %>% filter(Group=="Fibroblast-like")

wilcox.test(KLF4~Group_PA_AC,b)

t.test(KLF4~Group_PA_AC,b)


c <- vln.df %>% filter(Group=="Macrophage-like")
wilcox.test(KLF4~Group_PA_AC,c)

t.test(KLF4~Group_PA_AC,c)

d <- vln.df %>% filter(Group=="Fibroblast-like")

wilcox.test(KLF4~Group_PA_AC,d)

t.test(KLF4~Group_PA_AC,d)
