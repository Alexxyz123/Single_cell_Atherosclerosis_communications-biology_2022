setwd("./CellType")
library(Seurat)
library(RColorBrewer)
###加载颜色#####
color=brewer.pal(9,"Reds")
color
pie(rep(1,9),col = color,labels = color)
#####SMC cell marker####
#save(scRNA,file = "scRNA_cellmarker.RData")
scRNA <- get(load("scRNA_celltype.Rdata"))
####开始画图
SMC <- FeaturePlot(scRNA,
            features = c("TAGLN","ACTA2","MYH11"
                               ,"PDGFRB"),reduction = "umap",
            cols = c("grey","#EF3B2C","firebrick3"))+ guides(color=F)
ggsave(SMC,filename = "SMC_Feature.png",width = 20, height = 20)
#####EC cell marker####
EC <- FeaturePlot(scRNA,
            features = c("CD34","VWF","CLEC14A","PECAM1"),reduction = "umap",
            cols =  c("grey","#EF3B2C","firebrick3"))
ggsave(EC,filename = "EC_Feature.png",width = 20, height = 20)
#####T cell marker#####
T_cell <- FeaturePlot(scRNA,
            features = c("CD3E","CD4","CD8A","CD8B"),reduction = "umap",
            cols = c("grey","#EF3B2C","firebrick3") )
ggsave(T_cell,filename = "T_cell_Feature.png",width = 20, height = 20)
#####B cell marker#####
B_cell <- FeaturePlot(scRNA,
            features = c("CD79A","CD19"),reduction = "umap",
            cols = c("grey","#EF3B2C","firebrick3"))
ggsave(B_cell,filename = "B_cell_Feature.png",width = 20, height = 10)
######monocyte/Macrophage cell marker####
MC <- FeaturePlot(scRNA,
            features = c("CD68","CD14","CX3CR1","LYZ"),reduction = "umap",
            cols = c("grey","#EF3B2C","firebrick3"))
ggsave(MC,filename = "MyeloidCells_Feature.png",width = 20, height = 20)
####NK cell marker#####
NK_cell <- FeaturePlot(scRNA,
            features = c("KLRB1","GZMB","NKG7","CD160"),reduction = "umap",
            cols =c("grey","#EF3B2C","firebrick3"))
ggsave(NK_cell,filename = "NK_cell_Feature.png",width = 20, height = 20)
####DC cell marker######
DC_cell <- FeaturePlot(scRNA,
            features = c("CD1C"),reduction = "umap",
            cols = c("grey","#EF3B2C","firebrick3"))
ggsave(DC_cell,filename = "DC_cell_Feature.png",width = 10, height = 10)

######KIT+ stem cells######
STEM <- FeaturePlot(scRNA,
            features = c("KIT"),reduction = "umap",
            cols = c("grey","#EF3B2C","firebrick3"))
ggsave(STEM,filename = "stem_cell_Feature.png",width = 10, height = 10)
######pDC cell markers#####
pDC <- FeaturePlot(scRNA,
            features = c("IRF8","TCF4","LILRA4","IL3RA","NRP1"),reduction = "umap",
            cols = c("grey","#EF3B2C","firebrick3"))
ggsave(pDC,filename = "pDC_cell_Feature.png",width = 20, height = 20)

DefaultAssay(scRNA) <- "RNA"
DefaultAssay(scRNA) <- "integrated"


pdf("cell_feature.pdf",12,6)
FeaturePlot(scRNA,features = c("TAGLN","ACTA2","MYH11",
                               "PDGFRB","CD68","CD14","CX3CR1","LYZ"),
            ncol = 4,
            reduction = "umap",cols = c("grey","#EF3B2C","firebrick3"))

dev.off()




