
####一、Macrophage重分组 cluster 0 Trem2#####
rm(list = ls())
library(tidyverse)
cluster_marker <- readxl::read_xlsx("Macrophage_cluster_marker.xlsx")
deg_0 <- cluster_marker %>% filter(cluster==0) %>%  arrange(p_val) %>%
  pull(gene) %>% head(100)
####参考数据集 human
cluster_zernecke <- readxl::read_xlsx("Zernecke_macrophage.xlsx")
table(cluster_zernecke$cluster)
intersect_gene_human <- matrix(nrow = 1,ncol = 16)

colnames(intersect_gene_human) <- names(table(cluster_zernecke$cluster))[-12]
for (j in 1:length(colnames(intersect_gene_human))) {
  deg <- cluster_zernecke %>% 
    filter(cluster==colnames(intersect_gene_human)[j]) %>% 
    pull(gene...8) %>% toupper()
  #取交集
  deg_0 <- deg_0 %>% toupper() 
  intersect_gene_human[1,j] <- length(intersect(deg_0,deg))
}

#取交集
intersect_gene_human <- intersect_gene_human %>% t() %>% as.data.frame()
colnames(intersect_gene_human) <- "num_gene"
intersect_gene_human$celltype <- rownames(intersect_gene_human)
#write.csv(intersect_gene_human,file = "Macrophage_a0.csv")
data <- intersect_gene_human
data <- data %>% arrange(num_gene)
#柱状图可视化
library(ggplot2)
library(ggpubr)
library(ggsci)
data$celltype <- factor(data$celltype,
                        levels = unique(data$celltype) )
#转换数据
p1 <- ggplot(data = data,aes(x=celltype,y=num_gene)) +  
  geom_bar(data = data, 
           aes(x = celltype, y = num_gene,
               fill=celltype),
           stat = "identity",
           width = 0.4, 
           position = position_dodge(width = 0.9)) +
  scale_y_continuous(expand=c(0,0),limits = c(0,34))+
  #scale_fill_brewer(palette="Set3")+
  guides(fill=F)+xlab("Cell type")+ylab("Number of genes")+
  theme(plot.title = element_text(size=0.5,hjust = 0.5))+  #也就加上这一行
  labs(title= "TREM2-hi macrophage")+
  theme_minimal() +
  coord_flip()
pdf("06_zernecke_cluster0.pdf",8,6)
p1
dev.off()

####二、Macrophage重分组 cluster 1 resident-like#####

library(tidyverse)
cluster_marker <- readxl::read_xlsx("Macrophage_cluster_marker.xlsx")
deg_1 <- cluster_marker %>% filter(cluster==1) %>%  arrange(p_val) %>%
  pull(gene) %>% head(100)
####参考数据集 human
cluster_zernecke <- readxl::read_xlsx("Zernecke_macrophage.xlsx")
table(cluster_zernecke$cluster)
intersect_gene_human <- matrix(nrow = 1,ncol = 16)

colnames(intersect_gene_human) <- names(table(cluster_zernecke$cluster))[-12]
for (j in 1:length(colnames(intersect_gene_human))) {
  deg <- cluster_zernecke %>% 
    filter(cluster==colnames(intersect_gene_human)[j]) %>% 
    pull(gene...8) %>% toupper()
  #取交集
  deg_1 <- deg_1 %>% toupper() 
  intersect_gene_human[1,j] <- length(intersect(deg_1,deg))
}

#取交集
intersect_gene_human <- intersect_gene_human %>% t() %>% as.data.frame()
colnames(intersect_gene_human) <- "num_gene"
intersect_gene_human$celltype <- rownames(intersect_gene_human)
#write.csv(intersect_gene_human,file = "Macrophage_a1.csv")
data <- intersect_gene_human
data <- data %>% arrange(num_gene)
#柱状图可视化
library(ggplot2)
library(ggpubr)
library(ggsci)
data$celltype <- factor(data$celltype,
                        levels = unique(data$celltype) )
#转换数据
p2 <- ggplot(data = data,aes(x=celltype,y=num_gene)) +  
  geom_bar(data = data, 
           aes(x = celltype, y = num_gene,
               fill=celltype),
           stat = "identity",
           width = 0.4, 
           position = position_dodge(width = 0.9)) +
  scale_y_continuous(expand=c(0,0),limits = c(0,34))+
  #scale_fill_brewer(palette="Set3")+
  guides(fill=F)+xlab("Cell type")+ylab("Number of genes")+
  theme(plot.title = element_text(size=0.5,hjust = 0.5))+  #也就加上这一行
  labs(title= "Resident-like macrophage")+
  theme_minimal() +
  coord_flip()
pdf("06_zernecke_cluster1.pdf",8,6)
p2
dev.off()

####三、Macrophage重分组 cluster 2 Inflammatory like#####

library(tidyverse)
cluster_marker <- readxl::read_xlsx("Macrophage_cluster_marker.xlsx")
deg_2 <- cluster_marker %>% filter(cluster==2) %>%  arrange(p_val) %>%
  pull(gene) %>% head(100)
####参考数据集 human
cluster_zernecke <- readxl::read_xlsx("Zernecke_macrophage.xlsx")
table(cluster_zernecke$cluster)
intersect_gene_human <- matrix(nrow = 1,ncol = 16)

colnames(intersect_gene_human) <- names(table(cluster_zernecke$cluster))[-12]
for (j in 1:length(colnames(intersect_gene_human))) {
  deg <- cluster_zernecke %>% 
    filter(cluster==colnames(intersect_gene_human)[j]) %>% 
    pull(gene...8) %>% toupper()
  #取交集
  deg_2 <- deg_2 %>% toupper() 
  intersect_gene_human[1,j] <- length(intersect(deg_2,deg))
}

#取交集
intersect_gene_human <- intersect_gene_human %>% t() %>% as.data.frame()
colnames(intersect_gene_human) <- "num_gene"
intersect_gene_human$celltype <- rownames(intersect_gene_human)
#write.csv(intersect_gene_human,file = "Macrophage_a2.csv")
data <- intersect_gene_human
data <- data %>% arrange(num_gene)
#柱状图可视化
library(ggplot2)
library(ggpubr)
library(ggsci)
data$celltype <- factor(data$celltype,
                        levels = unique(data$celltype) )
#转换数据
p3 <- ggplot(data = data,aes(x=celltype,y=num_gene)) +  
  geom_bar(data = data, 
           aes(x = celltype, y = num_gene,
               fill=celltype),
           stat = "identity",
           width = 0.4, 
           position = position_dodge(width = 0.9)) +
  scale_y_continuous(expand=c(0,0),limits = c(0,34))+
  #scale_fill_brewer(palette="Set3")+
  guides(fill=F)+xlab("Cell type")+ylab("Number of genes")+
  theme(plot.title = element_text(size=0.5,hjust = 0.5))+  #也就加上这一行
  labs(title= "Inflamamtory macrophage")+
  theme_minimal() +
  coord_flip()
pdf("06_zernecke_cluster2.pdf",8,6)
p3
dev.off()

library(cowplot)
pdf("06_zernecke_macrophage.pdf",10,6)
plot_grid(p1,p2,p3,ncol = 2,labels = LETTERS[1:3])
dev.off()



