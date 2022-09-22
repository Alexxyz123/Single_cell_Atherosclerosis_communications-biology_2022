####一、VSMC亚群 cluster 4#####
rm(list = ls())
setwd("D:/Study/single cell/Atherosclerosis/communications Biology/大群分类")
library(tidyverse)
cluster_marker <- read_csv("ClusterMarker.csv")
deg_4 <- cluster_marker %>% filter(cluster==4) %>%  arrange(p_val) %>%
  pull(gene) %>% head(100)
####参考数据集 human
cluster_wirka <- read_csv("TCF21_human.csv")
table(cluster_wirka$Cluster)
intersect_gene_human <- matrix(nrow = 1,ncol = 15)

colnames(intersect_gene_human) <- c(
  "B Cell", "Endothelial","Fibromyocyte","Fibroblast" , "Pericyte 2",   
  "Fibroblast 2", "Macrophage","Mast Cell","Neuron", "NK Cell",  "Pericyte 1",   
  "Plasma Cell 1", "Plasma Cell 2","Tcell","SMC" )
for (j in 1:length(colnames(intersect_gene_human))) {
  deg <- cluster_wirka %>% 
    filter(Cluster==colnames(intersect_gene_human)[j]) %>% 
    pull(gene) %>% toupper()
  #取交集
  deg_4 <- deg_4 %>% toupper() 
  intersect_gene_human[1,j] <- length(intersect(deg_4,deg))
}

#取交集
intersect_gene_human <- intersect_gene_human %>% t()

write.csv(intersect_gene_human,file = "a4.csv")

####参考数据集 mouse_base
dataset <- c("mouse_base","mouse_8w","mouse_16w")
#构架矩阵，输出结果
intersect_gene <- matrix(nrow = 3,ncol = 11)
rownames(intersect_gene) <- dataset
colnames(intersect_gene) <- c(
  "Endothelial 1", "Endothelial 2","Epithelial-like","Fibroblast" ,    
  "Fibroblast 2", "Macrophage","Modulated SMC","Neuron",      
  "SMC1", "SMC2","Tcell" )
for (i in 1:3) {
  file <- paste0(paste0("TCF21_",dataset[i]),".csv")
  cluster_wirka <- read_csv(file)
  for (j in 1:length(colnames(intersect_gene))) {
    deg <- cluster_wirka %>% 
      filter(Cluster==colnames(intersect_gene)[j]) %>% 
      pull(gene) %>% toupper()
    #取交集
    deg_4 <- deg_4 %>% toupper() 
    intersect_gene[i,j] <- length(intersect(deg_4,deg))
  }
}

#柱状图可视化
library(ggplot2)
library(ggpubr)
library(ggsci)
#转换数据
data <- data.frame(celltype=colnames(intersect_gene),
                   mean=1,sd=1,lower_sd=1,
                   higher_sd=1)
sd=c()
for (i in 1:length(colnames(intersect_gene))) {
  data$mean[i] <-round(mean(intersect_gene[,i]),2)
  sd[i] <- round(sd(intersect_gene[,i]),2)
  data$lower_sd[i] <-  round(data$mean[i]-sd[i],2)
  data$higher_sd[i] <- round(data$mean[i]+sd[i],2)
}
data <- intersect_gene %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "celltype") %>% 
  pivot_longer(.,cols =mouse_base:mouse_16w,
               names_to ='time_point',
               values_to = "num_gene") %>% 
  inner_join(data,by="celltype") %>% arrange(mean)
class(data$celltype)
table(data$celltype)
data$celltype <- factor(data$celltype,
                        levels = unique(data$celltype) )
my_comparisons <- list(
  c("Modulated SMC", "Fibroblast")
)
p2 <- ggplot(data = data,aes(x=celltype,y=num_gene)) +  
  geom_bar(data = data, 
           aes(x = celltype, y = mean,
               fill=celltype),
           stat = "identity",
           width = 0.4, 
           position = position_dodge(width = 0.9)) +
  geom_jitter(width = 0.15,size=1)+
  geom_errorbar(data = data,aes(x=celltype,
                                ymin =lower_sd,
                                ymax = higher_sd),
                width = 0.2,
                position = position_dodge(width = 0.9))+
  scale_y_continuous(expand=c(0,0),limits = c(0,34))+
  scale_fill_brewer(palette="Set3")+
  guides(fill=F)+xlab("Cell type")+ylab("Number of genes")+
  theme(plot.title = element_text(size=3))+  #也就加上这一行
  labs(title= "The intersection of top 100 markers in cluster 4 with different clusters in Wirka's data")+
  theme_minimal() +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     method.args = list(alternative = "greater"),
                     label = "p.signif"
                     )+
  coord_flip()
pdf("Wirka_cluster4.pdf",8,6)
p2
dev.off()

                     
####二、VSMC亚群 cluster 3#####
rm(list = ls())
setwd("D:/Study/single cell/Atherosclerosis/communications Biology/大群分类")
library(tidyverse)
cluster_marker <- read_csv("ClusterMarker.csv")
deg_3 <- cluster_marker %>% filter(cluster==3) %>%  arrange(p_val) %>%
  pull(gene) %>% head(100)
####参考数据集 human
cluster_wirka <- read_csv("TCF21_human.csv")
table(cluster_wirka$Cluster)
intersect_gene_human <- matrix(nrow = 1,ncol = 15)

colnames(intersect_gene_human) <- c(
  "B Cell", "Endothelial","Fibromyocyte","Fibroblast" , "Pericyte 2",   
  "Fibroblast 2", "Macrophage","Mast Cell","Neuron", "NK Cell",  "Pericyte 1",   
  "Plasma Cell 1", "Plasma Cell 2","Tcell","SMC" )
for (j in 1:length(colnames(intersect_gene_human))) {
  deg <- cluster_wirka %>% 
    filter(Cluster==colnames(intersect_gene_human)[j]) %>% 
    pull(gene) %>% toupper()
  #取交集
  deg_3 <- deg_3 %>% toupper() 
  intersect_gene_human[1,j] <- length(intersect(deg_3,deg))
}

#取交集
intersect_gene_human <- intersect_gene_human %>% t()

write.csv(intersect_gene_human,file = "a3.csv")

####参考数据集 mouse_base
dataset <- c("mouse_base","mouse_8w","mouse_16w")
#构架矩阵，输出结果
intersect_gene <- matrix(nrow = 3,ncol = 11)
rownames(intersect_gene) <- dataset
colnames(intersect_gene) <- c(
  "Endothelial 1", "Endothelial 2","Epithelial-like","Fibroblast" ,    
  "Fibroblast 2", "Macrophage","Modulated SMC","Neuron",      
  "SMC1", "SMC2","Tcell" )
for (i in 1:3) {
  file <- paste0(paste0("TCF21_",dataset[i]),".csv")
  cluster_wirka <- read_csv(file)
  for (j in 1:length(colnames(intersect_gene))) {
    deg <- cluster_wirka %>% 
      filter(Cluster==colnames(intersect_gene)[j]) %>% 
      pull(gene) %>% toupper()
    #取交集
    deg_3 <- deg_3 %>% toupper() 
    intersect_gene[i,j] <- length(intersect(deg_3,deg))
  }
}

#柱状图可视化
library(ggplot2)
library(ggpubr)
library(ggsci)
#转换数据
data <- data.frame(celltype=colnames(intersect_gene),
                   mean=1,sd=1,lower_sd=1,
                   higher_sd=1)
sd=c()
for (i in 1:length(colnames(intersect_gene))) {
  data$mean[i] <-round(mean(intersect_gene[,i]),2)
  sd[i] <- round(sd(intersect_gene[,i]),2)
  data$lower_sd[i] <-  round(data$mean[i]-sd[i],2)
  data$higher_sd[i] <- round(data$mean[i]+sd[i],2)
}
data <- intersect_gene %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "celltype") %>% 
  pivot_longer(.,cols =mouse_base:mouse_16w,
               names_to ='time_point',
               values_to = "num_gene") %>% 
  inner_join(data,by="celltype") %>% arrange(mean)
class(data$celltype)
table(data$celltype)
data$celltype <- factor(data$celltype,
                        levels = unique(data$celltype) )
#
my_comparisons <- list(
  c("Modulated SMC", "SMC1")
)
p1 <- ggplot(data = data,aes(x=celltype,y=num_gene)) +  
  geom_bar(data = data, 
           aes(x = celltype, y = mean,
               fill=celltype),
           stat = "identity",
           width = 0.4, 
           position = position_dodge(width = 0.9)) +
  geom_jitter(width = 0.15,size=1)+
  geom_errorbar(data = data,aes(x=celltype,
                                ymin =lower_sd,
                                ymax = higher_sd),
                width = 0.2,
                position = position_dodge(width = 0.9))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_brewer(palette="Set3")+
  theme(plot.title = element_text(size = 3))+  #也就加上这一行
  labs(title= "The intersection of top 100 markers in cluster 3 with different clusters in Wirka's data")+
  guides(fill=F)+xlab("Cell type")+ylab("Number of genes")+
  theme_minimal() +
  scale_y_continuous(limits = c(0,42))+
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                         symbols = c("****", "***", "**", "*", "*")),
                     method.args = list(alternative = "greater"),
                     label = "p.signif"
  )+
  coord_flip()
pdf("Wirka_cluster3.pdf",8,6)
p1
dev.off()
####三、VSMC亚群 cluster 15#####
rm(list = ls())
setwd("D:/Study/single cell/Atherosclerosis/communications Biology/大群分类")
library(tidyverse)
cluster_marker <- read_csv("ClusterMarker.csv")
deg15 <- cluster_marker %>% filter(cluster==15) %>%  arrange(p_val) %>%
  pull(gene) %>% head(100)
####参考数据集 human
cluster_wirka <- read_csv("TCF21_human.csv")
table(cluster_wirka$Cluster)
intersect_gene_human <- matrix(nrow = 1,ncol = 15)

colnames(intersect_gene_human) <- c(
  "B Cell", "Endothelial","Fibromyocyte","Fibroblast" , "Pericyte 2",   
  "Fibroblast 2", "Macrophage","Mast Cell","Neuron", "NK Cell",  "Pericyte 1",   
  "Plasma Cell 1", "Plasma Cell 2","Tcell","SMC" )
for (j in 1:length(colnames(intersect_gene_human))) {
  deg <- cluster_wirka %>% 
    filter(Cluster==colnames(intersect_gene_human)[j]) %>% 
    pull(gene) %>% toupper()
  #取交集
  deg15 <- deg15 %>% toupper() 
  intersect_gene_human[1,j] <- length(intersect(deg15,deg))
}

#取交集
intersect_gene_human <- intersect_gene_human %>% t()

write.csv(intersect_gene_human,file = "a15.csv")
####参考数据集 mouse_base
dataset <- c("mouse_base","mouse_8w","mouse_16w")
#构架矩阵，输出结果
intersect_gene <- matrix(nrow = 3,ncol = 11)
rownames(intersect_gene) <- dataset
colnames(intersect_gene) <- c(
  "Endothelial 1", "Endothelial 2","Epithelial-like","Fibroblast" ,    
  "Fibroblast 2", "Macrophage","Modulated SMC","Neuron",      
  "SMC1", "SMC2","Tcell" )
for (i in 1:3) {
  file <- paste0(paste0("TCF21_",dataset[i]),".csv")
  cluster_wirka <- read_csv(file)
  for (j in 1:length(colnames(intersect_gene))) {
    deg <- cluster_wirka %>% 
      filter(Cluster==colnames(intersect_gene)[j]) %>% 
      pull(gene) %>% toupper()
    #取交集
    deg15 <- deg15 %>% toupper() 
    intersect_gene[i,j] <- length(intersect(deg15,deg))
  }
}

#柱状图可视化
library(ggplot2)
library(ggpubr)
library(ggsci)
#转换数据
data <- data.frame(celltype=colnames(intersect_gene),
                   mean=1,sd=1,lower_sd=1,
                   higher_sd=1)
sd=c()
for (i in 1:length(colnames(intersect_gene))) {
  data$mean[i] <-round(mean(intersect_gene[,i]),2)
  sd[i] <- round(sd(intersect_gene[,i]),2)
  data$lower_sd[i] <-  round(data$mean[i]-sd[i],2)
  data$higher_sd[i] <- round(data$mean[i]+sd[i],2)
}
data <- intersect_gene %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "celltype") %>% 
  pivot_longer(.,cols =mouse_base:mouse_16w,
               names_to ='time_point',
               values_to = "num_gene") %>% 
  inner_join(data,by="celltype") %>% arrange(mean)
class(data$celltype)
table(data$celltype)
data$celltype <- factor(data$celltype,
                        levels = unique(data$celltype) )
#
my_comparisons <- list(
  c("Modulated SMC", "SMC1")
)
p3 <- ggplot(data = data,aes(x=celltype,y=num_gene)) +  
  geom_bar(data = data, 
           aes(x = celltype, y = mean,
               fill=celltype),
           stat = "identity",
           width = 0.4, 
           position = position_dodge(width = 0.9)) +
  geom_jitter(width = 0.15,size=1)+
  geom_errorbar(data = data,aes(x=celltype,
                                ymin =lower_sd,
                                ymax = higher_sd),
                width = 0.2,
                position = position_dodge(width = 0.9))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_brewer(palette="Set3")+
  theme(plot.title = element_text(size = 3))+  #也就加上这一行
  labs(title= "The intersection of top 100 markers in cluster 15 with different clusters in Wirka's data")+
  guides(fill=F)+xlab("Cell type")+ylab("Number of genes")+
  theme_minimal()+
  scale_y_continuous(limits = c(0,17))+
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     method.args = list(alternative = "greater"),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "*")),
                     label = "p.signif") +
  coord_flip()
pdf("Wirka_cluster15.pdf",8,6)
p1
dev.off()
####四、巨噬细胞亚群 cluster 5#####
rm(list = ls())
library(tidyverse)
cluster_marker <- read_csv("ClusterMarker.csv")
deg_5 <- cluster_marker %>% filter(cluster==5) %>%  arrange(p_val) %>%
  pull(gene) %>% head(100)
####参考数据集 human
cluster_wirka <- read_csv("TCF21_human.csv")
table(cluster_wirka$Cluster)
intersect_gene_human <- matrix(nrow = 1,ncol = 15)

colnames(intersect_gene_human) <- c(
  "B Cell", "Endothelial","Fibromyocyte","Fibroblast" , "Pericyte 2",   
  "Fibroblast 2", "Macrophage","Mast Cell","Neuron", "NK Cell",  "Pericyte 1",   
  "Plasma Cell 1", "Plasma Cell 2","Tcell","SMC" )
for (j in 1:length(colnames(intersect_gene_human))) {
  deg <- cluster_wirka %>% 
    filter(Cluster==colnames(intersect_gene_human)[j]) %>% 
    pull(gene) %>% toupper()
  #取交集
  deg_5 <- deg_5 %>% toupper() 
  intersect_gene_human[1,j] <- length(intersect(deg_5,deg))
}


#取交集
intersect_gene_human <- intersect_gene_human %>% t()

write.csv(intersect_gene_human,file = "a5.csv")
####参考数据集 mouse_base#
dataset <- c("mouse_base","mouse_8w","mouse_16w")
#构架矩阵，输出结果
intersect_gene <- matrix(nrow = 3,ncol = 11)
rownames(intersect_gene) <- dataset
colnames(intersect_gene) <- c(
  "Endothelial 1", "Endothelial 2","Epithelial-like","Fibroblast" ,    
  "Fibroblast 2", "Macrophage","Modulated SMC","Neuron",      
  "SMC1", "SMC2","Tcell" )
for (i in 1:3) {
  file <- paste0(paste0("TCF21_",dataset[i]),".csv")
  cluster_wirka <- read_csv(file)
  for (j in 1:length(colnames(intersect_gene))) {
    deg <- cluster_wirka %>% 
      filter(Cluster==colnames(intersect_gene)[j]) %>% 
      pull(gene) %>% toupper()
    #取交集
    deg_5 <- deg_5 %>% toupper() 
    intersect_gene[i,j] <- length(intersect(deg_5,deg))
  }
}

#柱状图可视化
library(ggplot2)
library(ggpubr)
library(ggsci)
#转换数据
data <- data.frame(celltype=colnames(intersect_gene),
                   mean=1,sd=1,lower_sd=1,
                   higher_sd=1)
sd=c()
for (i in 1:length(colnames(intersect_gene))) {
  data$mean[i] <-round(mean(intersect_gene[,i]),2)
  sd[i] <- round(sd(intersect_gene[,i]),2)
  data$lower_sd[i] <-  round(data$mean[i]-sd[i],2)
  data$higher_sd[i] <- round(data$mean[i]+sd[i],2)
}
data <- intersect_gene %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "celltype") %>% 
  pivot_longer(.,cols =mouse_base:mouse_16w,
               names_to ='time_point',
               values_to = "num_gene") %>% 
  inner_join(data,by="celltype") %>% arrange(mean)
class(data$celltype)
table(data$celltype)
data$celltype <- factor(data$celltype,
                        levels = unique(data$celltype) )
my_comparisons <- list(
  c("Macrophage", "Modulated SMC")
)
p4 <- ggplot(data = data,aes(x=celltype,y=num_gene)) +  
  geom_bar(data = data, 
           aes(x = celltype, y = mean,
               fill=celltype),
           stat = "identity",
           width = 0.4, 
           position = position_dodge(width = 0.9)) +
  geom_jitter(width = 0.15,size=1)+
  geom_errorbar(data = data,aes(x=celltype,
                                ymin =lower_sd,
                                ymax = higher_sd),
                width = 0.2,
                position = position_dodge(width = 0.9))+
  scale_y_continuous(expand=c(0,0),limits = c(0,40))+
  scale_fill_brewer(palette="Set3")+
  guides(fill=F)+xlab("Cell type")+ylab("Number of genes")+
  theme(plot.title = element_text(size=3))+  #也就加上这一行
  labs(title= "The intersection of top 100 markers in cluster 5 with different clusters in Wirka's data")+
  theme_minimal() +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     method.args = list(alternative = "greater"),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "*")),
                     label = "p.signif"
  )+
  coord_flip()
pdf("Wirka_cluster5.pdf",8,6)
p1
dev.off()

####五、巨噬细胞亚群 cluster 7#####
rm(list = ls())
library(tidyverse)
cluster_marker <- read_csv("ClusterMarker.csv")
deg_7 <- cluster_marker %>% filter(cluster==7) %>%  arrange(p_val) %>%
  pull(gene) %>% head(100)
####参考数据集 human
cluster_wirka <- read_csv("TCF21_human.csv")
table(cluster_wirka$Cluster)
intersect_gene_human <- matrix(nrow = 1,ncol = 15)

colnames(intersect_gene_human) <- c(
  "B Cell", "Endothelial","Fibromyocyte","Fibroblast" , "Pericyte 2",   
  "Fibroblast 2", "Macrophage","Mast Cell","Neuron", "NK Cell",  "Pericyte 1",   
  "Plasma Cell 1", "Plasma Cell 2","Tcell","SMC" )
for (j in 1:length(colnames(intersect_gene_human))) {
  deg <- cluster_wirka %>% 
    filter(Cluster==colnames(intersect_gene_human)[j]) %>% 
    pull(gene) %>% toupper()
  #取交集
  deg_7 <- deg_7 %>% toupper() 
  intersect_gene_human[1,j] <- length(intersect(deg_7,deg))
}

#取交集
intersect_gene_human <- intersect_gene_human %>% t()

write.csv(intersect_gene_human,file = "a7.csv")

####参考数据集 mouse_base#
dataset <- c("mouse_base","mouse_8w","mouse_16w")
#构架矩阵，输出结果
intersect_gene <- matrix(nrow = 3,ncol = 11)
rownames(intersect_gene) <- dataset
colnames(intersect_gene) <- c(
  "Endothelial 1", "Endothelial 2","Epithelial-like","Fibroblast" ,    
  "Fibroblast 2", "Macrophage","Modulated SMC","Neuron",      
  "SMC1", "SMC2","Tcell" )
for (i in 1:3) {
  file <- paste0(paste0("TCF21_",dataset[i]),".csv")
  cluster_wirka <- read_csv(file)
  for (j in 1:length(colnames(intersect_gene))) {
    deg <- cluster_wirka %>% 
      filter(Cluster==colnames(intersect_gene)[j]) %>% 
      pull(gene) %>% toupper()
    #取交集
    deg_7 <- deg_7 %>% toupper() 
    intersect_gene[i,j] <- length(intersect(deg_7,deg))
  }
}

#柱状图可视化
library(ggplot2)
library(ggpubr)
library(ggsci)
#转换数据
data <- data.frame(celltype=colnames(intersect_gene),
                   mean=1,sd=1,lower_sd=1,
                   higher_sd=1)
sd=c()
for (i in 1:length(colnames(intersect_gene))) {
  data$mean[i] <-round(mean(intersect_gene[,i]),2)
  sd[i] <- round(sd(intersect_gene[,i]),2)
  data$lower_sd[i] <-  round(data$mean[i]-sd[i],2)
  data$higher_sd[i] <- round(data$mean[i]+sd[i],2)
}
data <- intersect_gene %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "celltype") %>% 
  pivot_longer(.,cols =mouse_base:mouse_16w,
               names_to ='time_point',
               values_to = "num_gene") %>% 
  inner_join(data,by="celltype") %>% arrange(mean)
class(data$celltype)
table(data$celltype)
data$celltype <- factor(data$celltype,
                        levels = unique(data$celltype) )
my_comparisons <- list(
  c("Macrophage", "Modulated SMC")
)
p5 <- ggplot(data = data,aes(x=celltype,y=num_gene)) +  
  geom_bar(data = data, 
           aes(x = celltype, y = mean,
               fill=celltype),
           stat = "identity",
           width = 0.4, 
           position = position_dodge(width = 0.9)) +
  geom_jitter(width = 0.15,size=1)+
  geom_errorbar(data = data,aes(x=celltype,
                                ymin =lower_sd,
                                ymax = higher_sd),
                width = 0.2,
                position = position_dodge(width = 0.9))+
  scale_y_continuous(expand=c(0,0),limits = c(0,21))+
  scale_fill_brewer(palette="Set3")+
  guides(fill=F)+xlab("Cell type")+ylab("Number of genes")+
  theme(plot.title = element_text(size=3))+  #也就加上这一行
  labs(title= "The intersection of top 100 markers in cluster 7 with different clusters in Wirka's data")+
  theme_minimal() +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     method.args = list(alternative = "greater"),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "*")),
                     label = "p.signif"
  )+
  coord_flip()
pdf("Wirka_cluster7.pdf",8,6)
p1
dev.off()
####六、巨噬细胞亚群 cluster 8#####
rm(list = ls())
library(tidyverse)
cluster_marker <- read_csv("ClusterMarker.csv")
deg_8 <- cluster_marker %>% filter(cluster==8) %>%  arrange(p_val) %>%
  pull(gene) %>% head(100)
####参考数据集 human
cluster_wirka <- read_csv("TCF21_human.csv")
table(cluster_wirka$Cluster)
intersect_gene_human <- matrix(nrow = 1,ncol = 15)

colnames(intersect_gene_human) <- c(
  "B Cell", "Endothelial","Fibromyocyte","Fibroblast" , "Pericyte 2",   
  "Fibroblast 2", "Macrophage","Mast Cell","Neuron", "NK Cell",  "Pericyte 1",   
  "Plasma Cell 1", "Plasma Cell 2","Tcell","SMC" )
for (j in 1:length(colnames(intersect_gene_human))) {
  deg <- cluster_wirka %>% 
    filter(Cluster==colnames(intersect_gene_human)[j]) %>% 
    pull(gene) %>% toupper()
  #取交集
  deg_8 <- deg_8 %>% toupper() 
  intersect_gene_human[1,j] <- length(intersect(deg_8,deg))
}

#取交集
intersect_gene_human <- intersect_gene_human %>% t()

write.csv(intersect_gene_human,file = "a8.csv")
####参考数据集 mouse_base#
dataset <- c("mouse_base","mouse_8w","mouse_16w")
#构架矩阵，输出结果
intersect_gene <- matrix(nrow = 3,ncol = 11)
rownames(intersect_gene) <- dataset
colnames(intersect_gene) <- c(
  "Endothelial 1", "Endothelial 2","Epithelial-like","Fibroblast" ,    
  "Fibroblast 2", "Macrophage","Modulated SMC","Neuron",      
  "SMC1", "SMC2","Tcell" )
for (i in 1:3) {
  file <- paste0(paste0("TCF21_",dataset[i]),".csv")
  cluster_wirka <- read_csv(file)
  for (j in 1:length(colnames(intersect_gene))) {
    deg <- cluster_wirka %>% 
      filter(Cluster==colnames(intersect_gene)[j]) %>% 
      pull(gene) %>% toupper()
    #取交集
    deg_8 <- deg_8 %>% toupper() 
    intersect_gene[i,j] <- length(intersect(deg_8,deg))
  }
}

#柱状图可视化
library(ggplot2)
library(ggpubr)
library(ggsci)
#转换数据
data <- data.frame(celltype=colnames(intersect_gene),
                   mean=1,sd=1,lower_sd=1,
                   higher_sd=1)
sd=c()
for (i in 1:length(colnames(intersect_gene))) {
  data$mean[i] <-round(mean(intersect_gene[,i]),2)
  sd[i] <- round(sd(intersect_gene[,i]),2)
  data$lower_sd[i] <-  round(data$mean[i]-sd[i],2)
  data$higher_sd[i] <- round(data$mean[i]+sd[i],2)
}
data <- intersect_gene %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "celltype") %>% 
  pivot_longer(.,cols =mouse_base:mouse_16w,
               names_to ='time_point',
               values_to = "num_gene") %>% 
  inner_join(data,by="celltype") %>% arrange(mean)
class(data$celltype)
table(data$celltype)
data$celltype <- factor(data$celltype,
                        levels = unique(data$celltype) )
my_comparisons <- list(
  c("Macrophage", "Tcell")
)
p6 <- ggplot(data = data,aes(x=celltype,y=num_gene)) +  
  geom_bar(data = data, 
           aes(x = celltype, y = mean,
               fill=celltype),
           stat = "identity",
           width = 0.4, 
           position = position_dodge(width = 0.9)) +
  geom_jitter(width = 0.15,size=1)+
  geom_errorbar(data = data,aes(x=celltype,
                                ymin =lower_sd,
                                ymax = higher_sd),
                width = 0.2,
                position = position_dodge(width = 0.9))+
  scale_y_continuous(expand=c(0,0),limits = c(0,25))+
  scale_fill_brewer(palette="Set3")+
  guides(fill=F)+xlab("Cell type")+ylab("Number of genes")+
  theme(plot.title = element_text(size=3))+  #也就加上这一行
  labs(title= "The intersection of top 100 markers in cluster 8 with different clusters in Wirka's data")+
  theme_minimal() +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     method.args = list(alternative = "greater"),
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "*")),
                     label = "p.signif"
  )+
  coord_flip()
pdf("Wirka_cluster8.pdf",8,6)
p1
dev.off()
#####绘图###
pdf("01_mouse_intersect.pdf",16,12)
p1 + p2 + p3 + p4 +p5+p6+
  plot_layout(guides = 'collect',ncol = 2)
dev.off()


#####Microanatomy_Marie######
rm(list = ls())
setwd("D:/Study/single cell/Atherosclerosis/communications Biology/大群分类")
library(tidyverse)
cluster_marker <- read_csv("ClusterMarker.csv")
#构架矩阵，输出结果
intersect_gene <- matrix(nrow = 1,ncol = 19)
rownames(intersect_gene) <- "num_gene"
colnames(intersect_gene) <- paste0("Cluster ",0:18)
Marie_marker <- read_csv("Microanatomy_Marie.csv")
deg_Marie <- Marie_marker %>%  arrange(p_val) %>%
  pull(gene) %>% head(100)
for (i in 0:18) {
  deg <- cluster_marker %>% filter(cluster==i) %>% 
    arrange(p_val) %>%
    pull(gene) %>% head(100)
  #取交集
  intersect_gene[1,i+1] <- length(intersect(deg,deg_Marie))
}

data <- intersect_gene %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "cluster") %>% arrange(num_gene)
#柱状图可视化
library(ggplot2)
library(ggpubr)
library(ggsci)
#转换数据
data$cluster <- factor(data$cluster,
                       levels = unique(data$cluster) )
cols <- c("Cluster 3"="#E64B35FF",
          "Cluster 4"="#4DBBD5FF",
          "Cluster 15"="#00A087FF",
          "Cluster 12"="#3C5488FF",
          "Cluster 2"="#F39B7FFF",
          "Cluster 7"="#8491B4FF",
          "Cluster 8"="#91D1C2FF")
p1 <- ggplot() +  
  geom_bar(data = data, 
           aes(x = cluster, y = num_gene,
               fill=cluster),
           stat = "identity",
           width = 0.4, 
           position = position_dodge(width = 0.9)) +
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values =cols)+
  guides(fill=F)+
  theme(plot.title = element_text(hjust = 0.5))+  #也就加上这一行
  labs(title= "The intersection of top 100 markers in different clusters with SMC subset in Depuydt's data",
       y="Number of genes",x="")+
  theme_minimal() +
  coord_flip()

pdf(file = "Depuydt_SMC.pdf",6,10)
p1
dev.off()

