
####一、VSMC重分组 cluster 0#####
rm(list = ls())
library(tidyverse)
cluster_marker <- readxl::read_xlsx("VSMC_cluster_marker.xlsx")
deg_0 <- cluster_marker %>% filter(cluster==0) %>%  arrange(p_val) %>%
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
  deg_0 <- deg_0 %>% toupper() 
  intersect_gene_human[1,j] <- length(intersect(deg_0,deg))
}

#取交集
intersect_gene_human <- intersect_gene_human %>% t()

write.csv(intersect_gene_human,file = "VSMC_a0.csv")

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
    deg_0 <- deg_0 %>% toupper() 
    intersect_gene[i,j] <- length(intersect(deg_0,deg))
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
  c("SMC1", "Fibroblast 2")
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
  scale_y_continuous(expand=c(0,0),limits = c(0,34))+
  scale_fill_brewer(palette="Set3")+
  guides(fill=F)+xlab("Cell type")+ylab("Number of genes")+
  theme(plot.title = element_text(size=3))+  #也就加上这一行
  labs(title= "The intersection of top 100 markers in cluster 0 with different clusters in Wirka's data")+
  theme_minimal() +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     method.args = list(alternative = "greater"),
                     label = "p.signif"
  )+
  coord_flip()
pdf("Wirka_cluster0.pdf",8,6)
p1
dev.off()

####二、VSMC重分组 cluster 1#####
rm(list = ls())
library(tidyverse)
cluster_marker <- readxl::read_xlsx("VSMC_cluster_marker.xlsx")
deg_1 <- cluster_marker %>% filter(cluster==1) %>%  arrange(p_val) %>%
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
  deg_1 <- deg_1 %>% toupper() 
  intersect_gene_human[1,j] <- length(intersect(deg_1,deg))
}

#取交集
intersect_gene_human <- intersect_gene_human %>% t()

write.csv(intersect_gene_human,file = "VSMC_a1.csv")

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
    deg_1 <- deg_1 %>% toupper() 
    intersect_gene[i,j] <- length(intersect(deg_1,deg))
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
  c("SMC1", "Fibroblast 2")
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
  scale_y_continuous(expand=c(0,0),limits = c(0,34))+
  scale_fill_brewer(palette="Set3")+
  guides(fill=F)+xlab("Cell type")+ylab("Number of genes")+
  theme(plot.title = element_text(size=3))+  #也就加上这一行
  labs(title= "The intersection of top 100 markers in cluster 1 with different clusters in Wirka's data")+
  theme_minimal() +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     method.args = list(alternative = "greater"),
                     label = "p.signif"
  )+
  coord_flip()
pdf("Wirka_cluster0.pdf",8,6)
p1
dev.off()

####一、VSMC重分组 cluster 2#####
rm(list = ls())
library(tidyverse)
cluster_marker <- readxl::read_xlsx("VSMC_cluster_marker.xlsx")
deg_2 <- cluster_marker %>% filter(cluster==2) %>%  arrange(p_val) %>%
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
  deg_2 <- deg_2 %>% toupper() 
  intersect_gene_human[1,j] <- length(intersect(deg_2,deg))
}

#取交集
intersect_gene_human <- intersect_gene_human %>% t()

write.csv(intersect_gene_human,file = "VSMC_a2.csv")

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
    deg_2 <- deg_2 %>% toupper() 
    intersect_gene[i,j] <- length(intersect(deg_2,deg))
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
  c("SMC1", "Fibroblast 2")
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
  scale_y_continuous(expand=c(0,0),limits = c(0,34))+
  scale_fill_brewer(palette="Set3")+
  guides(fill=F)+xlab("Cell type")+ylab("Number of genes")+
  theme(plot.title = element_text(size=3))+  #也就加上这一行
  labs(title= "The intersection of top 100 markers in cluster 2 with different clusters in Wirka's data")+
  theme_minimal() +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     method.args = list(alternative = "greater"),
                     label = "p.signif"
  )+
  coord_flip()
pdf("Wirka_cluster0.pdf",8,6)
p1
dev.off()

####一、VSMC重分组 cluster 3#####
rm(list = ls())
library(tidyverse)
cluster_marker <- readxl::read_xlsx("VSMC_cluster_marker.xlsx")
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

write.csv(intersect_gene_human,file = "VSMC_a3.csv")

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
my_comparisons <- list(
  c("SMC1", "Fibroblast 2")
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
  scale_y_continuous(expand=c(0,0),limits = c(0,34))+
  scale_fill_brewer(palette="Set3")+
  guides(fill=F)+xlab("Cell type")+ylab("Number of genes")+
  theme(plot.title = element_text(size=3))+  #也就加上这一行
  labs(title= "The intersection of top 100 markers in cluster 3 with different clusters in Wirka's data")+
  theme_minimal() +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     method.args = list(alternative = "greater"),
                     label = "p.signif"
  )+
  coord_flip()
pdf("Wirka_cluster0.pdf",8,6)
p1
dev.off()



