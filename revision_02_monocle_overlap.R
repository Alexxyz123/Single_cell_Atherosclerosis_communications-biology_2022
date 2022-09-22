library(monocle)
library(Seurat)
library(tidyverse)
####1.获取state的差异marker基因#####
#加载数据
load("scRNAsub_classify.Rdata")
#VSMC monocle数据
load("mycds.Rdata")
###获得表型信息
p_data <- subset(phenoData(mycds)@data,select='State')
scRNAsub <- AddMetaData(scRNAsub,p_data,col.name = 'State')
####查看PA组中的macrophage-like细胞在state4中一共有多少个
Idents(scRNAsub) <- "State"
ClusterMarker <- FindAllMarkers(scRNAsub, assay = "integrated", 
                                slot = "data", only.pos = T,
                                logfc.threshold = 0.1, min.pct = 0.25)
#保存结果
write.csv(ClusterMarker,file = "Clustermarker_state.csv")
####2.取交集####
ClusterMarker <- read_csv("Clustermarker_state.csv")
#构架矩阵，输出结果
intersect_gene <- matrix(nrow = 3,ncol = 12)
dataset <- c("mouse_base","mouse_8w","mouse_16w")
rownames(intersect_gene) <- dataset
colnames(intersect_gene) <- names(table(cluster_wirka$Cluster))
list <- list()
list_p <- list()
p <- intersect_gene
for (j in 1:5) {
  for (i in 1:3) {
    file <- paste0(paste0("TCF21_",dataset[i]),".csv")
    cluster_wirka <- read_csv(file)
    for(k in 1:length(names(table(cluster_wirka$Cluster)))){
    deg_wirka <- cluster_wirka %>% 
      filter(Cluster==names(table(cluster_wirka$Cluster))[k]) %>% 
      pull(gene) %>% 
      toupper()
    deg <- ClusterMarker %>% 
      filter(cluster==j) %>% arrange(p_val) %>%
      pull(gene) %>% head(100) %>% toupper()
    intersect_gene[i,k] <- length(intersect(deg,deg_wirka))
    #超几何分布分析
    p[i,k] <- phyper( intersect_gene[i,k]-1 ,100,20000-100,100,
                 lower.tail = F)
    
    }
  }
  list_p[j] = p
  list[[j]] = intersect_gene
  names(list[[j]]) = paste0("State",j)
}



#转换数据
#转换数据
for (j in 1:5) {
  data <- data.frame(celltype=colnames(list[[j]]),
                     mean=1,sd=1,lower_sd=1,
                     higher_sd=1)
  sd=c()
  for (i in 1:length(colnames(list[[j]]))) {
    data$mean[i] <-round(mean(list[[j]][,i]),2)
    sd[i] <- round(sd(list[[j]][,i]),2)
    data$lower_sd[i] <-  round(data$mean[i]-sd[i],2)
    data$higher_sd[i] <- round(data$mean[i]+sd[i],2)
  }
  data <- list[[j]] %>% t() %>% as.data.frame() %>% 
    rownames_to_column(var = "celltype") %>% 
    pivot_longer(.,cols =mouse_base:mouse_16w,
                 names_to ='time_point',
                 values_to = "num_gene") %>% 
    inner_join(data,by="celltype") %>% arrange(mean)
  class(data$celltype)
  table(data$celltype)
  data$celltype <- factor(data$celltype,
                          levels = unique(data$celltype))
  
  library(ggplot2)
  library(ggpubr)
  library(ggsci)
  list_p[[j]] <- ggplot(data = data,aes(x=celltype,y=num_gene)) +  
    geom_bar(data = data, 
             aes(x = celltype, y = mean,
                 fill=celltype),
             stat = "identity",
             width = 0.3, 
             position = position_dodge(width = 0.9)) +
    geom_jitter(width = 0.1,size=1)+
    geom_errorbar(data = data,aes(x=celltype,
                                  ymin =lower_sd,
                                  ymax = higher_sd),
                  width = 0.15,
                  position = position_dodge(width = 0.9))+
    scale_y_continuous(expand=c(0,0))+
    scale_fill_brewer(palette="Set3")+
    guides(fill=F)+xlab("celltype")+ylab("Number of genes")+
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5) )+
    labs(title =paste0("State ",j))+
    stat_compare_means(comparisons = my_comparisons,
                       method = "wilcox.test",
                       method.args = list(alternative = "greater"),
                       label = "p.signif"
    )+
    coord_flip()
}
library(patchwork)


#####不同得monocle state和modulated SMC交集分析#######
#构架矩阵，输出结果
intersect_gene_state <- matrix(nrow = 3,ncol = 5)
dataset <- c("mouse_base","mouse_8w","mouse_16w")
rownames(intersect_gene_state) <- dataset
colnames(intersect_gene_state) <- paste0("State ",1:5)
p <- intersect_gene_state
for (j in 1:5) {
  for (i in 1:3) {
    file <- paste0(paste0("TCF21_",dataset[i]),".csv")
    cluster_wirka <- read_csv(file)
    deg_wirka <- cluster_wirka %>% 
      filter(Cluster=="Modulated SMC") %>% 
      pull(gene) %>% 
      toupper()
    deg <- ClusterMarker %>% 
      filter(cluster==j) %>% arrange(p_val) %>%
      pull(gene) %>% head(100) %>% toupper()
    intersect_gene_state[i,j] <- length(intersect(deg,deg_wirka))
    #超几何分布分析
    p[i,j] <- phyper( intersect_gene_state[i,j]-1 ,100,20000-100,100,
                      lower.tail = F)
    
  }
}


my_comparisons <- list(
  c("State 2", "State 4"),c("State 4", "State 3")
)

data <- data.frame(State=paste0("State ",1:5),
                   mean=1,sd=1,lower_sd=1,
                   higher_sd=1)
sd=c()
for (i in 1:length(colnames(intersect_gene_state))) {
  data$mean[i] <-round(mean(intersect_gene_state[,i]),2)
  sd[i] <- round(sd(intersect_gene_state[,i]),2)
  data$lower_sd[i] <-  round(data$mean[i]-sd[i],2)
  data$higher_sd[i] <- round(data$mean[i]+sd[i],2)
}
data <- intersect_gene_state %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "State") %>% 
  pivot_longer(.,cols =mouse_base:mouse_16w,
               names_to ='time_point',
               values_to = "num_gene") %>% 
  inner_join(data,by="State") %>% arrange(mean)
class(data$State)
table(data$State)
data$State <- factor(data$State,
                     levels = unique(data$State))

library(ggplot2)
library(ggpubr)
library(ggsci)
p1 <- ggplot(data = data,aes(x=State,y=num_gene)) +  
  geom_bar(data = data, 
           aes(x = State, y = mean,
               fill=State),
           stat = "identity",
           width = 0.3, 
           position = position_dodge(width = 0.9)) +
  geom_jitter(width = 0.1,size=1)+
  geom_errorbar(data = data,aes(x=State,
                                ymin =lower_sd,
                                ymax = higher_sd),
                width = 0.15,
                position = position_dodge(width = 0.9))+
  scale_y_continuous(expand=c(0,0),limits = c(0,30))+
  scale_fill_brewer(palette="Set3")+
  guides(fill=F)+xlab("State")+ylab("Number of genes")+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,size = 10) )+
  labs(title =c("The intersection of different states with Modulated SMC"))+
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     method.args = list(alternative = "greater"),
                     label = "p.signif"
  )+
  coord_flip()
p1
library("cowplot")
pdf("04_State_celltype_intersection.pdf",12,8)


plot_grid(list_p[[1]],list_p[[2]],list_p[[3]],list_p[[4]],list_p[[5]], p1,
          labels = c("A", "B", "C","D","E","F"),
          ncol = 3, nrow = 2)
dev.off()



