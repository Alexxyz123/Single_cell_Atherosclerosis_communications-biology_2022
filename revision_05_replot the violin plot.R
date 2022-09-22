col_data <- data.frame(FetchData(object = scRNA,
                              vars = c("orig.ident",
                                       "seurat_clusters","celltype"))) %>% 
  rownames_to_column(var = "ID")

  
marker_gene <- c("TAGLN","CD34","PECAM1",
                 "CD3E","CD79A",
                 "CD68","CD1C",
                 "KLRB1","KIT")

exprs <- scRNA@assays$SCT@scale.data[marker_gene,]
exp <- exprs %>% t() %>% as.data.frame() %>% 
  rownames_to_column(var = "ID")%>% inner_join(col_data)
#获取更多颜色
getPalette <- colorRampPalette(brewer.pal(8,"Set1"))(19)
show_col(getPalette)




##画图顺序
list=list()
for(i in 1:9){
  exp$seurat_clusters <- factor(exp$seurat_clusters, 
                                levels = names(rev(sort(tapply(X = exp[,marker_gene[i]], INDEX = exp$seurat_clusters,FUN = mean),
                                                        decreasing = F))))
  
list[[i]] <- 
ggplot(exp,aes_string("seurat_clusters",marker_gene[i],
                      fill="seurat_clusters",
                      color="seurat_clusters"))+
  geom_violin(scale = "width", adjust =1, trim = TRUE,color="white")+ 
  geom_jitter(width = 0.1,size=0.1e-10)+
  scale_fill_manual(values=getPalette)+
  theme_minimal()+
  ggtitle(marker_gene[i])+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90,hjust = 1))+ 
  xlab(label = "")+
  guides(fill=FALSE,color=FALSE)+
  ylab(label = "Expression level")

}
pdf(file = "04_Cellmarker_type_gg.pdf",8,8)
plot_grid(plotlist = list,ncol = 3,
          align = "h")

dev.off()
