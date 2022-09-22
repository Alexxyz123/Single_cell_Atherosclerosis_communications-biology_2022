####计算组间细胞百分比差异
tmp <-dplyr::select(scRNA@meta.data,c("orig.ident", "celltype"))
df <- data.frame()
for(i in unique(tmp$orig.ident)){
  df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(celltype) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("sample", "celltype", "value")
  df <- rbind(df, df_i)
}
df_i=spread(df,key = sample,value = value)
df_i=as.data.frame(df_i)
options(digits=2)
patient <- tmp$orig.ident %>% unique() %>% as.character()
for (i in 1:6) {
  total <- df_i[1:8,patient[i]] %>%as.vector() %>%sum()
  for (j in 1:8) {
    num <- df_i[j,patient[i]]
    df_i[j+8,patient[i]] <- num/total
  }
}
df_i[9:16,1] <- df_i[1:8,1]
df_i <- df_i[9:16,]
write.csv(df_i,file = "percent_overall.csv")
