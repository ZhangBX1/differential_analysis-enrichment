library(ggplot2)  
library(DESeq2)
library(clusterProfiler)
library(org.At.tair.db)
library(enrichplot)
countData <- as.matrix(read.csv("D://data//bac//4h.csv",row.names="Geneid")) # RNA-seq处理后导出的matrix文件

countData <- countData[rowMeans(countData)>1,]

condition <- factor(c(rep("mock",3),rep("infection",3)))

colData <- data.frame(row.names=colnames(countData), condition)

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE) 

res <- results(dds1)
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

res1_up<- res1[which(res1$log2FoldChange >= 2 & res1$pvalue < 0.05),]      # 表达量显著上升的基因
res1_down<- res1[which(res1$log2FoldChange <= -2 & res1$pvalue < 0.05),]    # 表达量显著下降的基因
res1_total <- rbind(res1_up,res1_down)

write.csv(res1_down,"D://data//bac//down.csv") # 下调基因导出的文件位置
write.csv(res1_up,"D://data//bac//up.csv") # 上调基因导出的文件位置
write.csv(res1_total,"D://data//bac//total.csv") # 差异基因导出的文件位置

genes <- read.table('D://data//bac//geneid.txt', stringsAsFactors = FALSE)[,1]

enrich.go <- enrichGO(gene = genes,  
                      OrgDb = org.At.tair.db,  
                      keyType = 'ENTREZID', 
                      ont = 'ALL',  
                      pAdjustMethod = 'fdr',  
                      pvalueCutoff = 0.05,  
                      qvalueCutoff = 0.2,  
                      readable = FALSE)

##可视化--点图
dotplot(enrich.go,title="EnrichmentGO_dot")#点图，按富集的数从大到小的
##可视化--条形图
barplot(enrich.go, showCategory=20,title="EnrichmentGO_MF")#条状图，按p从小到大排，绘制前20个Term
dotplot(enrich.go, showCategory=30) + ggtitle("dotplot for ORA")

edox2 <- pairwise_termsim(enrich.go)
p1 <- treeplot(edox2,hexpand=0.3,offset=15,offset_tiplab=1)
p2 <- treeplot(edox2, hclust_method = "average",hexpand=0.3,offset=15,offset_tiplab=1)
aplot::plot_list(p1, p2, tag_levels='A')

barplot(enrich.go, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图
barplot(KEGG,showCategory = 40,title = 'KEGG Pathway')
dotplot(enrich.go, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图
dotplot(KEGG)

enrich.kegg <- enrichKEGG(
  gene = genes,  #基因列表文件中的基因名称
  keyType = 'ncbi-geneid',  #KEGG 富集
  organism = 'ath',  #例如，ath 代表拟南芥，其它物种更改这行即可
  pAdjustMethod = 'fdr',  #指定 p 值校正方法
  pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
  qvalueCutoff = 0.2)  #指定 q 值阈值（可指定 1 以输出全部）

dotplot(enrich.kegg, showCategory=20)
