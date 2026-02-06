library(DESeq2)
library(aPEAR)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(ggplot2)
library(viridis)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(GSVA) 
library(GSEABase)
library(msigdbr)
library(limma)
library(BiocParallel)

setwd("E:/IN/smart-seq IN共培养/HTSeq_out")
dat <- read.csv("TPM_filtered.csv",header = T, row.names = 1)
file_path <- c("WH2410025556_counts.txt","WH2410025557_counts.txt","WH2410025558_counts.txt","WH2410025559_counts.txt","WH2410025560_counts.txt","WH2410025561_counts.txt","WH2411017023_counts.txt","WH2411017024_counts.txt","WH2411017025_counts.txt","WH2411017029_counts.txt","WH2411017030_counts.txt","WH2411017031_counts.txt")
count_data <- lapply(file_path, function(file) {
  data <- read.table(
    file, 
    header = FALSE, 
    sep = "\t", 
    stringsAsFactors = FALSE,
    quote = ""  
  )
  data <- data[!grepl("^__", data$V1), ]
  counts <- as.data.frame(data$V2)
  rownames(counts) <- data$V1
  colnames(counts) <- gsub("_counts.txt", "", basename(file))
  return(counts)
})
merged_counts <- do.call(cbind, count_data)
colnames(merged_counts) <- c("GG1","GG2","GG3", "GIN1","GIN2","GIN3",
                             "GG4","GG5","GG6", "GIN4","GIN5","GIN6")
# 去除低表达基因（例如，至少在一个样本中计数大于 10）
keep <- rowSums(merged_counts > 3) >= 1
filtered_counts <- merged_counts[keep, ]

sample_info <- data.frame(
  sample_id = colnames(merged_counts),
  group = factor(c("GG","GG","GG", "GIN","GIN","GIN",
                   "GG","GG","GG", "GIN","GIN","GIN"), 
                 levels = c("GG", "GIN"))
)
rownames(sample_info) <- sample_info$sample_id  

ensembl_ids <- gsub("\\..*", "", rownames(filtered_counts))

gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = ensembl_ids,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"  # 遇到多个匹配时取第一个
)

gene_map <- data.frame(
  ensembl_id = rownames(filtered_counts),
  gene_symbol = gene_symbols,
  stringsAsFactors = FALSE
)

valid_genes <- !is.na(gene_map$gene_symbol)
gene_map <- gene_map[valid_genes, ]
filtered_counts <- filtered_counts[valid_genes, ]

if (any(duplicated(gene_map$gene_symbol))) {
  warning("存在重复的基因名，将保留第一个出现的条目")
  unique_indices <- !duplicated(gene_map$gene_symbol)
  gene_map <- gene_map[unique_indices, ]
  filtered_counts <- filtered_counts[unique_indices, ]
}
rownames(filtered_counts) <- gene_map$gene_symbol

head(filtered_counts)

# 构建DESeq2对象
dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts[, sample_info$group %in% c("GG", "GIN")],
  colData = sample_info[sample_info$group %in% c("GG", "GIN"), ],
  design = ~ group
)

# 差异分析
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "GG", "GIN"))
res <- res[!is.na(res$padj), ]  # 去除NA

head(res)
# 筛选差异基因（例如padj < 0.05 & |log2FC| > 1）
# GG组高表达的（GG组表达的倍数是GIN组表达倍数的2倍及以上）
diff_genes_GG <- rownames(res)[res$padj < 0.05 & res$log2FoldChange > 0]
# GG组高表达的（GIN表达的倍数是GG组表达倍数的2倍及以上）
diff_genes_GIN <- rownames(res)[res$padj < 0.05 & res$log2FoldChange < 0]
length(diff_genes_GG)  # 605
length(diff_genes_GIN)  # 755
diff_genes_all <- rownames(res)[res$padj < 0.05]
head(diff_genes_GG)


go_enrich <- enrichGO(gene = diff_genes_sorted, 
                      OrgDb = "org.Mm.eg.db", 
                      keyType = "SYMBOL", 
                      ont = "ALL")

#gseGO富集
# GSEA分析的前提就是构造一个按照log2FoldChange降序排列的基因列表
# feature 1: numeric vector
geneList = res[,"log2FoldChange"]
names(geneList) = as.character(row.names(res))
geneList = sort(geneList, decreasing = TRUE)
geneList <- geneList[names(geneList) !=""] #删除基因列表中的空基因名geneList <- geneList[!duplicated(names(geneList))] #删除基因列表中的重复基因
geneList <- geneList[!duplicated(names(geneList))] #删除基因列表中的重复基因
head(geneList)
tail(geneList)
genelist_all <- names(geneList)

ego <- gseGO(geneList = geneList,
             OrgDb = org.Mm.eg.db,
             ont = "ALL",
             keyType = "SYMBOL",
             minGSSize = 10,
             maxGSSize = 500,
             eps = 0,
             pvalueCutoff = 0.05,
             verbose = FALSE) 

result_df <- as.data.frame(ego@result)
result_MGE_up <- result_df[result_df$NES < 0,]
result_MGE_down <- result_df[result_df$NES > 0,]
result1 <- result_MGE_up[-which(row.names(result_MGE_up)=="GO:0010959"),]

p <- enrichmentNetwork(result_MGE_up, drawEllipses = TRUE, fontSize = 4)
p2 <- p + scale_color_gradient(low = "steelblue", high = "lightblue")
print(p2)

p1 <- enrichmentNetwork(result1, drawEllipses = TRUE, fontSize = 4)

write.csv(result_df,file = "bubblenetwork.csv")
