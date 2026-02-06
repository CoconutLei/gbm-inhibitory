
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

rm(list=ls())
file_paths <- c("D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2410025556_counts.txt", 
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2410025557_counts.txt", 
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2410025558_counts.txt",
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2410025559_counts.txt", 
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2410025560_counts.txt", 
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2410025561_counts.txt",
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2411017023_counts.txt",
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2411017024_counts.txt",
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2411017025_counts.txt",
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2411017028_counts.txt",
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2411017029_counts.txt", 
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2411017030_counts.txt", 
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2411017031_counts.txt", 
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2411017033_counts.txt", 
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2411017034_counts.txt",
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2411017035_counts.txt", 
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2411017037_counts.txt", 
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2411017038_counts.txt",
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2411017039_counts.txt", 
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2411017040_counts.txt", 
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2411017041_counts.txt"
)

count_data <- lapply(file_paths, function(file) {
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
                             "GG4","GG5","GG6", "GIN4","GIN5","GIN6",
                             "NIN1","NIN2","NIN3", "GEN1","GEN2","GEN3",
                             "NEN1","NEN2","NEN3")

keep <- rowSums(merged_counts > 3) >= 1
filtered_counts <- merged_counts[keep, ]

sample_info <- data.frame(
  sample_id = colnames(merged_counts),
  group = factor(c("GG","GG","GG", "GIN","GIN","GIN",
                   "GG","GG","GG", "GIN","GIN","GIN",
                   "NIN","NIN","NIN", "GEN","GEN","GEN",
                   "NEN","NEN","NEN"), 
                 levels = c("GG", "GIN","NIN","GEN","NEN"))
)
rownames(sample_info) <- sample_info$sample_id  

# 加载所需包
library(org.Mm.eg.db)
library(clusterProfiler)

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
  warning("duplicated")
  unique_indices <- !duplicated(gene_map$gene_symbol)
  gene_map <- gene_map[unique_indices, ]
  filtered_counts <- filtered_counts[unique_indices, ]
}
rownames(filtered_counts) <- gene_map$gene_symbol

head(filtered_counts)

library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts[, sample_info$group %in% c("GG", "GIN")],
  colData = sample_info[sample_info$group %in% c("GG", "GIN"), ],
  design = ~ group
)

dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "GG", "GIN"))
res <- res[!is.na(res$padj), ]  # 去除NA

head(res)
diff_genes_GG <- rownames(res)[res$padj < 0.05 & res$log2FoldChange > 0.5]
diff_genes_GIN <- rownames(res)[res$padj < 0.05 & res$log2FoldChange < -0.5]
length(diff_genes_GG)  # 605
length(diff_genes_GIN)  # 755

(1) enrichGO富集
# GG组
library(aPEAR)
head(diff_genes_GG)
go_enrich <- enrichGO(gene = diff_genes_GG, 
                      OrgDb = "org.Mm.eg.db", 
                      keyType = "SYMBOL", 
                      ont = "ALL")

head(go_enrich)

library(aPEAR)
library(dplyr)
go_enrich_df <- as.data.frame(go_enrich@result)
p <- enrichmentNetwork(go_enrich_df, drawEllipses = TRUE, fontSize = 4)
p

p1 <- p + scale_color_gradient(low = "lightblue", high = "steelblue")
p1
dir3 <- "D:/Desktop/GBM/250510/NewResult/Bulk_enrichment_bubble/"
ggsave(paste(dir3,"01enrichment1_GG.pdf", sep=""), p1, width = 8, height = 8)


# GIN组
head(diff_genes_GIN)
go_enrich <- enrichGO(gene = diff_genes_GIN, 
                      OrgDb = "org.Mm.eg.db", 
                      keyType = "SYMBOL", 
                      ont = "ALL")

head(go_enrich)

library(aPEAR)
library(dplyr)
go_enrich_df <- as.data.frame(go_enrich@result)
p <- enrichmentNetwork(go_enrich_df, drawEllipses = TRUE, fontSize = 4)
# p

p2 <- p + scale_color_gradient(low = "lightblue", high = "steelblue")
p2

dir3 <- "D:/Desktop/GBM/250510/NewResult/Bulk_enrichment_bubble/"
ggsave(paste(dir3,"01enrichment1_GIN.pdf", sep=""), p2, width = 8, height = 8)

plot_data_filtered <- go_enrich_df %>%
  filter(ID %in% target_ids) 

p <- enrichmentNetwork(
  plot_data_filtered,
  drawEllipses = TRUE, 
  fontSize = 4)
p

p + scale_color_viridis(discrete = FALSE, option = "D") 
p3 <- p + scale_color_gradient(low = "lightblue", high = "steelblue")
p3
dir3 <- "D:/Desktop/GBM/250510/NewResult/Bulk_enrichment_bubble/"
ggsave(paste(dir3,"01enrichment1_GIN3.pdf", sep=""), p3, width = 8, height = 8)


(2)gseGO富集
geneList = res[,"log2FoldChange"]
names(geneList) = as.character(row.names(res))
geneList = sort(geneList, decreasing = TRUE)
geneList <- geneList[names(geneList) !=""]
geneList <- geneList[!duplicated(names(geneList))]
head(geneList)
tail(geneList)
length(geneList)
geneList_GG <- geneList[geneList > 1]
length(geneList_GG)
geneList_GIN <- geneList[geneList < -1]
length(geneList_GIN)

ego <- gseGO(geneList = geneList_GG,
             OrgDb = org.Mm.eg.db,
             ont = "ALL",
             keyType = "SYMBOL",
             minGSSize = 10,
             maxGSSize = 500,
             eps = 0,
             pvalueCutoff = 0.05,
             verbose = FALSE) 

ego <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)


ego <- gseGO(geneList = geneList,
             OrgDb = org.Mm.eg.db,
             ont = "ALL",
             keyType = "SYMBOL",
             minGSSize = 10,
             maxGSSize = 500,
             eps = 0,
             pvalueCutoff = 0.05,
             verbose = FALSE) |> 
  clusterProfiler::simplify()

head(ego)

plot_data <- ego@result

p <- enrichmentNetwork(plot_data, drawEllipses = TRUE, fontSize = 4)
p


plot_data_filtered <- plot_data %>%
  filter(!ID %in% target_ids) 

p <- enrichmentNetwork(
  plot_data_filtered,
  drawEllipses = TRUE, 
  fontSize = 4)
p


p + scale_color_viridis(discrete = FALSE, option = "D") 
p3 <- p + scale_color_gradient(low = "lightblue", high = "steelblue")
p3
dir3 <- "D:/Desktop/GBM/250510/NewResult/Bulk_enrichment_bubble/"
ggsave(paste(dir3,"01enrichment2_GG.pdf", sep=""), p3, width = 8, height = 8)


# GIN组
ego <- gseGO(geneList = geneList_GIN,
             OrgDb = org.Mm.eg.db,
             ont = "ALL",
             keyType = "SYMBOL",
             minGSSize = 10,
             maxGSSize = 500,
             eps = 0,
             pvalueCutoff = 0.05,
             verbose = FALSE) 
ego <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
plot_data=ego@result

ego <- gseGO(geneList = geneList,
             OrgDb = org.Mm.eg.db,
             ont = "ALL",
             keyType = "SYMBOL",
             minGSSize = 10,
             maxGSSize = 500,
             eps = 0,
             pvalueCutoff = 0.05,
             verbose = FALSE) |> 
  clusterProfiler::simplify()

# 查看富集结果
head(ego)

plot_data <- ego@result

plot_data <- plot_data %>%
  mutate(
    enrichmentScore = abs(enrichmentScore),
    NES = abs(NES)
  )

p <- enrichmentNetwork(plot_data, drawEllipses = TRUE, fontSize = 4)
p

p + scale_color_viridis(discrete = FALSE, option = "D") 
p3 <- p + scale_color_gradient(low = "lightblue", high = "steelblue")
p3
dir3 <- "D:/Desktop/GBM/250510/NewResult/Bulk_enrichment_bubble/"
ggsave(paste(dir3,"01enrichment2_GIN_orange.pdf", sep=""), p, width = 8, height = 8)
ggsave(paste(dir3,"01enrichment2_GIN_blue.pdf", sep=""), p3, width = 8, height = 8)

