library(DESeq2)
library(tidyverse)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)  
library(rtracklayer)   
library(ggpubr)
library(extrafont)
library(dplyr)
loadfonts(device = "win") 

rm(list = ls())

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
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2411017030_counts.txt"
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
                             "GG4","GG5","GG6", "GIN4","GIN5","GIN6"
                             )

sample_info <- data.frame(
  sample_id = colnames(merged_counts),
  group = factor(c("GG","GG","GG", "GIN","GIN","GIN",
                   "GG","GG","GG", "GIN","GIN","GIN"), 
                 levels = c("GG", "GIN"))
)
rownames(sample_info) <- sample_info$sample_id  

gtf_path <- "D:/Desktop/GBM/250510/ref/gencode.vM38.annotation.gtf"  

gtf <- rtracklayer::import(gtf_path)

head(gtf)
gtf_gene <- gtf[gtf$type == "gene"]
gene_len_df <- data.frame(
  ensembl_id = gtf_gene$gene_id,          
  gene_length = width(gtf_gene),         
  stringsAsFactors = FALSE
) %>% distinct(ensembl_id, .keep_all = TRUE)  

head(gene_len_df)
head(merged_counts)

count_ensembl <- gsub("\\..*", "", rownames(merged_counts))  
gene_len_df$ensembl_id <- gsub("\\..*", "", gene_len_df$ensembl_id)

head(count_ensembl)
head(gene_len_df$ensembl_id)

gene_length <- gene_len_df$gene_length[match(count_ensembl, gene_len_df$ensembl_id)]
names(gene_length) <- rownames(merged_counts)
head(gene_length)

valid_len <- !is.na(gene_length)
counts_with_len <- merged_counts[valid_len, ]  
gene_length_filtered <- gene_length[valid_len] 
cat(nrow(counts_with_len))

tpm_matrix <- apply(counts_with_len, 2, function(sample_counts) {
  rpkm <- sample_counts / gene_length_filtered
  scaling_factor <- sum(rpkm) / 1e6
  tpm <- rpkm / scaling_factor
  return(tpm)
})
cat(dim(tpm_matrix)) 

keep <- rowSums(tpm_matrix > 1) >= 1
tpm_filtered <- tpm_matrix[keep, ]
cat(dim(tpm_filtered))


ensembl_ids_filtered <- gsub("\\..*", "", rownames(tpm_filtered))
gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = ensembl_ids_filtered,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first" 
)
gene_map <- data.frame(
  ensembl_id = rownames(tpm_filtered),
  gene_symbol = gene_symbols,
  stringsAsFactors = FALSE
)
valid_genes <- !is.na(gene_map$gene_symbol)
gene_map <- gene_map[valid_genes, ]
tpm_filtered <- tpm_filtered[valid_genes, ]

if (any(duplicated(gene_map$gene_symbol))) {
  warning("duplicated")
  unique_indices <- !duplicated(gene_map$gene_symbol)
  gene_map <- gene_map[unique_indices, ]
  tpm_filtered <- tpm_filtered[unique_indices, ]
}
rownames(tpm_filtered) <- gene_map$gene_symbol

counts_for_deseq <- merged_counts[rownames(counts_with_len)[keep], ]
dds <- DESeqDataSetFromMatrix(
  countData = counts_for_deseq,
  colData = sample_info,
  design = ~ group  
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "GIN", "GG"))
res_df <- as.data.frame(res) %>%
  rownames_to_column("ensembl_id") %>%
  mutate(ensembl_id = gsub("\\..*", "", ensembl_id)) %>%  
  left_join(gene_map, by = "ensembl_id") %>%             
  dplyr::select(gene_symbol, everything()) %>%           
  arrange(padj) 
sig_genes <- res_df %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)

sig_gene_symbols <- sig_genes$gene_symbol[!is.na(sig_genes$gene_symbol)]
if (length(sig_gene_symbols) == 0) {
  warning("skip")
} else {
  go_enrich <- enrichGO(
    gene = sig_gene_symbols,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",          
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  kegg_enrich <- enrichKEGG(
    gene = mapIds(org.Mm.eg.db, sig_gene_symbols, "ENTREZID", "SYMBOL", multiVals = "first"),
    organism = "mmu",
    pvalueCutoff = 0.05
  )
  go_plot <- dotplot(go_enrich, showCategory = 10, title = "GO Biological Process Enrichment") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("D:/Desktop/GBM/GO_enrichment.pdf", go_plot, width = 10, height = 6, dpi = 300)
  kegg_plot <- dotplot(kegg_enrich, showCategory = 10, title = "KEGG Pathway Enrichment") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("D:/Desktop/GBM/KEGG_enrichment.pdf", kegg_plot, width = 10, height = 6, dpi = 300)
}

pathway_list <- list(
  CAMKK_AMPK = c("Prkaa1", "Camkk2", "Irgm"),
  Calcium_Channel = c("P2rx4", "Trpm2", "Trpm5"),
  CREB_Signaling = c("Eif2ak4", "Sik1", "Ddit3")
)
pathway_scores <- map_dfc(names(pathway_list), function(pathway_name) {
  valid_genes <- intersect(pathway_list[[pathway_name]], rownames(tpm_filtered))
  if (length(valid_genes) == 0) {
    warning(paste("pathway", pathway_name, "NA！"))
    return(rep(NA, ncol(tpm_filtered)))
  }
  colMeans(tpm_filtered[valid_genes, , drop = FALSE])
}) %>%
  set_names(names(pathway_list)) %>%
  mutate(sample_id = colnames(tpm_filtered)) %>%
  left_join(sample_info, by = "sample_id") %>%
  filter(group %in% c("GG", "GIN")) 
if (nrow(pathway_scores) == 0) {
  stop("GG/GIN no pathway score")
}
heatmap_matrix <- pathway_scores %>%
  dplyr::select(-sample_id, -group) %>%
  as.matrix()
rownames(heatmap_matrix) <- pathway_scores$sample_id
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))
sample_anno <- pathway_scores %>%
  dplyr::select(sample_id, group) %>%
  column_to_rownames("sample_id")
pheatmap(
  mat = heatmap_matrix_scaled,
  annotation_col = sample_anno,
  annotation_colors = list(group = c(GG = "#2E86AB", GIN = "#A23B72")),
  scale = "none",          
  show_rownames = TRUE,    
  show_colnames = TRUE,    
  fontsize = 9,
  cellwidth = 35,
  cellheight = 15,
  color = colorRampPalette(c("#2166AC", "#F7FBFF", "#B2182B"))(100),
  main = "Pathway Score Heatmap (TPM)",
  filename = "D:/Desktop/GBM/Pathway_Heatmap_TPM.pdf",
  width = 9,
  height = 7,
  dpi = 300
)

write.csv(tpm_filtered, "D:/Desktop/GBM/TPM_filtered.csv", row.names = TRUE)
write.csv(res_df, "D:/Desktop/GBM/Differential_Genes.csv", row.names = FALSE)
if (nrow(sig_genes) > 0) {
  write.csv(sig_genes, "D:/Desktop/GBM/Significant_Genes.csv", row.names = FALSE)
}
