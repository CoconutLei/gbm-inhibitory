library(DESeq2)
library(tidyverse)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggpubr)
library(extrafont)
loadfonts(device = "win")
library(dplyr)

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
                "D:/Desktop/GBM/00RawData/smart_seq2/HTSeq_out/WH2411017030_counts.txt"
)

count_data <- lapply(file_paths, function(file) {
  data <- read.table(file, 
                     header = FALSE,  
                     sep = "\t",       
                     stringsAsFactors = FALSE,
                     skip = 0)         
  
  data <- data[!grepl("^__", data$V1), ]
  counts <- as.data.frame(data$V2)  
  rownames(counts) <- data$V1       
  colnames(counts) <- basename(file) %>% gsub("_counts.txt", "", .)  
  
  return(counts)
})

merged_counts <- do.call(cbind, count_data)
colnames(merged_counts) <- c("WH2410025556","WH2410025557","WH2410025558",
                             "WH2410025559","WH2410025560","WH2410025561",
                             "WH2411017023","WH2411017024","WH2411017025",
                             "WH2411017028","WH2411017029","WH2411017030")
dim(merged_counts)

ensembl_ids <- gsub("\\..*", "", rownames(merged_counts))
gene_length <- mapIds(
  org.Mm.eg.db,
  keys = ensembl_ids,
  keytype = "ENSEMBL",
  column = "TXLEN",  
  multiVals = "first"
)

names(gene_length) <- rownames(merged_counts)

merged_counts_with_len <- merged_counts[!is.na(gene_length), ]
gene_length_filtered <- gene_length[!is.na(gene_length)]

tpm_calc <- function(count_matrix, length_vector) {
  tpm_matrix <- apply(count_matrix, 2, function(sample_counts) {
    rpkm <- sample_counts / length_vector
    scaling_factor <- sum(rpkm) / 1e6
    tpm <- rpkm / scaling_factor
    return(tpm)
  })
  return(tpm_matrix)
}

tpm_matrix <- tpm_calc(merged_counts_with_len, gene_length_filtered)

keep <- rowSums(tpm_matrix > 1) >= 1
filtered_counts <- tpm_matrix[keep, ]  

sample_info <- data.frame(
  sample_id = colnames(filtered_counts),
  group = c("GG", "GG", "GG", "GIN", "GIN", "GIN",
            "GG", "GG", "GG", "GIN", "GIN", "GIN",
            "NIN", "NIN", "NIN", "GEN", "GEN", "GEN",
            "NEN", "NEN", "NEN")
)
rownames(sample_info) <- sample_info$sample_id

ensembl_ids <- gsub("\\..*", "", rownames(filtered_counts))
head(ensembl_ids)

gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = ensembl_ids,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"  
)
head(gene_symbols)

gene_map <- data.frame(
  ensembl_id = rownames(filtered_counts),
  gene_symbol = gene_symbols,
  stringsAsFactors = FALSE
)
head(gene_map)

valid_genes <- !is.na(gene_map$gene_symbol)
gene_map <- gene_map[valid_genes, ]
dim(gene_map)
head(filtered_counts)
head(gene_map)
filtered_counts <- filtered_counts[valid_genes, ]

if (any(duplicated(gene_map$gene_symbol))) {
  warning("duplicated")
  unique_indices <- !duplicated(gene_map$gene_symbol)
  gene_map <- gene_map[unique_indices, ]
  filtered_counts <- filtered_counts[unique_indices, ]
}
head(filtered_counts)
head(gene_map)
identical(rownames(filtered_counts),gene_map$ensembl_id)
rownames(filtered_counts) <- gene_map$gene_symbol

dim(filtered_counts)
head(filtered_counts)

CAMKK_AMPK_signaling_cascade <- c("Prkaa1","Zmpste24","Camkk2","Irgm")
ligand_gated_calcium_channel_activity <- c("P2rx4","Trpm2","Trpm5","Trpm4","Mcoln1","Trpm8")
negative_regulation_of_calcium_mediated_signaling <- c("Cd22","Itpr1","Plek","Sla2","Slc24a4")
positive_regulation_of_calcium_mediated_signaling <- c("Ada","Cd3e","Cd4","Cdh13","Neurod2","P2rx3","P2rx4",
                                                       "P2rx5","Plcg2","Ccl3","Ccl4","Syk","Zap70","P2rx2","Trat1",
                                                       "Trem2","Clec7a","Lacrt")
regulation_of_calcium_mediated_signaling <- c("Bst1","Calm1","Calm2","Calm3","Cmklr1","Ptk2b","Gbp1","Mapt",
                                              "Gpr143","Pdk2","Rit2","Rgn","Tmbim4","Tmem100","Selenok")
negative_regulation_of_CREB_transcription_factor_activity <- c("Eif2ak4","Sik1","Ddit3","Foxp3","Msx2","Adgrg3")

calcium_mediated_signaling <- c("Lmcd1","Pomc","Grin1","Akap6","Ackr3","Mctp1","Slc9a1","Tmem100")

positive_regulation_of_CREB_transcription_factor_activity <- c("Adcy1","Adcy8","Epha5","Cd200","Oprd1","Prkd1","Reln",
                                                               "Vegfa","Lrp8","Rps6ka4","Rps6ka5","Crtc1",
                                                               "Prkd2","Camk1d","Lpar5","Crtc3","Crtc2",
                                                               "Adgrf1","Tssk4","Syt14p1")

pathway_list <- list(
  CAMKK_AMPK_signaling_cascade = CAMKK_AMPK_signaling_cascade,
  ligand_gated_calcium_channel_activity = ligand_gated_calcium_channel_activity,
  calcium_mediated_signaling = calcium_mediated_signaling,
  negative_regulation_of_calcium_mediated_signaling = negative_regulation_of_calcium_mediated_signaling,
  positive_regulation_of_calcium_mediated_signaling = positive_regulation_of_calcium_mediated_signaling,
  regulation_of_calcium_mediated_signaling = regulation_of_calcium_mediated_signaling,
  negative_regulation_of_CREB_transcription_factor_activity = negative_regulation_of_CREB_transcription_factor_activity,
  positive_regulation_of_CREB_transcription_factor_activity = positive_regulation_of_CREB_transcription_factor_activity
)

pathway_scores <- map_dfc(names(pathway_list), function(pathway_name) {
  valid_genes <- intersect(pathway_list[[pathway_name]], rownames(filtered_counts))
  if (length(valid_genes) == 0) {
    warning(paste("pathway", pathway_name, "skip"))
    return(rep(NA, ncol(filtered_counts)))
  }
  colMeans(filtered_counts[valid_genes, , drop = FALSE])  
}) %>%
  set_names(names(pathway_list)) %>%
  mutate(sample_id = colnames(filtered_counts)) %>%
  left_join(sample_info, by = "sample_id") %>%
  filter(group %in% c("GG", "GIN"))

heatmap_matrix <- pathway_scores %>%
  dplyr::select(-sample_id, -group) %>%
  as.matrix()
rownames(heatmap_matrix) <- pathway_scores$sample_id
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))  

library(tibble)
sample_anno <- pathway_scores %>%
  dplyr::select(sample_id, group) %>%
  tibble::column_to_rownames("sample_id")

output_path <- "D:/Desktop/GBM/241113/GG_vs_GIN_Ca_Heatmap_TPM.tiff"  # 标注TPM
pheatmap(
  mat = heatmap_matrix_scaled,
  annotation_col = sample_anno,
  annotation_colors = list(group = c("GG" = "#2E86AB", "GIN" = "#A23B72")),
  scale = "none",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 9,
  cellwidth = 35,
  cellheight = 15,
  color = colorRampPalette(c("#2166AC", "#F7FBFF", "#B2182B"))(100),
  main = "GG vs GIN Pathway Score Heatmap (TPM)",  # 标题标注TPM
  filename = output_path,  # 用指定路径保存
  width = 9,
  height = 7,
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

output_path <- "D:/Desktop/GBM/241113/GG_vs_GIN_Pathway_Heatmap3_TPM.pdf"
pheatmap(
  mat = heatmap_matrix_scaled,
  annotation_col = sample_anno,
  annotation_colors = list(group = c(GG = "#2E86AB", GIN = "#A23B72")),  # 显式指定命名向量
  scale = "none",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 9,
  cellwidth = 35,
  cellheight = 15,
  color = colorRampPalette(c("#2166AC", "#F7FBFF", "#B2182B"))(100),
  main = "GG vs GIN Pathway Score Heatmap (TPM)",
  filename = output_path,
  width = 9,
  height = 7,
  dpi = 300,
  device = "pdf",
  compression = "lzw"
)

output_path <- "D:/Desktop/GBM/241113/GG_vs_GIN_Pathway_Heatmap4_TPM.pdf"

sorted_samples <- pathway_scores %>%
  arrange(group) %>%  
  pull(sample_id)     

heatmap_matrix_scaled_sorted <- heatmap_matrix_scaled[sorted_samples, ]

pheatmap(
  mat = heatmap_matrix_scaled_sorted,  
  annotation_col = sample_anno,
  annotation_colors = list(group = c(GG = "#2E86AB", GIN = "#A23B72")),
  scale = "none",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 9,
  cellwidth = 35,
  cellheight = 15,
  color = colorRampPalette(c("#2166AC", "#F7FBFF", "#B2182B"))(100),
  main = "GG vs GIN Pathway Score Heatmap (TPM)",
  filename = output_path,
  width = 9,
  height = 7,
  dpi = 300,
  device = "pdf",
  cluster_rows = FALSE, 
  cluster_cols = TRUE,   
  compression = "lzw"
)