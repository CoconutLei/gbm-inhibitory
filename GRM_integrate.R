# 221020 合并counts矩阵+创建Seurat对象+标准化处理
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(tidyverse)

options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")
rm(list=ls())
gc()

dir <- "D:/Desktop/GBM/00Overview/rdata/"
dir2 <- "D:/Desktop/GBM/00RawData/scCount_human/"

sample_names <- c("H01", "H02", "H03", "H04", "H05_1", "H05_2", 
                  "H06", "H07", "H08", "H09", "H10", "H11", 
                  "H12", "H13", "H14", "H15")

process_single_sample <- function(sample_id, data_dir) {
  data_path <- paste0(data_dir, sample_id, "/filtered_feature_bc_matrix/")
  count_mat <- Read10X(data.dir = data_path)
  new_cell_ids <- paste0(sample_id, "_", 1:ncol(count_mat))
  colnames(count_mat) <- new_cell_ids
  seurat_obj <- CreateSeuratObject(counts = count_mat, min.cells = 5)
  seurat_obj$sample_name <- sample_id
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & 
                         percent.mt < 10 & nCount_RNA < 30000)
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  return(seurat_obj)
}

seurat_list <- lapply(sample_names, function(x) {
  cat("正在处理样本：", x, "\n") # 打印进度，便于查看
  process_single_sample(sample_id = x, data_dir = dir2)
})
names(seurat_list) <- sample_names

list2env(seurat_list, envir = globalenv())

count_mat_list <- lapply(seurat_list, function(x) x@assays$RNA@counts)
save(list = names(count_mat_list), 
     file = paste0(dir, "00_Human_RAWData_Seperate.RData"),
     envir = as.environment(count_mat_list))

gbm_seurat_merged <- merge(seurat_list[[1]], 
                           y = seurat_list[-1], 
                           add.cell.ids = sample_names, 
                           project = "GBM_scRNA")

print(gbm_seurat_merged)

H01$Pathology <- "glioma"
H02$Pathology <- "glioma"
H03$Pathology <- "glioma"
H04$Pathology <- "glioma"
H05_1$Pathology <- "glioma"
H05_2$Pathology <- "glioma"
H06$Pathology <- "glioma"
H07$Pathology <- "glioma"
H08$Pathology <- "glioma"
H09$Pathology <- "glioma"
H10$Pathology <- "glioma"
H11$Pathology <- "glioma"
H12$Pathology <- "glioma"
H13$Pathology <- "glioma"
H14$Pathology <- "glioma"
H15$Pathology <- "glioma"


H01$WHO <- "3"
H02$WHO <- "4"
H03$WHO <- "4"
H04$WHO <- "4"
H05_1$WHO <- "4"
H05_2$WHO <- "4"
H06$WHO <- "4"
H07$WHO <- "3"
H08$WHO <- "4"
H09$WHO <- "4"
H10$WHO <- "4"
H11$WHO <- "4"
H12$WHO <- "4"
H13$WHO <- "3"
H14$WHO <- "4"
H15$WHO <- "4"

head(H02@meta.data)
head(H12@meta.data)
head(H01@meta.data)
head(H07@meta.data)
head(H13@meta.data)


#*******************************************************#
#Save both objects as RData object for quicker load
save(list=c("H01","H02","H03","H04","H05_1","H05_2","H06","H07",
            "H08","H09","H10","H11","H12","H13","H14","H15"), 
     file=paste(dir,"00_Human_Seperate.RData", sep=""))


#*******************************************************#
library(Matrix)  #自动加载的是Matrix1.5-1版本
library(Seurat) #SCT develop版本
library(patchwork)
library(dplyr)
library(ggplot2)

rm(list=ls())
gc()
dir <- "D:/Desktop/GBM/00Overview/rdata/"
dir2 <- "D:/Desktop/GBM/00RawData/scCount_human/"

load(file=paste(dir,"00_Human_Seperate.RData", sep=""))

# 10号暂且不加入分群
# 1 7 13 少突胶质细胞瘤
# AA
gbm.patient.list <- list(H01 = H01, H02 = H02, H03 = H03, H04 = H04, 
                         H05_1 = H05_1, H05_2 = H05_2, H06 = H06, H07 = H07, 
                         H08 = H08, H09 = H09, H11 = H11, H12 = H12, 
                         H13 = H13, H14 = H14, H15 = H15)

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = gbm.patient.list)

######### for RPCA
gbm.patient.list <- lapply(X = gbm.patient.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
#########
# conda install -c conda-forge r-base=4.2.3
# 通过调整k.anchor的值可以调整整合的强度，值越大，整合强度越大
gbm.anchors <- FindIntegrationAnchors(object.list = gbm.patient.list, 
                                      anchor.features = features,
                                      reduction = "rpca")
# save(gbm.anchors, file=paste(dir,"01_Human_Anchors_rpca.RData", sep=""))
# load(file=paste(dir,"01_Human_Anchors_rpca.RData", sep=""))
# this command creates an 'integrated' data assay
gbm.combined <- IntegrateData(anchorset = gbm.anchors)
# save(gbm.combined, file=paste(dir,"01_250418SC_GBM_Integrated_rpca.RData", sep=""))

load(file=paste(dir,"01_250418SC_GBM_Integrated_rpca.RData", sep=""))
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(gbm.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
gbm.combined <- ScaleData(gbm.combined, verbose = FALSE)
gbm.combined <- RunPCA(gbm.combined, npcs = 30, verbose = FALSE)
ElbowPlot(gbm.combined)
gbm.combined <- RunUMAP(gbm.combined, reduction = "pca", dims = 1:20)
gbm.combined <- RunTSNE(gbm.combined, dims = 1:20)
gbm.combined <- FindNeighbors(gbm.combined, reduction = "pca", dims = 1:20)
# save(gbm.combined, file=paste(dir,"01_Human_FindNeighbors.RData", sep=""))
# load(file=paste(dir,"01_Human_FindNeighbors.RData", sep=""))
gbm.combined <- FindClusters(gbm.combined, resolution = 0.8)

# Visualization
p1 <- DimPlot(gbm.combined, reduction = "umap", group.by = "sample_name", raster=FALSE)
p2 <- DimPlot(gbm.combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, raster=FALSE, cols = colplot)
p3 <- DimPlot(gbm.combined, reduction = "umap", group.by = "orig.ident", raster=FALSE)
p2
p4 <- DimPlot(gbm.combined, reduction = "umap", group.by = "sample_name", raster=FALSE)
p5 <- DimPlot(gbm.combined, reduction = "umap", group.by = "seurat_clusters",label = TRUE, repel = TRUE, raster=FALSE)
p6 <- DimPlot(gbm.combined, reduction = "umap", group.by = "orig.ident", raster=FALSE)

p1 + p2
p4 + p5
p3 + p6

DimPlot(gbm.combined, reduction = "umap", split.by = "sample_name",raster=FALSE)

dir4 <- "D:/Desktop/GBM/03Myeloid/rdata/"
load(file=paste(dir4,"00_Microglia_aa.Rdata", sep=""))
dir5 <- "D:/Desktop/GBM/02TandNK/rdata/"
load(file=paste(dir5,"00_TandNK_aa.Rdata", sep=""))
dir6 <- "D:/Desktop/GBM/04Fibroblast/rdata/"
load(file=paste(dir6,"00_Fibroblast_aa.Rdata", sep=""))
dir8 <- "D:/Desktop/GBM/05Endothelial/rdata/"
load(file=paste(dir8,"00_Endothelial_aa.Rdata", sep=""))
dir9 <- "D:/Desktop/GBM/06Neuron/rdata/"
load(file=paste(dir9,"00_Neuron_aa.Rdata", sep=""))
dir10 <- "D:/Desktop/GBM/07Oligodendrocyte/rdata/"
load(file=paste(dir10,"00_Oligodendrocyte_aa.Rdata", sep=""))

DimPlot(gbm.combined, reduction = "umap", cells.highlight = WhichCells(microglia.aa))
DimPlot(gbm.combined, reduction = "umap", cells.highlight = WhichCells(tandNK.aa))
DimPlot(gbm.combined, reduction = "umap", cells.highlight = WhichCells(tandNK.aa))
DimPlot(gbm.combined, reduction = "umap", cells.highlight = WhichCells(endothelial.aa))
DimPlot(gbm.combined, reduction = "umap", cells.highlight = WhichCells(neuron.aa))
DimPlot(gbm.combined, reduction = "umap", cells.highlight = WhichCells(oligodendro.aa))


dir7 <- "D:/Desktop/GBM/01Tumor/rdata/"
load(file=paste(dir7,"00_tumor.gbm_aa.Rdata", sep=""))
load(file=paste(dir7,"00_tumor.OPC_aa.Rdata", sep=""))
load(file=paste(dir7,"00_tumor.Neuron_aa.Rdata", sep=""))
load(file=paste(dir7,"00_tumor.VEGFA_aa.Rdata", sep=""))
load(file=paste(dir7,"00_tumor.gbm.OPC.Rdata", sep=""))
dir9 <- "D:/Desktop/GBM/06Neuron/rdata/"
load(file=paste(dir9,"00_Neuron.TumorAssociated_aa.Rdata", sep=""))

DimPlot(gbm.combined, reduction = "umap", cells.highlight = WhichCells(tumor.gbm))
DimPlot(gbm.combined, reduction = "umap", cells.highlight = WhichCells(tumor.OPC))
DimPlot(gbm.combined, reduction = "umap", cells.highlight = WhichCells(tumor.Neuron))
DimPlot(gbm.combined, reduction = "umap", cells.highlight = WhichCells(tumor.VEGFA))
DimPlot(gbm.combined, reduction = "umap", cells.highlight = WhichCells(tumor.gbm.OPC))
DimPlot(gbm.combined, reduction = "umap", cells.highlight = WhichCells(Neuron.TumorAssociated))


# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(gbm.combined) <- "RNA"

# Cluster
aa.markers <- FindAllMarkers(gbm.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# dir3 <- "D:/Desktop/GBM/00Overview/result/"
# write.table(aa.markers, file=paste(dir3,"Human_FindAllMarkers_reso_1.5.csv", sep=""))
#top 10 markers
require(tidyverse)
require(dplyr)
top_50_markers <- aa.markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)
top_100_markers  <- aa.markers %>% group_by(cluster) %>% top_n(100, avg_log2FC)
genes_to_check <- unique(top_50_markers$gene)

# Create Dotplot 
# DotPlot(gbm.combined, features = genes_to_check) + coord_flip()
# DoHeatmap(gbm.combined, features = genes_to_check, raster = F) + NoLegend()
# Barplot of patients per cluster 

#Save Rdata
# save(gbm.combined, file=paste(dir,"02_250422Cluster_gbm.combined_rpca.Rdata", sep=""))
load(file=paste(dir,"02_250422Cluster_gbm.combined_rpca.Rdata", sep="")) #23155  6328
# save(aa.markers, file=paste(dir,"03_250422GBM_FindAllMarkers_rpca.Rdata", sep=""))
load(file=paste(dir,"03_250422GBM_FindAllMarkers_rpca.Rdata", sep=""))

VlnPlot(gbm.combined,
        features = c("nFeature_RNA", "nCount_RNA"),
        ncol = 3,split.by = 'seurat_clusters',pt.size = 0,y.max = NULL)

# T3/3Microglia/Astrocyte5/OPC3/oligodendrocyte1/Neuron/Fibroblasts3/Endothelial
genes_to_check = c()
# All on Dotplot 
p10 <- DotPlot(gbm.combined, features = genes_to_check, cluster.idents = T) + coord_flip()
ggsave(paste(dir2,"plot_out/S03/03integrate/31dotplot_of_gene_markers.pdf", sep=""), p10, width = 22, height = 10)

#T3/3Monocytes and macrophages/NK3
FeaturePlot(gbm.combined, reduction = "umap", features = c("CD3D","CD3G","CD3E", "CD68","CD163","CD14","FGFBP2", "KLRF1", "KLRD1"), raster=FALSE, min.cutoff = "q9", order=TRUE, label=TRUE)
FeaturePlot(gbm.combined, reduction = "umap", features = c("CD3D","CD3G","CD3E", "CD68","CD163","CD14"), raster=FALSE, min.cutoff = "q9", order=TRUE, label=TRUE)
VlnPlot(gbm.combined, features = c("CD3D","CD3G","CD3E", "CD68","CD163","CD14","FGFBP2", "KLRF1", "KLRD1"),pt.size = 0) 
#3Plasma/B3/3Fibroblasts
FeaturePlot(gbm.combined, reduction = "umap", features = c("MZB1","SDC1","CD79A", "CD19","CD79A","MS4A1",  "FGF7","LUM","DCN"), raster=FALSE,
            min.cutoff = "q9", order=TRUE, label=TRUE)
VlnPlot(gbm.combined, features = c("MZB1","SDC1","CD79A", "CD19","CD79A","MS4A1",  "FGF7","LUM","DCN"),pt.size = 0)
#epi or tumor5/Endothelial cells3/stromal1
FeaturePlot(gbm.combined, reduction = "umap", features = c("EPCAM","KRT19","PROM1","ALDH1A1","CD24", "PECAM1","VWF","RAMP2", "MME"), raster=FALSE,
            min.cutoff = "q9", order=TRUE, label=TRUE)
FeaturePlot(gbm.combined, reduction = "umap", features = c("PECAM1","VWF","RAMP2"), raster=FALSE,
            min.cutoff = "q9", order=TRUE, label=TRUE)

VlnPlot(gbm.combined, features = c("EPCAM","KRT19","PROM1","ALDH1A1","CD24", "PECAM1","VWF","RAMP2", "MME"),pt.size = 0)
#Mast2/DC2/Myeloid2/Melanocytes2
FeaturePlot(gbm.combined, reduction = "umap", features = c("TPSB2","TPSAB1", "FCER1A", "CD1C", "AIF1","LYZ", "MLANA","PMEL"), 
            min.cutoff = "q9", order=TRUE, label=TRUE, raster=FALSE)
VlnPlot(gbm.combined, features = c("TPSB2","TPSAB1", "FCER1A", "CD1C", "AIF1","LYZ", "MLANA","PMEL"),pt.size = 0)
#Neutrophil6/liver
FeaturePlot(gbm.combined, reduction = "umap", features = c("FCGR3B", "CXCR1", "CXCR2", "CMTM2","PROK2", "MMP25"), 
            min.cutoff = "q9", order=TRUE, label=TRUE,raster=FALSE)
VlnPlot(gbm.combined, features = c("FCGR3B", "CXCR1", "CXCR2", "CMTM2","PROK2", "MMP25"),pt.size = 0)

FeaturePlot(gbm.combined, reduction = "umap", features = c("SNAP25", "GAD1", "SLC17A7"), 
            min.cutoff = "q9", order=TRUE, label=TRUE,raster=FALSE)
VlnPlot(gbm.combined, features = c("SNAP25", "GAD1", "SLC17A7"),pt.size = 0)

FeaturePlot(gbm.combined, reduction = "umap", features = c("HBB", "HBN"), 
            min.cutoff = "q9", order=TRUE, label=TRUE,raster=FALSE)
VlnPlot(gbm.combined, features = c("HBB", "HBN"),pt.size = 0)

# 四种状态GBM
# NPC-like
FeaturePlot(gbm.combined, reduction = "umap", features = c("CD24"), 
            min.cutoff = "q9", order=TRUE, label=TRUE,raster=FALSE)
VlnPlot(gbm.combined, features = c("CD24"),pt.size = 0)


# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
# 进一步确定亚群特征
dedoub <- subset(gbm.combined, subset = seurat_clusters == "10")

DefaultAssay(dedoub) <- "integrated"
dedoub <- ScaleData(dedoub)
# dedoub <- FindVariableFeatures(object = dedoub)
dedoub <- RunPCA(dedoub, features = VariableFeatures(object = dedoub))
ElbowPlot(dedoub)
dedoub <- RunUMAP(dedoub, reduction = "pca", dims = 1:20,verbose = FALSE)
dedoub <- Runumap(dedoub, dims = 1:20)
dedoub <- FindNeighbors(dedoub, dims = 1:20, verbose = FALSE)
# save(dedoub, file=paste(dir,"04Dedoub_pre_26.Rdata", sep=""))
# load(file=paste(dir,"04Dedoub_pre_26.Rdata", sep=""))
dedoub <- FindClusters(dedoub, resolution = 0.8, verbose = FALSE)

p4 <- DimPlot(dedoub, reduction = "umap", group.by = "sample_name")
p5 <- DimPlot(dedoub, reduction = "umap", group.by = "seurat_clusters",label = TRUE, repel = TRUE)
p4 + p5
p6 <- DimPlot(dedoub, reduction = "umap", group.by = "orig.ident")

DefaultAssay(dedoub) <- "RNA"
DimPlot(gbm.combined, reduction = "umap", cells.highlight = WhichCells(dedoub, idents = 3))

microglia <- subset(dedoub,subset = seurat_clusters != "2")
dim(microglia) #31878  1225
microglia_cell_list <- colnames(microglia)

# save(microglia_cell_list, file=paste(dir,"05_250510GBM_microglia_cell_list.Rdata", sep=""))
load(file=paste(dir,"05_250510GBM_microglia_cell_list.Rdata", sep=""))

# write.csv(panNeuron_cell_list,"230318_gbm.combined_panNeuron.csv")
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

Idents(gbm.combined) <- "celltype"
FeaturePlot(gbm.combined, reduction = "umap", 
            features = c("GABRA1","GABRA2","GABRA3","GABRA4",
                         "GABRA5","GABRA6","GABRE","GABRB1"
            ), 
            min.cutoff = "q9", order=TRUE, label=TRUE,raster=FALSE)
FeaturePlot(gbm.combined, reduction = "umap", 
            features = c("GABRA1","GABRA2","GABRA3","GABRA4",
                         "GABRB1","GABRG2"), 
            min.cutoff = "q9", order=TRUE, label=TRUE,raster=FALSE)
VlnPlot(gbm.combined, features = c("GABRA1","GABRA2","GABRA3","GABRA4",
                                   "GABRB1","GABRG2"),pt.size = 0)

FeaturePlot(gbm.combined, reduction = "umap", 
            features = c("GABRG2","GABBR1","GABRR2","GABRR3",
                         "GABRP","GABRD","GABRP","GABRQ","GABRE"), 
            min.cutoff = "q9", order=TRUE, label=TRUE,raster=FALSE)


# T_cell = c(39) 
# # Microglia = c(11,20,30,40,44,28) 
# Microglia = c(11,20,30) 
# Endothelial = c(34,47)
# Fibroblast = c(37)
# Oligodendrocyte = c(14,27,45,46)
# Excitatory_Neuron = c(25,33,35)
# Inhibitory_Neuron = c(32,41)
# # Tumor_Astrocyte = c(0,1,3,5,6,10,12,14, 9)
# CNV = c(13,19,23,15,48,21,31,22,17,12,10,4,38,3,16,36,24,26,0,9,2,5,18,
#         7,6,29,43,8,1,42)
# Neuron = c(21)


# 250510分群
# T_cell = c(23) 
# # Microglia = c(11,20,30,40,44,28) 
# Microglia = c(14,18) 
# Endothelial = c(20,27)
# Fibroblast = c(22)
# Oligodendrocyte = c(13,17,26)
# Excitatory_Neuron = c(16,21,24)
# Inhibitory_Neuron = c(19)
# # Tumor_Astrocyte = c(0,1,3,5,6,10,12,14, 9)
# Tumor = c(0,2,3,4,5,6,7,8,9,11,12,15,25)
# 
# 
# 
# current.cluster.ids <- c(T_cell,
#                          Microglia,
#                          Endothelial,
#                          Fibroblast,
#                          Oligodendrocyte,
#                          Excitatory_Neuron,
#                          Inhibitory_Neuron,
#                          Tumor
# )
# 
# new.cluster.ids <- c(rep("T_cell",length(T_cell)),
#                      rep("Microglia",length(Microglia)),
#                      rep("Endothelial",length(Endothelial)),
#                      rep("Fibroblast",length(Fibroblast)),
#                      rep("Oligodendrocyte",length(Oligodendrocyte)),
#                      rep("Excitatory_Neuron",length(Excitatory_Neuron)),
#                      rep("Inhibitory_Neuron",length(Inhibitory_Neuron)),
#                      rep("Tumor",length(Tumor))
# )

# 250510分群--外加一个GABA receptor tumor
# T_cell = c(23) 
# # Microglia = c(11,20,30,40,44,28) 
# Microglia = c(14,18) 
# Endothelial = c(20,27)
# Fibroblast = c(22)
# Oligodendrocyte = c(13,17,26)
# Excitatory_Neuron = c(16,21,24)
# Inhibitory_Neuron = c(19)
# # Tumor_Astrocyte = c(0,1,3,5,6,10,12,14, 9)
# Tumor = c(0,2,3,4,5,6,7,8,9,11,12,25)
# Tumor_GABRB1 = c(15)

T_cell = c(23) 
# Microglia = c(11,20,30,40,44,28) 
Microglia = c(14,18) 
Endothelial = c(20,27)
Fibroblast = c(22)
Oligodendrocyte = c(13,17,26)
Excitatory_Neuron = c(16,21,24)
Inhibitory_Neuron = c(19)
# Tumor_Astrocyte = c(0,1,3,5,6,10,12,14, 9)
Tumor = c(0,2,4,5,6,7,8,9,11,12,15)
Tumor_proliferation = c(3,25)

current.cluster.ids <- c(T_cell,
                         Microglia,
                         Endothelial,
                         Fibroblast,
                         Oligodendrocyte,
                         Excitatory_Neuron,
                         Inhibitory_Neuron,
                         Tumor,
                         Tumor_proliferation
)

new.cluster.ids <- c(rep("T_cell",length(T_cell)),
                     rep("Microglia",length(Microglia)),
                     rep("Endothelial",length(Endothelial)),
                     rep("Fibroblast",length(Fibroblast)),
                     rep("Oligodendrocyte",length(Oligodendrocyte)),
                     rep("Excitatory_Neuron",length(Excitatory_Neuron)),
                     rep("Inhibitory_Neuron",length(Inhibitory_Neuron)),
                     rep("Tumor",length(Tumor)),
                     rep("Tumor_proliferation",length(Tumor_proliferation))
)


head(gbm.combined@meta.data)

gbm.combined@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(gbm.combined@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)

head(gbm.combined@meta.data)
table(gbm.combined@meta.data$celltype)

1                10       Endothelial Excitatory_Neuron 
13736              7018              1425              5609 
Fibroblast Inhibitory_Neuron         Microglia   Oligodendrocyte 
1055              2276              6627              7129 
T_cell             Tumor 
993            105795

# 28                40                44               CNV 
# 1886               995               616            120479 
# Endothelial Excitatory_Neuron        Fibroblast Inhibitory_Neuron 
# 1425              4652              1063              2276 
# Microglia   Oligodendrocyte            T_cell 
# 10139              7136               996
# 
# CNV       Endothelial Excitatory_Neuron 
# 120479              1425              4652 
# Fibroblast Inhibitory_Neuron         Microglia 
# 1063              2276             13636 
# Oligodendrocyte            T_cell 
# 7136               996

# gbm.combined$celltype[gbm.combined$cell_id %in% tumor_P01_list] <- "Tumor"


gbm.combined$cell_id <- colnames(gbm.combined)

gbm.combined$celltype[gbm.combined$cell_id %in% microglia_cell_list] <- "Microglia"
table(gbm.combined$celltype)
1                10       Endothelial Excitatory_Neuron 
13736               766              1425              5609 
Fibroblast Inhibitory_Neuron         Microglia   Oligodendrocyte 
1055              2276             12879              7129 
T_cell             Tumor 
993            105795 
colplot <- c("#87A2AC","#62B99D","#E5C956","#2470A4",
             "#ADD269","#D0B7D4","#EC893E","#A3C9DE","#55488E")
# colplot <- c("#D0B7D4","#62B99D","#2470A4","#E5C956",
#              "#ADD269","#EC893E")

DimPlot(gbm.combined, reduction = "umap", cells.highlight = WhichCells(gbm.combined, idents = 10))

gbm.combined@meta.data$QC <- rownames(gbm.combined@meta.data)
plotCB=as.data.frame(gbm.combined@meta.data %>% filter(celltype != "1" & celltype != "10"))[,"QC"]
p9 = DimPlot(gbm.combined, reduction = "umap", label = F, group.by = "celltype", pt.size=0.1,cells = plotCB,raster=FALSE, cols = colplot)
p9
p10 = DimPlot(gbm.combined, reduction = "umap", label = T, group.by = "celltype", pt.size=0.1,raster=FALSE, cols = colplot)
p10
# p9

dir2 <- "D:/Desktop/GBM/250510/"
ggsave(paste(dir2,"UMAP_GBM.pdf", sep=""), p9, width = 8, height = 6)

p11 = DimPlot(gbm.combined, reduction = "tsne", label = F, group.by = "celltype", pt.size=0.1,cells = plotCB,raster=FALSE, cols = colplot)
p11


#--------------------------保存分群之后的Rdata------------------------------
# dim(gbm.combined) # 31878 115865
# length(plotCB) # 107821
# gbm.combined.prefilter <- gbm.combined
# save(gbm.combined.prefilter, file=paste(dir,"05_413ClusterDefined_aa_prefilter.Rdata", sep=""))
# load(file=paste(dir,"05_413ClusterDefined_aa_prefilter.Rdata", sep=""))

gbm.combined <- subset(gbm.combined, subset = celltype !="28" & celltype !="40" & celltype !="44")
dim(gbm.combined) # 31878 107821
# 31878 115726


library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")
library(Matrix)  
library(dplyr)
library(tidyverse)

# rm(list=ls())
gc()
dir <- "D:/Desktop/GBM/00Overview/rdata/"
dir2 <- "D:/Desktop/GBM/00RawData/scCount_human/"

dim(gbm.combined) #32161 151663
gbm.combined <- gbm.combined[,plotCB]
dim(gbm.combined) #32161 137161

# save(gbm.combined, file=paste(dir,"05_250510ClusterDefined_GBM.Rdata", sep=""))
load(file=paste(dir,"05_250510ClusterDefined_GBM.Rdata", sep=""))
table(gbm.combined@meta.data$celltype)

Endothelial Excitatory_Neuron        Fibroblast Inhibitory_Neuron 
1425              5609              1055              2276 
Microglia   Oligodendrocyte            T_cell             Tumor 
12879              7129               993            105795

# save(gbm.combined, file=paste(dir,"05_250512Map_to_Spatial_GBM.Rdata", sep=""))
load(file=paste(dir,"05_250512Map_to_Spatial_GBM.Rdata", sep=""))


dir3 <- "D:/Desktop/GBM/250510/rdata/SingleCell/"
IN_neuron <- subset(gbm.combined, subset = celltype == "Inhibitory_Neuron")
dim(IN_neuron)
# save(IN_neuron, file=paste(dir3,"05_250611_IN_neuron.Rdata", sep=""))
load(file=paste(dir3,"05_250611_IN_neuron.Rdata", sep=""))


# CNV       Endothelial Excitatory_Neuron        Fibroblast 
# 120479              1425              4652              1063 
# Inhibitory_Neuron         Microglia   Oligodendrocyte            T_cell 
# 2276             13636              7136               996

FeaturePlot(gbm.combined, reduction = "tsne", 
            features = c("GABRA1","GABRA2","GABRA3","GABRA4",
                         "GABRB1","GABRG2",ncol = 3), 
            min.cutoff = "q9", order=TRUE, label=FALSE,raster=FALSE, cols = c("grey", "#903323"))
FeaturePlot(gbm.combined, reduction = "umap", 
            features = c("MKI67","TOP2A"), 
            min.cutoff = "q9", order=TRUE, label=FALSE,raster=FALSE, cols = c("grey", "#903323"))

p1 <- FeaturePlot(gbm.combined, reduction = "umap", 
                  features = c("GABRA1"), 
                  min.cutoff = "q9", order=TRUE, label=FALSE,raster=FALSE, cols = c("grey", "#903323"))
p2 <- FeaturePlot(gbm.combined, reduction = "umap", 
                  features = c("GABRA2"), 
                  min.cutoff = "q9", order=TRUE, label=FALSE,raster=FALSE, cols = c("grey", "#903323"))
p3 <- FeaturePlot(gbm.combined, reduction = "umap", 
                  features = c("GABRA3"), 
                  min.cutoff = "q9", order=TRUE, label=FALSE,raster=FALSE, cols = c("grey", "#903323"))
p4 <- FeaturePlot(gbm.combined, reduction = "umap", 
                  features = c("GABRA4"), 
                  min.cutoff = "q9", order=TRUE, label=FALSE,raster=FALSE, cols = c("grey", "#903323"))
p5 <- FeaturePlot(gbm.combined, reduction = "umap", 
                  features = c("GABRB1"), 
                  min.cutoff = "q9", order=TRUE, label=FALSE,raster=FALSE, cols = c("grey", "#903323"))
p6 <- FeaturePlot(gbm.combined, reduction = "umap", 
                  features = c("GABRG2"), 
                  min.cutoff = "q9", order=TRUE, label=FALSE,raster=FALSE, cols = c("grey", "#903323"))
p7 <- FeaturePlot(gbm.combined, reduction = "umap", 
                  features = c("MKI67"), 
                  min.cutoff = "q9", order=TRUE, label=FALSE,raster=FALSE, cols = c("grey", "#903323"))
p8 <- FeaturePlot(gbm.combined, reduction = "umap", 
                  features = c("TOP2A"), 
                  min.cutoff = "q9", order=TRUE, label=FALSE,raster=FALSE, cols = c("grey", "#903323"))

dir2 <- "D:/Desktop/GBM/250510/"
ggsave(paste(dir2,"01Featureplot_GABRA1.pdf", sep=""), p1, width = 7, height = 6)
ggsave(paste(dir2,"01Featureplot_GABRA2.pdf", sep=""), p2, width = 7, height = 6)
ggsave(paste(dir2,"01Featureplot_GABRA3.pdf", sep=""), p3, width = 7, height = 6)
ggsave(paste(dir2,"01Featureplot_GABRA4.pdf", sep=""), p4, width = 7, height = 6)
ggsave(paste(dir2,"01Featureplot_GABRB1.pdf", sep=""), p5, width = 7, height = 6)
ggsave(paste(dir2,"01Featureplot_GABRG2.pdf", sep=""), p6, width = 7, height = 6)
ggsave(paste(dir2,"01Featureplot_MKI67.pdf", sep=""), p7, width = 7, height = 6)
ggsave(paste(dir2,"01Featureplot_TOP2A.pdf", sep=""), p8, width = 7, height = 6)



proliferation_score <- c("TOP2A","MKI67") %>% list()
gbm.combined <- AddModuleScore(object = gbm.combined, features = proliferation_score, name = 'Proliferation_score', assay = "RNA") 
Idents(gbm.combined) <- "celltype"
gbm.combined.tumor <- subset(gbm.combined,
                             subset = celltype == "Tumor" | celltype == "Tumor_GABRB1")
RidgePlot(gbm.combined.tumor, features = c("Proliferation_score1"), ncol = 1, cols = colplot)

plotdf <- gbm.combined.tumor@meta.data[,c("celltype","Proliferation_score1")]
# plotdf <- plotdf[plotdf$celltype == "TF" | plotdf$celltype == "TS", ]


#--------------------------保存分群之后的各subset---------------------------
microglia.gbm <- subset(gbm.combined,
                        subset = (celltype == "Microglia")) #32161 13636
tandNK.gbm <- subset(gbm.combined,
                     subset = (celltype == "T_cell")) #32161   996
fibroblast.gbm <- subset(gbm.combined,
                         subset = (celltype == "Fibroblast")) #32161  1063
tumor.gbm <- subset(gbm.combined,
                    subset = (celltype == "CNV")) #32161 120479
endothelial.gbm <- subset(gbm.combined,
                          subset = (celltype == "Endothelial")) #32161  1425
neuron.gbm <- subset(gbm.combined,
                     subset = (celltype == "Excitatory_Neuron" | celltype == "Inhibitory_Neuron")) #32161  6928
oligodendro.gbm <- subset(gbm.combined,
                          subset = (celltype == "Oligodendrocyte")) #32161  7136


dir4 <- "D:/Desktop/GBM/03Myeloid/rdata/"
save(microglia.gbm, file=paste(dir4,"00_Microglia_gbm.Rdata", sep=""))
dir5 <- "D:/Desktop/GBM/02TandNK/rdata/"
save(tandNK.gbm, file=paste(dir5,"00_TandNK_gbm.Rdata", sep=""))
dir6 <- "D:/Desktop/GBM/04Fibroblast/rdata/"
save(fibroblast.gbm, file=paste(dir6,"00_Fibroblast_gbm.Rdata", sep=""))
dir7 <- "D:/Desktop/GBM/01Tumor/rdata/"
save(tumor.gbm, file=paste(dir7,"00_Tumor_gbm.Rdata", sep=""))
dir8 <- "D:/Desktop/GBM/05Endothelial/rdata/"
save(endothelial.gbm, file=paste(dir8,"00_Endothelial_gbm.Rdata", sep=""))
dir9 <- "D:/Desktop/GBM/06Neuron/rdata/"
save(neuron.gbm, file=paste(dir9,"00_Neuron_gbm.Rdata", sep=""))
dir10 <- "D:/Desktop/GBM/07Oligodendrocyte/rdata/"
save(oligodendro.gbm, file=paste(dir10,"00_Oligodendrocyte_gbm.Rdata", sep=""))





#-------------------------继分群之后的进一步过滤----------------------------

gbm.combined.backup <- gbm.combined
gbm.combined$cell_id <- row.names(gbm.combined@meta.data)

gbm.combined$celltype[gbm.combined$cell_id %in% microglia.celist] <- "RealMicroglia"
table(gbm.combined$celltype)
plotCB=as.data.frame(gbm.combined@meta.data %>% filter(celltype !="26" &
                                                         celltype !="Microglia"))[,"QC"]
# 5724
p9 = DimPlot(gbm.combined, reduction = "umap", label = T, group.by = "celltype", pt.size=0.1,cells = plotCB,raster=FALSE)

gbm.combined <- gbm.combined[,plotCB]
# gbm.combined$celltype[gbm.combined$celltype == "RealMicroglia"] <- "Microglia"
# table(gbm.combined$celltype)
# Astrocyte     Endothelial     Fibroblasts       Microglia          Neuron Oligodendrocyte 
# 116651            1414             907            9624            6328            4874 
# T_and_NK 
# 1265
#最终保存的应该是去掉无关细胞的这张
p8 = DimPlot(gbm.combined, reduction = "umap", label = T, group.by = "celltype", pt.size=0.1,cells = plotCB,raster=FALSE)
save(gbm.combined, file=paste(dir,"06_Filter1_Human.Rdata", sep=""))
load(file=paste(dir,"06_Filter1_Human.Rdata", sep=""))



