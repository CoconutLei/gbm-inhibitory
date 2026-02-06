# inferCNV persample
library(infercnv)
library(Seurat)
library(Matrix)
library(Cairo)
library(grDevices)

rm(list=ls())
dir <- "D:/Desktop/GBM/00Overview/rdata/"
dir2 <- "D:/Desktop/GBM/00RawData/scCount_human/"
dir3 <- "D:/Desktop/GBM/250510/rdata/SingleCell/inferCNV/"
# dir4 <- "D:/Desktop/GBM/01Tumor/result/01HumanInferCNV/"

load(file=paste(dir,"05_250510ClusterDefined_GBM.Rdata", sep=""))
table(gbm.combined$celltype)
# Endothelial Excitatory_Neuron        Fibroblast Inhibitory_Neuron         Microglia   Oligodendrocyte 
# 1425              5609              1055              2276             12879              7129 
# T_cell             Tumor 
# 993            105795
table(gbm.combined$sample_name)
# H01   H02   H03   H04 H05_1 H05_2   H06   H07   H08   H09   H11   H12   H13   H14   H15 
# 11643  4857 14093  6606 13695  7740  6911 13006  8134 10086 10681 10432  8264  8782  2231

gbm.combined.filter <- subset(gbm.combined, subset = sample_name != "H01" & sample_name != "H07" &
                                sample_name != "H13")
gbm1w <- gbm.combined.filter[,sample(ncol(gbm.combined.filter), size = 10000)]
# gbm1w <- subset(gbm.combined, subset = sample_name == "1w")
gbm_expresscount_1w <- as.matrix(gbm1w@assays$RNA@counts)
table(gbm1w$sample_name)
# H02   H03   H04 H05_1 H05_2   H06   H08   H09   H11   H12   H14   H15 
# 492  1353   646  1310   736   684   780   932  1031   969   830   237
#----------------------------------------------------------------------------------
#生成gene位置信息文件
library(AnnoProbe)
geneInfor=annoGene(rownames(gbm_expresscount_1w),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[,c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)
## 这里可以去除性染色体
# 也可以把染色体排序方式改变
# gbm_expresscount_1w=gbm_expresscount_1w[rownames(gbm_expresscount_1w) %in% geneInfor[,1],]
# gbm_expresscount_1w=gbm_expresscount_1w[match( geneInfor[,1], rownames(gbm_expresscount_1w) ),] 

dim(gbm_expresscount_1w)  # 32161 10000
write.table(x = gbm_expresscount_1w, file =  paste(dir3, "01inferCNV_gbm1w_ExpressCount.txt", sep = ""), quote = F,sep = '\t')
write.table(geneInfor,file = paste(dir3, "gbm1w_encode_gene_pos.txt", sep = ""),
            quote = F,col.names = F,row.names = F,sep = '\t')
#---------------------------------------------------------------------------

table(gbm.combined$celltype)
# Endothelial Excitatory_Neuron        Fibroblast Inhibitory_Neuron         Microglia 
# 1425              5609              1055              2276             12879 
# Oligodendrocyte            T_cell             Tumor 
# 7129               993            105795

table(gbm1w$celltype)
# Endothelial Excitatory_Neuron        Fibroblast Inhibitory_Neuron         Microglia 
# 92               272                52               109              1174 
# Oligodendrocyte            T_cell             Tumor 
# 484                79              7738

# Epithelial cells 
temp.1 <- gbm1w@meta.data[gbm1w@meta.data$celltype == "Tumor", c("celltype")] 
temp.1.rowname <- row.names(gbm1w@meta.data[gbm1w@meta.data$celltype == "Tumor",])

temp.1 <- cbind(as.data.frame(temp.1.rowname), as.data.frame(temp.1))
colnames(temp.1) <- c("V1","V2")
# row.names(temp.1) <- temp.1[,1]

# Endothelial cells and fibroblasts 
#*
# temp.2.cell_id <- gbm.combined@meta.data[gbm.combined@meta.data$celltype %in% "T_cell" & gbm.combined@meta.data$orig.ident == "MyData", c("cell_id")]
# temp.2.celltype <- gbm.combined@meta.data[gbm.combined@meta.data$celltype %in% "T_cell" & gbm.combined@meta.data$orig.ident == "MyData", c("celltype")]
# temp.2 <- cbind(as.data.frame(temp.2.cell_id), as.data.frame(temp.2.celltype))
# colnames(temp.2) <- c("V1","V2")
######################################到此开始运行

#************************用Microglia细胞为对照***************************
temp.3.cell_id <- gbm1w@meta.data[gbm1w@meta.data$celltype %in% "Microglia", c("cell_id")]
temp.3.celltype <- gbm1w@meta.data[gbm1w@meta.data$celltype %in% "Microglia", c("celltype")]
temp.3 <- cbind(as.data.frame(temp.3.cell_id), as.data.frame(temp.3.celltype))
colnames(temp.3) <- c("V1","V2")
# Keep a total of 500 cells from both 
# temp.2 <- temp.2[sample(nrow(temp.2), size = 198),]
# temp.3<- temp.3[sample(nrow(temp.3), size = 800),]

# Convert some of the fibroblasts and endos to spikeins for the clustering 
# temp.2$V2 <- paste(as.character(temp.2$V2), "_normal", sep="")
# temp.2$V2[sample(which(temp.2$V2=="endothelial_normal"),300)] <- "endothelial"
temp.3$V2 <- paste(as.character(temp.3$V2), "_normal", sep="")
# temp.3$V2[sample(which(temp.3$V2=="Fibroblasts_normal"),300)] <- "Fibroblasts"
# Combine 
inferCNV.annotation <- rbind(temp.1, temp.3)
rm("temp.1", "temp.3","temp.1.rowname","temp.3.cell_id","temp.3.celltype")
inferCNV.annotation[,1] <- paste(inferCNV.annotation[,1],inferCNV.annotation[,2], sep="\t")
gbm1w.inferCNV.annotation <- inferCNV.annotation[,-2]
# Write table
write.table(gbm1w.inferCNV.annotation, file = paste(dir3,"01inferCNV_annotation_gbm1w.txt", sep=""),row.names = F, col.names = F, quote=F, sep="\t")


# express_count <- read.table("D:/Desktop/GBM/01Tumor/rdata/01infercnv_gbm1w_ExpressCount.txt")
# annotation <- read.csv("D:/Desktop/GBM/01Tumor/rdata/01inferCNV_annotation_aaH01.txt")


#1
infercnv_gbm1w = CreateInfercnvObject(raw_counts_matrix="D:/Desktop/GBM/250510/rdata/SingleCell/inferCNV/01inferCNV_gbm1w_ExpressCount.txt",
                                     annotations_file="D:/Desktop/GBM/250510/rdata/SingleCell/inferCNV/01inferCNV_annotation_gbm1w.txt",
                                     delim="\t",
                                     gene_order_file="D:/Desktop/GBM/250510/rdata/SingleCell/inferCNV/gbm1w_encode_gene_pos.txt",
                                     ref_group_names=c("Microglia_normal") 
)


#2
infercnv_gbm1w = infercnv::run(infercnv_gbm1w,
                              cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              out_dir="D:/Desktop/GBM/250510/result/inferCNV2",
                              cluster_by_groups=TRUE, 
                              #analysis_mode="subclusters", #默认是"samples"
                              denoise=TRUE,
                              HMM=TRUE
                              #num_threads=4
)
save(infercnv_gbm1w, file = "D:/Desktop/GBM/250510/result/inferCNVinfercnv_gbm1w.Rdata")


#plot
# library(RColorBrewer)
# infercnv::plot_cnv(infercnv_meta,
#                    plot_chr_scale = T,
#                    output_filename = "inferCNV",output_format = "pdf",
#                    custom_color_pal =  color.palette(c("#972332","white","#14497E"), c(2, 2)))
# 
# infercnv::plot_cnv(infercnv_meta,
#                    plot_chr_scale = T,
#                    output_filename = "inferCNV",output_format = "pdf", c(2, 2))
# 
# saveRDS(object = infercnv_meta, file = "D:/Desktop/Sm2LungCanc/result/plot_out/S04/infer_CNV_new_FINAL.infercnv_meta")

