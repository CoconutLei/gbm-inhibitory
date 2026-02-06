# ssGSEA的包
library(ggplot2)
library(tinyarray)
library(GSVA)
library(dplyr)
library(Hmisc)
library(pheatmap)
library(ggpubr)
#输入表达矩阵 不要取log（标化后十几内可以），不可以有负值和缺失值---所有基因全表达矩阵
setwd("E:/CGGA")
expr <- read.table("CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)   
expr<- as.matrix(expr)    #注意将表达谱的data.frame转化为matrix，否则后续的gsva分析会报错
#输入基因集
cellMarker <- list("EN" = c("SEPT1","STMN1","RPH3A","GRIN2D","DLG4","CNIH2","HOMER2","PSD"),"IN" = c("GABRA4", "GABBR1", "GPHN"," GAD1", "NLGN2", "GAD2"))
#进行gsva分析
# 创建ssGSEA专属参数对象
params <- ssgseaParam(
  exprData = expr,              # 强制转换为矩阵格式
  geneSets = cellMarker,        # 基因集列表（必须为list格式）
  alpha = 0.25,                 # 默认权重参数
  normalize = TRUE,             # 启用结果标准化（推荐）
)
gsva_data <- gsva(params, verbose = TRUE)

# SINGLEGENEEXPR <- expr[c("BMP2"),,drop=F]
TOP10 <- expr[c("ARHGEF39","PODNL1","CD37","TOE1","MFAP2","HOXA10","KIF20A","STIL","DMRTA2","PARS2","OSR2"),,drop=F]

cor.test(as.numeric(OSR2),as.numeric(BMP2),alternative = c("less"),method=c("spearman"))