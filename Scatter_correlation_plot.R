library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
rm(list = ls())

brainH15 <- Load10X_Spatial(data.dir = "D:/Desktop/GBM/00RawData/SpatialCount_human/SH15", # 存放cellranger结果的目录
                            filename = "filtered_feature_bc_matrix.h5",  # h5格式的表达矩阵名称
                            assay = "Spatial", # seurat对象中存放表达矩阵的assay名称
                            slice = "sliceH15"   # seurat对象中设置的空转样本ID，多样本合并分析时有用
)

### 数据质控
# 计算SPOT的线粒体基因比例
brainH15[["percent.mt"]] <- PercentageFeatureSet(brainH15, pattern = "^MT-")
# 计算SPOT的核糖体基因比例
brainH15[["percent.rb"]] <- PercentageFeatureSet(brainH15, pattern = "^RP[LS]")

# 降维聚类
# SPOTs聚类的结果无法与细胞类型对应起来，它通常与组织的不同结构区域或病理区域匹配度更高
brainH15 <- SCTransform(brainH15, assay = "Spatial") # 注意assay要指定为Spatial
brainH15 <- RunPCA(brainH15, assay = "SCT", verbose = FALSE)
brainH15 <- FindNeighbors(brainH15, reduction = "pca", dims = 1:30)
brainH15 <- FindClusters(brainH15, verbose = FALSE)
brainH15 <- RunUMAP(brainH15, reduction = "pca", dims = 1:30)

# dir <- "D:/Desktop/GBM/250510/rdata/"
# load(file=paste(dir,"01brainH15_ClusterDefined.Rdata", sep=""))

# Glutamatergic_presynaptic_marker <- c("EFNA1","EFNA2","EFNA3","EFNA4","EFNA5","EFNB1","EFNB2","EFNB3", #神经钙粘素
#                                       "SLC17A6","SLC17A7","SLC17A8", #把谷氨酸包装进囊泡
#                                       "GLS","GLS2" #合成谷氨酸
# ) %>% list()
# Glutamatergic_postsynaptic_marker <- c("GRIP1", #神经钙粘素受体
#                                        "NLGN1","NLGN3", #NLGN与Neurexin连接
#                                        "GRIA2","GRIA3","GRIA4", #AMPAR
#                                        "GRIN1","GRIN2B","GRIN2A", #NMDAR
#                                        "GRIK1","GRIK3","GRIK4","GRIK5",
#                                        "DLG4","HOMER1","HOMER2","HOMER3", #突触后致密区骨架
#                                        "SHANK1","SHANK2","SHANK3" #SHANK
# ) %>% list()
Glutamatergic_synaptic_marker <- c("EFNA1","EFNA2","EFNA3","EFNA4","EFNA5","EFNB1","EFNB2","EFNB3", #神经钙粘素
                                   "SLC17A6","SLC17A7","SLC17A8", #把谷氨酸包装进囊泡
                                   "GLS","GLS2", #合成谷氨酸
                                   "GRIP1", #神经钙粘素受体
                                   "NLGN1","NLGN3", #NLGN与Neurexin连接
                                   "GRIA2","GRIA3","GRIA4", #AMPAR
                                   "GRIN1","GRIN2B","GRIN2A", #NMDAR
                                   "GRIK1","GRIK3","GRIK4","GRIK5",
                                   "DLG4","HOMER1","HOMER2","HOMER3", #突触后致密区骨架
                                   "SHANK1","SHANK2","SHANK3" #SHANK
                                   
) %>% list()
# GABAergic_presynaptic_marker <- c("GAD1","GAD2") %>% list()
# GABAergic_postsynaptic_marker <- c("NLGN2","GABRA4","GABBR1","GPHN") %>% list()
# GABAergic_postsynaptic_marker <- c("GPHN","GABRA1","GABRA2","GABRA3",
#                                    "GABRA4","GABRA5","GABRB1",
#                                    "GABRB2","GABRB3") %>% list()
GABAergic_synaptic_marker <- c("GAD1","GAD2",
                               "GPHN","GABRA1","GABRA2","GABRA3",
                               "GABRA4","GABRA5","GABRB1",
                               "GABRB2","GABRB3") %>% list()
# Cholinergic_postsynaptic_marker <- c("CHRM1","CHRM2","CHRM3","CHRM4","CHRM5",
#                                      "CHRNA7","CHRNA4","CHRNB2") %>% list()
Cholinergic_synaptic_marker <- c("SLC18A3","CHAT","CHRNA4","CHRNA7","CHRNB4",
                                 "CHRM1","CHRM2","CHRM3","CHRM4","CHRM5",
                                 "CHRNA7","CHRNA4","CHRNB2") %>% list()
Proliferation_marker <- c("MKI67","PCNA","CDK1") %>% list()

brainH15 <- AddModuleScore(object = brainH15, features = Glutamatergic_synaptic_marker, name = 'Glutamatergic_synaptic', assay = "SCT") 
brainH15 <- AddModuleScore(object = brainH15, features = GABAergic_synaptic_marker, name = 'GABAergic_synaptic', assay = "SCT") 
brainH15 <- AddModuleScore(object = brainH15, features = Cholinergic_synaptic_marker, name = 'Cholinergic_synaptic', assay = "SCT") 
brainH15 <- AddModuleScore(object = brainH15, features = Proliferation_marker, name = 'Proliferation', assay = "SCT") 

# Proliferation_score <- c("MKI67","TOP2A","PCNA","MCM2","PROM1") %>% list()
Apoptosis <- c("BOK","PTGIS","TP53","TP53BP2","DLC1","SIRT2","ZC3H15A") %>% list()
# RTK_RAS_PI3K_pathway <- c("EGFR","PDGFRA","FGFR1","FGFR2","MET","VEGFR2","IGF1R",
#                           "KRAS","NRAS","HRAS","PIK3CA","PIK3R1","AKT1","AKT2",
#                           "AKT3","MTOR","PTEN") %>% list()
cell_cycle <- c("AHR","BIRC5","ARF6","CCNH","CCNT1","CCNT2","CDK7","CDKN1C","CKS1B",
                "CLTA","MAP3K8","CREBL2","CSNK1A1","CSNK2A1","CSNK2A2","CYLD",
                "BRINP1","DDIT3","DUSP1","E2F2","E4F1","EP300","EPB41","EPB41L2",
                "ERH","GAK","GAS1","GAS2","GNAI1","GNAI2","GNAI3","GPER1","NR3C1",
                "HCFC1","HELLS","HNRNPU","ING1","INSM1","KRT18","LIG3","LIG4","MKI67",
                "MLF1","MYOG","NEDD9","NEK1","PIK3C3","PIM1","PIN1","PPP1CA","PPP1CB",
                "PPP1CC","PRCC","PRKCD","PRKCE","PKN2","MAPK1","MAPK3","MAPK4","MAPK6",
                "MAPK7","MAPK13","PRNP","KLK10","RALA","RALB","RBBP4","RBL1","RBL2",
                "RFPL1","RGS2","RPS6KA1","RPS6KA3","MAPK12","SIAH1","SIAH2","SMARCB1",
                "SRC","TRIM21","STK10","STK11","SUV39H1","TAF1","DYNLT3","TERF2",
                "TFDP1","TFDP2","TP53BP2","TP73","TSG101","ZNF16","PTP4A1","EVI5",
                "CDK2AP1","CHAF1B","USP9X","DYRK3","RAE1","CAMK1","CDC14A","RUVBL1",
                "CCNK","KAT2B","CDC123","RNF8","TM4SF5","USP2","ARHGEF2","CCPG1",
                "KIF20B","WTAP","RAB11FIP3","RASSF2","RB1CC1","TLK1","KLHL21","DCLRE1A",
                "DMTF1","CASP8AP2","CHAF1A","PSME3","MAEA","PAK4","IKZF1","HMG20B",
                "ANAPC10","RACK1","TXNIP","CGRRF1","CGREF1","USP39","ZMYND11","ENTR1",
                "BLCAP","TXNL4A","TLK2","CNTRL","TRIOBP","KATNA1","FAM107A","WDR6",
                "RASSF1","PMF1","TUSC2","OIP5","VASH1","CEP164","CEP131","KLHL18",
                "USP22","TTC28","SPECC1L","ITGB3BP","CDK20","CCNDBP1","CD2AP","RABGAP1",
                "GSPT2","EID1","ANAPC13","ANAPC15","AHCTF1","FAM32A","APPL1",
                "SPAG8","MTBP","APEX2","SGSM3","RGCC","MCTS1","BRD7","UHRF1",
                "ANAPC2","ANAPC4","CUZD1","DYNC1LI1","ING4","ANAPC5","YTHDF2",
                "SPOUT1","ZC3HC1","CINP","MAP3K20","PELO","MIS18A","SPIN2A",
                "PIMREG","TET2","ERCC6L","TXNL4B","BANP","FIGN","THAP1","RIF1",
                "APPL2","ARL8B","RCBTB1","FANCI","STEAP3","MIS18BP1","HJURP",
                "STRADB","RBM38","PRR5","PRPF40A","URGCP","CCAR1","RIOK2",
                "TXLNG","KMT2E","RCC2","KLHL9","PARD3","PCNP","BIRC6","HACE1",
                "KLHL42","MICAL3","MARK4","BRINP2","CCAR2","AVPI1","CHTF18",
                "KIF13A","TSPYL2","TENT4B","MRPL41","NUP37","FSD1","CDC73",
                "SUV39H2","BORA","DBF4B","CABLES2","CDCA3","RHNO1","MARVELD1",
                "HMCN1","RASSF4","USP44","PDCD2L","MCM8","BEX2","TBRG1","AJUBA",
                "LRRCC1","LMLN","KLHL13","GADD45GIP1","NEK9","CABLES1","STRADA",
                "ACTR8","CCDC124","UHRF2","PARD3B","ANAPC16","NEDD1","CCSAP",
                "ARL8A","LIN54","POC5","MPLKIP","CTCFL","KCTD11","SPC24","PYHIN1",
                "SIK1","PPP1R1C","DAB2IP","CDCA2","ZBTB49","THAP5","CENPV","SENP5",
                "STOX1","DDIAS","HEPACAM","CDC26","BOD1L2","NSMCE2","LIN9","BRINP3",
                "TMPRSS11A","DCDC1","MCIDAS","NUP43","USP17L2","SPDYC","TRNP1","NUPR2",
                "PIM3","EIF2AK4","SPIN2B","GMNC","KLLN") %>% list()
cell_division <- c("BIRC5","ARF6","CCND1","BUB1","BUB1B","CCNA2","CCNB1","CCND2",
                   "CCND3","CCNE1","CCNF","CCNG1","CCNG2","CCNT1","CCNT2","CDK1",
                   "CDC6","CDC20","CDC25A","CDC25B","CDC25C","CDC27","CDK2","CDK3",
                   "CDK4","CDK5","CDK6","CDK7","CENPC","CENPE","CENPF","CETN1",
                   "CETN2","CETN3","RCC1","CKS1B","CKS2","CLTA","CLTC","CSNK1A1",
                   "DCTN1","DYNC1H1","E4F1","ENSA","EPB41","EPB41L2","GNAI1","GNAI2",
                   "GNAI3","NR3C1","HELLS","HNRNPU","KIF2A","KIF11","KIFC1","LIG1",
                   "LIG3","LIG4","LLGL2","MAD2L1","MAP4","NEDD9","NEK1","NEK2","NEK3",
                   "NUMA1","CHMP1A","CDK14","PIK3C3","PPP1CA","PPP1CB","PPP1CC",
                   "PRKCE","PKN2","RAD21","RALA","RALB","RAN","RB1","RBBP8","RPS3",
                   "NEK4","AURKA","AURKC","SYCP1","TACC1","DYNLT3","DYNLT1","TERF1",
                   "TPR","TSG101","UBE2I","VRK1","WEE1","ZNF16","ZNF207","EVI5",
                   "TUBA1A","HMGA2","USP9X","SMC1A","CDC7","MAD1L1","DYRK3","RAE1",
                   "CDC14A","RUVBL1","TNKS","BECN1","CDC23","CCNK","CDC123","CDC16",
                   "CCNA1","TIMELESS","WASL","BRSK2","RNF8","PRC1","LATS1","SMC3",
                   "CCNB2","CCNE2","CTDP1","ZW10","BUB3","PTTG1","RECQL5","BCAR1",
                   "BABAM2","KIF20B","RAB11FIP3","KNTC1","CKAP5","IST1","KLHL21",
                   "NCAPD2","KIF14","DCLRE1A","SMC4","STAG1","MAEA","KATNB1","CCNO",
                   "TUBA1B","SYCP2","ANAPC10","NDC80","MAD2L2","TACC3","CIB1","ZNRD2",
                   "SMC2","USP16","SPAG5","RGS14","USP39","NUDC","STAG2","ARPP19",
                   "NEK6","ENTR1","TXNL4A","MAPRE2","KIF2C","CNTRL","UBE2C","TPPP",
                   "TRIOBP","KATNA1","ZWINT","CHEK2","PMF1","OIP5","CEP164","MAPRE1",
                   "MAPRE3","SIRT2","TPX2","PDS5B","WAPL","SPART","SPART","CLASP2",
                   "POGZ","SMC5","ANKLE2","ABRAXAS2","FBXL7","PDS5A","KLHL18","NCAPD3",
                   "TTC28","CLASP1","HAUS5","MAU2","SPECC1L","NCAPH","ITGB3BP","CDK20",
                   "CD2AP","ANAPC13","ANAPC15","AHCTF1","NSL1","FBXO5","LATS2","VPS4A",
                   "UBE2S","EML4","BABAM1","ANAPC2","GPSM2","SAC3D1","ANAPC4","SYCP3",
                   "CUZD1","PARD6A","DYNC1LI1","FZR1","ANAPC5","ANAPC7","SPOUT1",
                   "ANAPC11","ZC3HC1","CINP","PELO","MIS18A","PIMREG","INO80","MAP10",
                   "HAUS6","NDE1","ERCC6L","NSUN2","NCAPG2","SPDL1","HAUS4","TIPIN",
                   "NRDE2","ZWILCH","FIGN","HAUS2","CDCA8","ARL8B","MIS18BP1","HAUS7",
                   "PRPF40A","INTS13","CHFR","CENPJ","PPP2R2D","RCC2","KLHL9","TEX14",
                   "KNL1","CHMP1B","SPC25","BIRC6","KLHL42","MICAL3","USP37","MARK4",
                   "KIF13A","NCAPG","TENT4B","ANAPC1","MIS12","NUP37","BRCC3","FSD1",
                   "HAUS3","BORA","MCMBP","DSN1","ANKRD53","CENPT","CEP63","REEP4",
                   "FAM83D","CDT1","CABLES2","SEH1L","CDCA3","NUF2","SETDB2","HMCN1",
                   "USP44","PARD6G","PARD6B","KIF2B","PSRC1","TUBA1C","KLHL22","MASTL",
                   "ZFYVE19","CCNB3","LRRCC1","LMLN","KLHL13","KNSTRN","BOD1","ZNF830",
                   "NEK9","CABLES1","HAUS8","SYCE1","ACTR8","CDCA5","CCDC124","HAUS1",
                   "PARD3B","ANAPC16","NEDD1","MISP","PLK5","CCSAP","ARL8A","DIS3L2",
                   "MPLKIP","KIF18B","SPC24","PHF13","PPP1R1C","SGO2","SGO1","SPICE1",
                   "CDCA2","SDE2","CENPV","CENPX","TUBB","SENP5","STOX1","CCNY","SKA1",
                   "SKA1","REEP3","SKA3","CDC26","HEPACAM2","SYCE2","EML3","BOD1L2",
                   "NSMCE2","DCDC1","SKA2","NUP43","CENPS","CENPW","KMT5A","SYCE3") %>% list()
mitotic_nuclear_division <- c("USP16","SMPD3","KLHDC8B","MAD2L2","NUDC","RCC1","CDCA8",
                              "CDC20","CENPF","NCAPH","BUB1","NME6","SPICE1","TACC3") %>% list()
TRPV1_expression <- c("TRPV1") %>% list()

brainH15 <- AddModuleScore(object = brainH15, features = Apoptosis, name = 'Apoptosis', assay = "SCT") 
# brainH15 <- AddModuleScore(object = brainH15, features = RTK_RAS_PI3K_pathway, name = 'RTK_RAS_PI3K_pathway', assay = "SCT") 
brainH15 <- AddModuleScore(object = brainH15, features = cell_cycle, name = 'cell_cycle', assay = "SCT") 
brainH15 <- AddModuleScore(object = brainH15, features = cell_division, name = 'cell_division', assay = "SCT") 
brainH15 <- AddModuleScore(object = brainH15, features = mitotic_nuclear_division, name = 'mitotic_nuclear_division', assay = "SCT") 
brainH15 <- AddModuleScore(object = brainH15, features = TRPV1_expression, name = 'TRPV1_expression', assay = "SCT")

library(ggplot2)
mydata <- brainH15@meta.data[, c("Glutamatergic_synaptic1", "TRPV1_expression1")]
dim(mydata)
mydata_filtered <- mydata %>%
  filter(Glutamatergic_synaptic1 > 0 & TRPV1_expression1 > 0)
# 查看过滤后的数据量
dim(mydata_filtered)  # 对比过滤前后的行数变化

ggplot(data = mydata_filtered, aes(Glutamatergic_synaptic1,TRPV1_expression1)) +
  geom_point(fill="black",colour="black",size=3,shape=21)+ # 绘制二维散点
  geom_smooth(method = "loess",span=0.4,se=TRUE,colour = "#00A5FF",fill="#00A5FF",alpha=0.2)
# 使用LOESS方法平滑数据，添加平滑曲线

# 加载必需的包
library(ggplot2)
library(ggpmisc)  # 提供stat_poly_eq函数

# 绘制散点图+loess拟合，并添加方程
ggplot(data = mydata_filtered, aes(x = Glutamatergic_synaptic1, y = TRPV1_expression1)) +  # 明确x和y的列名
  geom_point(fill = "black", colour = "black", size = 3, shape = 21) +  # 散点图
  geom_smooth(
    method = 'loess', 
    span = 0.4, 
    se = TRUE, 
    colour = "#00A5FF", 
    fill = "#00A5FF", 
    alpha = 0.2
  ) +  # loess拟合曲线
  # 添加loess拟合方程（关键代码）
  stat_poly_eq(
    formula = y ~ x,  # 公式（loess本质是局部多项式拟合，此处用y~x表示）
    method = "loess",  # 指定方法为loess
    span = 0.4,  # 与geom_smooth的span保持一致
    aes(label = after_stat(eq.label)),  # 显示方程
    parse = TRUE,  # 解析方程为数学表达式
    label.x.npc = 0.1,  # 方程在x轴的相对位置（0.1=左偏10%）
    label.y.npc = 0.9   # 方程在y轴的相对位置（0.9=上偏10%）
  ) +
  theme_bw()  # 美化主题

# (1)显示 loess 拟合的 R²; loess（局部加权回归）

library(ggplot2)
library(ggpubr)  # 提供stat_cor
library(ggpmisc)

p1 <- ggplot(data = mydata_filtered, 
             aes(x = Glutamatergic_synaptic1, y = TRPV1_expression1)) +
  geom_point(fill = "black", colour = "black", size = 1.8, shape = 21) +
  geom_smooth(
    method = 'loess', 
    span = 0.4, 
    se = TRUE, 
    colour = "#00A5FF", 
    fill = "#00A5FF", 
    alpha = 0.2
  ) +
  # 显示loess拟合的R²（关键：用stat_fit_glance提取）
  stat_fit_glance(
    method = "loess",
    span = 0.4,
    aes(label = sprintf("R² = %.3f", after_stat(r.squared))),  # 提取R²并格式化
    label.x.npc = 0.1,
    label.y.npc = 0.9,
    parse = FALSE  # 非方程形式，无需解析
  ) +
  # 可选：添加Pearson相关性（反映整体线性趋势）
  stat_cor(method = "pearson", label.x.npc = 0.1, label.y.npc = 0.8) +
  theme_bw()
p1

dir3 <- "D:/Desktop/GBM/250510/result/TRPV1/"
ggsave(paste(dir3,"SH15_Correlation_TRPV1.pdf", sep=""), p1, width = 4.5, height = 4)


(2)改用线性拟合（可稳定显示方程）
library(ggplot2)
library(ggpmisc)

ggplot(data = mydata_filtered, 
       aes(x = Glutamatergic_synaptic1, y = TRPV1_expression1)) +
  geom_point(fill = "black", colour = "black", size = 3, shape = 21) +
  geom_smooth(
    method = 'lm',  # 线性拟合
    se = TRUE, 
    colour = "#00A5FF", 
    fill = "#00A5FF", 
    alpha = 0.2
  ) +
  # 显示线性方程和R²
  stat_poly_eq(
    formula = y ~ x,
    method = "lm",  # 对应线性拟合
    aes(label = paste0(after_stat(eq.label), ", ", after_stat(rr.label))),  # 方程+R²
    parse = TRUE,
    label.x.npc = 0.1,
    label.y.npc = 0.9
  ) +
  theme_bw()

# 加载必要的库
library(dplyr)
library(corrplot)
head(brainH15@meta.data)

# 提取需要分析的列
cor_data1 <- brainH15@meta.data %>%
  select(
    # 突触相关评分
    Glutamatergic_synaptic1,
    GABAergic_synaptic1,
    Cholinergic_synaptic1
  )
cor_data2 <- brainH15@meta.data %>%
  select(
    # 细胞周期及通路相关指标
    Proliferation1,
    Apoptosis1,
    cell_cycle1,
    cell_division1,
    mitotic_nuclear_division1,
    TRPV1_expression1
  )

library(psych)
library(corrplot)
#计算基因表达量之间的pearson相关性；
ct1 <- corr.test(cor_data1, cor_data2, method = "pearson")
#提取相关性系数矩阵；
r1 <- ct1$r
#提取pvalue值矩阵；
p1 <- round(ct1$p,3)
#预览转置后的相关性系数矩阵和pvalue矩阵；
r2 <- t(r1)
p2 <- t(p1)
#使用显著性星号标记进行替换；
p2[p2>=0 & p2 < 0.001] <- "***"
p2[p2>=0.001 & p2 < 0.01] <- "**"
p2[p2>=0.01 & p2 < 0.05] <- "*"
p2[p2>=0.05 & p2 <= 1] <- ""

#载入pheatmap包；
library(pheatmap)
library(circlize)
#自定义颜色；
mycol<-colorRampPalette(c("#0f86a9", "white", "#ed8b10"))(200)
# col_fun1 = colorRamp2(c(-0.808, 0, 0.883), c("#0f86a9", "white", "#FC8452"))
# col_fun2 = colorRamp2(c(-0.808, 0, 0.883), c("#A5CC26", "white", "#FF7BAC"))
# col_fun3 = colorRamp2(c(-0.808, 0, 0.883), c("#3FA9F5", "white", "#FF931E"))
# col_fun4 = colorRamp2(c(-1, 0, 1), c("#ffa500", "white", "#B3A9EB"))
#绘制热图；
plot1 <- pheatmap(r2,scale = "none",
                  border_color ="white",
                  number_color="white",
                  fontsize_number=14,
                  fontsize_row=8,
                  fontsize_col=9,
                  cellwidth=15,
                  cellheight=15,
                  cluster_rows=T,
                  cluster_cols=T,
                  color = mycol,
                  display_numbers = p2,
                  show_rownames=T)

# 安装并加载包
library(RColorBrewer)

# 查看RdBu配色卡（包含不同梯度的选项）
display.brewer.pal(n = 11, name = "RdBu")  # n为颜色数量（3-11之间，奇数更对称）

# 提取配色（例如提取7种颜色）
my_colors <- rev(brewer.pal(n = 7, name = "RdBu"))
my_colors  # 输出：7种颜色的十六进制代码
plot1 <- pheatmap(r2,scale = "row",
                  # border_color ="white",
                  # number_color="white",
                  fontsize_number=14,
                  fontsize_row=8,
                  fontsize_col=9,
                  cellwidth=15,
                  cellheight=15,
                  cluster_rows=T,
                  cluster_cols=T,
                  color = my_colors,
                  # display_numbers = p2,
                  show_rownames=T)

dir3 <- "D:/Desktop/GBM/250510/result/Spatial_gaba_correaltion/"
ggsave(paste(dir3,"SH15_Correlation_Synpase_Spatial.pdf", sep=""), plot1, width = 4, height = 4)




