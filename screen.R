library(GSVA)
library(dplyr)
library(ggplot2)
setwd("E:/CGGA")
expr <- read.table("CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)   
expr<- as.matrix(expr)
dat<-read.table("CGGA.mRNAseq_693_clinical.20200506.txt", sep='\t', head=T,row.names = 1)
datall <- cbind(dat,as.data.frame(t(expr)))

setwd("E:/IN/fig3 data/筛选")
df <- read.csv("res_df.csv")
df <- na.omit(df)
df[,8] <- rank(df[,2])
colnames(df)[8] <- "FCrank"
df[,9] <- (max(df[,8])-df[,8]+1)
colnames(df)[9] <- "FCrevrank"

maligdf <- read.csv("maligdf.csv",header = T,row.names = 1)
surv <- read.csv("GBM_gene_log_rank_results.csv", header = T)
inferiorsur <- surv[surv$survival_effect=="reduce survival",1]
maligdf2 <- inferiorsur[inferiorsur %in% toupper(df[,7])]

calref <- read.csv("calcium related genes.csv",row.names = 1)
calrefmarker <- list("calref1" = calref[,1])
params <- ssgseaParam(
  exprData = expr,              # 强制转换为矩阵格式
  geneSets = calrefmarker,        # 基因集列表（必须为list格式）
  alpha = 0.25,                 # 默认权重参数
  normalize = TRUE,             # 启用结果标准化（推荐）
)
calrefgsva_data <- gsva(params, verbose = TRUE)
datall <- cbind(as.data.frame(t(calrefgsva_data)),datall)
dat4 <- datall[datall$Grade=="WHO IV",]
maligdf1 <- toupper(maligdf[,7])

result <- matrix(0:0,nrow = 177,ncol=3)
a = 1
while (a<178) {
  result[a,1] <- maligdf1[a]
  result[a,2] <- as.numeric((cor.test(as.numeric(dat4[,maligdf1[a]]),as.numeric(dat4[,"calref1"]),method=c("pearson"))[3]))
  result[a,3] <- as.numeric((cor.test(as.numeric(dat4[,maligdf1[a]]),as.numeric(dat4[,"calref1"]),method=c("pearson"))[4]))
  a=a+1
}
colnames(result) <- c("symbol","P-value","R-value")

candidates1 <- result[result[,2] <= 0.05,1]
candidates2 <- result[result[,2] >=1,1]
candidates <- c(candidates1,candidates2)
df <- df %>%
  mutate(IsCandidate = ifelse(toupper(SYMBOL) %in% candidates, "Candidate", "Other"))

result2 <- matrix(0:0,nrow = 1504,ncol=3)
a = 1
while (a<1505) {
  result2[a,1] <- maligdf2[a]
  result2[a,2] <- as.numeric((cor.test(as.numeric(dat4[,maligdf2[a]]),as.numeric(dat4[,"calref1"]),method=c("pearson"))[3]))
  result2[a,3] <- as.numeric((cor.test(as.numeric(dat4[,maligdf2[a]]),as.numeric(dat4[,"calref1"]),method=c("pearson"))[4]))
  a=a+1
}
colnames(result2) <- c("symbol","P-value","R-value")

candidates21 <- result2[result2[,2] <= 0.05,1]
candidates22 <- result2[result2[,2] >=1,1]
candidates222 <- c(candidates21,candidates22)

df <- df %>%
  mutate(IsCandidate = ifelse(toupper(SYMBOL) %in% candidates222, "Candidate", "Other"))

p <- ggplot(df, aes(x = FCrevrank, y = log2FoldChange, color = IsCandidate)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c("Candidate" = "red", "Other" = "gray60")) +
  labs(
    title = "Differential Gene Expression",
    x = "FC Rank",
    y = "log2(Fold Change)",
    color = "Gene Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )
print(p)
candcsv <- result2[(result2[,1] %in% candidates222),]
candcsvfinal <- merge(candcsv,surv,by = "gene")