library(ggplot2)
setwd("E:/CGGA")
plotdata <- read.csv(file = "CGGA_grade_OSR2去NA和异常大值.csv")
ggplot(plotdata, aes(x = `LOGOSR2.1`, fill = Grade)) +
  # 绘制半透明的直方图
  geom_histogram(aes(y = ..density..), 
                 alpha = 0.6, 
                 position = "identity",
                 bins = 30,
                 color = "white") +
  # 添加光滑的密度曲线
  geom_density(alpha = 0.5, size = 0.8) +
  # 使用经典主题
  theme_classic() +
  # 自定义颜色
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
  # 添加标签和标题
  labs(title = "基因A表达水平在不同疾病类型中的分布",
       x = "基因A表达水平",
       y = "sample count",
       fill = "疾病类型") +
  # 调整图例位置
  theme(legend.position = c(0.9, 0.8),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))








ggplot(plotdata, aes(x = LOGOSR2.1, fill = Grade)) +
  geom_histogram(aes(y = ..count..),
                 position = "identity",
                 bins = 30,
                 alpha = 0.6,
                 color = "white") +
  geom_line(stat = "density", aes(y = ..density.. * (n/3) * 2.5, color = Grade), 
            size = 1, adjust = 1.5) + # 调整缩放因子使曲线与直方图匹配
  facet_wrap(~ Grade, ncol = 1) +
  labs(x = "基因A表达水平", 
       y = "样本计数", 
       title = "三种疾病类型的基因A表达水平分布（频率显示）") +
  theme_minimal()
