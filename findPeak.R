#读取文件 每个细胞一列
data <- read.csv("E:/IN/fig3 data/calcium imaging/result/coactive/mono/coactiveB.csv")
frametime = 2 #每帧多少秒

#找到所有极大值点
peak_1 <- 0*data
peak_2 <- peak_1
cellno = length(data[1,])
frameno = length(data[,1])
currentcell = 0
currentframe = 0
repeat{
  currentcell = currentcell + 1
  currentframe = 0
  repeat{
    currentframe = currentframe + 1
    if((data[currentframe,currentcell]) <= (data[(currentframe+1),currentcell])){peak_1[currentframe,currentcell] <- 1}
    if((data[currentframe,currentcell]) > (data[(currentframe+1),currentcell])){peak_1[currentframe,currentcell] <- 0}
    if(currentframe >= (frameno-1)){break}
  }
  
  
  if(currentcell >= cellno){break}
}
# 筛出变化幅度大于20%
result <- peak_2[c(1,2),]
row.names(result) <- c("peak count","mHz")
currentcell = 0
currentframe = 0
repeat{
  currentcell = currentcell + 1
  currentframe = 0
  repeat{
    currentframe = currentframe + 1
    if(identical((peak_1[c(currentframe,(currentframe + 1)),currentcell]),(c(1,0))) & (data[(currentframe+1),currentcell]>1.2) ){peak_2[currentframe,currentcell] <- 1}
    if((identical((peak_1[c(currentframe,(currentframe + 1)),currentcell]),(c(1,0))) & (data[(currentframe+1),currentcell]>1.2)) == FALSE){peak_2[currentframe,currentcell] <- 0}
    if(currentframe >= (frameno-1)){
      result[1,currentcell] <- sum(peak_2[,currentcell])
      result[2,currentcell] = (1000*result[1,currentcell])/(frametime*frameno)
      break}
  }
  
  
  if(currentcell >= cellno){break}
}



