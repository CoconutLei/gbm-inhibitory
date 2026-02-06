
#定义模长公式modulus
modulus <- function(vector){
  return(sqrt(sum(vector^2)))
}

#定义相似值计算公式
Sindex <- function(Ca,Cb){
  return(2*sum((Ca*Cb))/(((modulus(Ca))^2)+(((modulus(Cb))^2))))
}


#bootstrapping procedure封装成函数
Bootstrapping <- function(Ca){
  framelength <- length(Ca)
  randomshift <- sample(1:(framelength-1),1)
  Cashift <- Ca
  Cashift[c(1:randomshift)] <-Ca[c((framelength-randomshift+1):framelength)]
  Cashift[c((randomshift+1):framelength)] <- Ca[c(1:(framelength-randomshift))]
  return(Cashift)
}

# cell pair sindex 分布

Sindexdistribute <- function(Cad,Cbd){
  repeattime = 1000
  currentrepeat = 0
  Sindexdistribution <- c(1:repeattime)
  repeat{
    currentrepeat = currentrepeat+1
    Sindexdistribution[currentrepeat] <- Sindex((Bootstrapping(Cad)),Cbd)
    if(currentrepeat >= repeattime) {
      break
    }
  }
  return(Sindexdistribution)
}

# 检测是否显著
Similarindexsignificance <- function(cella,cellb){
  presentage <- quantile(Sindexdistribute(cella,cellb),probs = Sindex(cella,cellb))
  if(presentage > 0.99){
    return(c("TRUE",unname(presentage),names(presentage)))
  }
  if(presentage <= 0.99){
    return(c("FALSE",unname(presentage),names(presentage)))
  }
}
#放缩函数
scale_values <- function (x){10*(x-min(x))/(max(x)-min(x))}

#读取文件 每个细胞一列 第一行是细胞编号
setwd("E:/IN/fig3 data/calcium imaging/result/coactive/MGE")
data <- read.csv("rmact.csv")
data <- scale_values(data)
result_sig <- matrix(numeric((length(data[1,]))^2),nrow = length(data[1,]),ncol = length(data[1,]),dimnames = list(colnames(data),colnames(data)))
result_presentage <- result_sig
result_Sindexvalue <- result_sig

repeats = length(data[1,])
cellthatshift = 0
cellthatnotshift = 0
repeat{
  
  cellthatshift = cellthatshift + 1
  cellthatnotshift = cellthatshift + 1
  repeat{
    
    
    result_all <- Similarindexsignificance(data[,cellthatshift],data[,cellthatnotshift])
    result_sig[cellthatshift,cellthatnotshift] <- result_all[1]
    result_Sindexvalue[cellthatshift,cellthatnotshift] <- result_all[2]
    result_presentage[cellthatshift,cellthatnotshift] <- result_all[3]
    cellthatnotshift = cellthatnotshift + 1
    if(cellthatnotshift > repeats) {
      break
    }
    
    
  }
  
  
  if(cellthatshift >= (repeats-1)) {
    break
  }
}                    
write.csv(result_presentage,file = "result_presentageRM.csv" )
write.csv(result_sig,file = "result_sigRM.csv" )
write.csv(result_Sindexvalue,file = "result_SindexvalueRM.csv" )
