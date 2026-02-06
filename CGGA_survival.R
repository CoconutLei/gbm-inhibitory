## load R package
library(survminer)
library(survival)
library(ggpubr)
library(gridExtra)
library(dplyr)
## import data
setwd("E:/CGGA")
dat<-read.table("CGGA.mRNAseq_693_clinical.20200506.txt", sep='\t', head=T,row.names = 1)

t_gsva <- t(gsva_data) #转置ssGSEA矩阵
df_gsva <- as.data.frame(t_gsva) #ssGSEA矩阵转数据框
dat <- cbind(dat, df_gsva)       #生存期合并ssGSEA
# dat <- cbind(dat,as.data.frame(t(BMP2)),as.data.frame(t(OSR2)))

dat <- dat %>% 
  rename(Expression = EN, PRS.type = PRS_type, Censor = Censor..alive.0..dead.1. )      #对应基因集改列名

#head(dat)
mat<-dat[!is.na(dat$PRS.type)&(dat$PRS.type=="Primary"|dat$PRS.type=="Recurrent")&
           !is.na(dat$Grade)&
           !is.na(dat$OS)&
           !is.na(dat$Censor),]

### all grade - Primary ###
matt<-mat[mat$PRS.type=='Primary',]
med.exp<-median(matt$Expression)
more.med.exp.index<-which(matt$Expression>=med.exp)
less.med.exp.index<-which(matt$Expression< med.exp)
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(OS, Censor) ~ status, data = matt)
s.diff<-survdiff(Surv(OS, Censor) ~ status, data = matt)

sdata.plot1<-ggsurvplot(s.fit,
                        data=matt,
                        pval = TRUE,
                        #conf.int = TRUE,
                        xlab = 'Time (day)',
                        ggtheme = theme_light(),
                        surv.median.line = 'hv',
                        title=paste0("All WHO grade survival (primary glioma)"))

### all grade - Recurrent ###
matt<-mat[mat$PRS.type=='Recurrent',]
med.exp<-median(matt$Expression)
more.med.exp.index<-which(matt$Expression>=med.exp)
less.med.exp.index<-which(matt$Expression< med.exp)
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(OS, Censor) ~ status, data = matt)
s.diff<-survdiff(Surv(OS, Censor) ~ status, data = matt)


sdata.plot2<-ggsurvplot(s.fit,
                        data=matt,
                        pval = TRUE,
                        #conf.int = TRUE,
                        xlab = 'Time (day)',
                        ggtheme = theme_light(),
                        surv.median.line = 'hv',
                        title=paste0("All WHO grade survival (recurrent glioma)"))

### WHO grade II - Primary ###
matt<-mat[mat$PRS.type=='Primary'&mat$Grade=="WHO II",]
med.exp<-median(matt$Expression)
more.med.exp.index<-which(matt$Expression>=med.exp)
less.med.exp.index<-which(matt$Expression< med.exp)
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(OS, Censor) ~ status, data = matt)
s.diff<-survdiff(Surv(OS, Censor) ~ status, data = matt)

sdata.plot3<-ggsurvplot(s.fit,
                        data=matt,
                        pval = TRUE,
                        #conf.int = TRUE,
                        xlab = 'Time (day)',
                        ggtheme = theme_light(),
                        surv.median.line = 'hv',
                        title=paste0("WHO grade II survival (primary glioma)"))

### WHO grade II - Recurrent ###
matt<-mat[mat$PRS.type=='Recurrent'&mat$Grade=="WHO II",]
med.exp<-median(matt$Expression)
more.med.exp.index<-which(matt$Expression>=med.exp)
less.med.exp.index<-which(matt$Expression< med.exp)
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(OS, Censor) ~ status, data = matt)
s.diff<-survdiff(Surv(OS, Censor) ~ status, data = matt)

sdata.plot4<-ggsurvplot(s.fit,
                        data=matt,
                        pval = TRUE,
                        #conf.int = TRUE,
                        xlab = 'Time (day)',
                        ggtheme = theme_light(),
                        surv.median.line = 'hv',
                        title=paste0("WHO grade II survival (recurrent glioma)"))

### WHO grade III - Primary ###
matt<-mat[mat$PRS.type=='Primary'&mat$Grade=="WHO III",]
med.exp<-median(matt$Expression)
more.med.exp.index<-which(matt$Expression>=med.exp)
less.med.exp.index<-which(matt$Expression< med.exp)
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(OS, Censor) ~ status, data = matt)
s.diff<-survdiff(Surv(OS, Censor) ~ status, data = matt)

sdata.plot5<-ggsurvplot(s.fit,
                        data=matt,
                        pval = TRUE,
                        #conf.int = TRUE,
                        xlab = 'Time (day)',
                        ggtheme = theme_light(),
                        surv.median.line = 'hv',
                        title=paste0("WHO grade III survival (primary glioma)"))

### WHO grade III - Recurrent ###
matt<-mat[mat$PRS.type=='Recurrent'&mat$Grade=="WHO III",]
med.exp<-median(matt$Expression)
more.med.exp.index<-which(matt$Expression>=med.exp)
less.med.exp.index<-which(matt$Expression< med.exp)
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(OS, Censor) ~ status, data = matt)
s.diff<-survdiff(Surv(OS, Censor) ~ status, data = matt)

sdata.plot6<-ggsurvplot(s.fit,
                        data=matt,
                        pval = TRUE,
                        #conf.int = TRUE,
                        xlab = 'Time (day)',
                        ggtheme = theme_light(),
                        surv.median.line = 'hv',
                        title=paste0("WHO grade III survival (recurrent glioma)"))

### WHO grade II - Primary ###
matt<-mat[mat$PRS.type=='Primary'&mat$Grade=="WHO IV",]
med.exp<-median(matt$Expression)
more.med.exp.index<-which(matt$Expression>=med.exp)
less.med.exp.index<-which(matt$Expression< med.exp)
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(OS, Censor) ~ status, data = matt)
s.diff<-survdiff(Surv(OS, Censor) ~ status, data = matt)

sdata.plot7<-ggsurvplot(s.fit,
                        data=matt,
                        pval = TRUE,
                        #conf.int = TRUE,
                        xlab = 'Time (day)',
                        ggtheme = theme_light(),
                        surv.median.line = 'hv',
                        title=paste0("WHO grade IV survival (primary glioma)"))

### WHO grade II - Recurrent ###
matt<-mat[mat$PRS.type=='Recurrent'&mat$Grade=="WHO IV",]
med.exp<-median(matt$Expression)
more.med.exp.index<-which(matt$Expression>=med.exp)
less.med.exp.index<-which(matt$Expression< med.exp)
matt$status<-NA
matt$status[more.med.exp.index]<-paste0('High (',length(more.med.exp.index),')')
matt$status[less.med.exp.index]<-paste0('Low (',length(less.med.exp.index),')')

s.fit<-survfit(Surv(OS, Censor) ~ status, data = matt)
s.diff<-survdiff(Surv(OS, Censor) ~ status, data = matt)

sdata.plot8<-ggsurvplot(s.fit,
                        data=matt,
                        pval = TRUE,
                        #conf.int = TRUE,
                        xlab = 'Time (day)',
                        ggtheme = theme_light(),
                        surv.median.line = 'hv',
                        title=paste0("WHO grade IV survival (recurrent glioma)"))

## output pdf 
splots <- list()
splots[[1]]<-sdata.plot1
splots[[5]]<-sdata.plot2
splots[[2]]<-sdata.plot3
splots[[6]]<-sdata.plot4
splots[[3]]<-sdata.plot5
splots[[7]]<-sdata.plot6
splots[[4]]<-sdata.plot7
splots[[8]]<-sdata.plot8
arrange_ggsurvplots(splots,nrow=4,ncol=2)


