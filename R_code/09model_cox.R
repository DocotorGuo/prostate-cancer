rm(list=ls())

library(survival)
library(tidyverse)

# Univariate Cox analysis
rt=read.table("03Cox/TCGA_PRAD_diff_MGs_exprSetTime.txt",header=T,check.names=F,row.names=1,sep = "\t") 

library(survival)
rt1 <- rt

(variable_names<-colnames(rt)[-c(1:7)])
sur<-Surv(time=rt$Disease_Free_Survival, event = rt$Disease_Free_Status)

pFilter=0.05
outTab=data.frame()

sigGenes = colnames(rt1)[1:7]
for(gene in variable_names){
  if(sd(rt[,gene])<0.001){next}
  if(grepl("-", gene)){next}
  cox=coxph(as.formula(paste0('sur~',gene))  ,data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    group=ifelse(rt[,gene]>median(rt[,gene]),"high","low")
    diff=survdiff(sur ~group,data = rt)
    pValue=1-pchisq(diff$chisq,df=1)
    if(pValue<0.05){
      sigGenes=c(sigGenes,gene)
      outTab=rbind(outTab,
                   cbind(gene=gene,
                         KM=pValue,
                         HR=coxSummary$conf.int[,"exp(coef)"],
                         HR.95L=coxSummary$conf.int[,"lower .95"],
                         HR.95H=coxSummary$conf.int[,"upper .95"],
                         coxPvalue=coxP) )
    }
  }
}


write.table(outTab,file="03Cox/trianSet_uniCox.txt",sep="\t",row.names=F,quote=F)   

surSigExp=rt1[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)

write.table(surSigExp,file="03Cox/trianSet_uniSigExp.txt",sep="\t",row.names=F,quote=F)

#Draw the forest map function
bioForest=function(coxFile=null,forestFile=null){
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <-  gsub("(.*?)\\|(.*?)\\|(.*?)","\\3\\", rownames(rt)) 
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$coxPvalue<0.001, "<0.001", sprintf("%.3f", rt$coxPvalue))

  pdf(file=forestFile, width = 5,height = 5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(5,2))
  
  xlim = c(0,3)
  par(mar=c(4,3,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.6-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.6-0.5*0.2,n+1,'Pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard Ratio',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,0,2,1),mgp=c(2,0.5,0))
  xlim = c(0 ,max(as.numeric(hrLow),as.numeric(hrHigh)+0.5))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard Ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=1.5)
  abline(v=1,col="black",lty=2,lwd=1.3)
  boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.8)
  axis(1)
  dev.off()
}


bioForest(coxFile="03Cox/trianSet_uniCox.txt",forestFile="03Cox/trian_uniForest.pdf")

####################################
# Multivariate Cox +LASSO penalty  
options(stringAsFactors=F)
library(survival)                                         #引用包
library(glmnet)

rt=read.table("03Cox/trianSet_uniSigExp.txt",header=T,sep="\t",check.names=F,row.names = 1)    #读取输入文件
rt1 = rt
str(rt)
set.seed(9)

x = as.matrix(rt[,-c(1:7)])
y = Surv(time=rt$Disease_Free_Survival, event = rt$Disease_Free_Status)
cvfit = cv.glmnet(x,y, family = "cox", nfold = 10)
pdf("03Cox/cvfit.pdf",width = 4.6,height = 5)
plot(cvfit)
# abline(v = c(log(cvfit$lambda.min), log(cvfit$lambda.1se)),lty=2)+
text(x = log(cvfit$lambda.min),y = 12.6,
     paste('Lambda.min\n',round(cvfit$lambda.min,4)),cex=1,adj=0.9)
text(x = log(cvfit$lambda.1se)+0.01,y = 12.7,
     paste('Lambda.lse\n',round(cvfit$lambda.1se,4)),cex=1,adj=0.9)
dev.off()

fit <- glmnet(x, y, family = "cox",nfold = 10)
# pdf(file="03Cox/fit.pdf",width=4.6,height=5)
# plot(fit, xvar = "lambda", label = TRUE)
# dev.off()
# Draw nice Lasso diagrams  
library("reshape")
library("ggsci")
x <- coef(fit)  
tmp <- as.data.frame(as.matrix(x))
tmp$coef <- row.names(tmp)
tmp <- reshape::melt(tmp, id = "coef")
tmp$variable <- as.numeric(gsub("s", "", tmp$variable))
tmp$coef <- gsub('_','-',tmp$coef)
tmp$lambda <- fit$lambda[tmp$variable+1] # extract the lambda values
tmp$norm <- apply(abs(x[-1,]), 2, sum)[tmp$variable+1] # compute L1 norm  
tmp$coef2 <- ifelse(tmp$norm==max(tmp$norm),tmp$coef,NA)

ggplot(tmp,aes(log(lambda),value,color = coef)) + 
  geom_vline(xintercept = log(cvfit$lambda.min),size=0.8,color='grey60',alpha=0.8,linetype=2)+
  geom_line(size=1) + 
  xlab("Log Lambda") + 
  #xlab("L1 norm")+
  ylab('Coefficients')+
  theme_bw()+ 
  scale_color_manual(values = c(pal_npg()(10),pal_d3()(10),pal_jco()(7)))+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0.01,0.01))+
  theme(panel.grid = element_blank(),
        # axis.title = element_text(size=15,color='black'),
        # axis.text = element_text(size=12,color='black'),
        legend.title = element_blank(),
        legend.text = element_text(size=12,color='black'),
        legend.position = 'right')+
  annotate('text',x = -4.3,y=0.05,label='lambda.min =  0.0219',color='black')+
  guides(col=guide_legend(ncol = 1),shape = guide_legend(override.aes = list(size = 0.5)))

ggsave("03Cox/fit.pdf",width = 5.5,height = 4)

myCoefs <- coef(cvfit, s=cvfit$lambda.min)
lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
lasso_fea
sigGenes= colnames(rt1)[1:7]
rt <- rt1[,c(sigGenes,lasso_fea)]
rt$riskScore <-  apply(rt[,lasso_fea], 1, function(x) {x %*% myCoefs@x})

rt$risk=as.vector(ifelse(rt$riskScore>median(rt$riskScore),"High Risk","Low Risk"))

write.table(cbind(id=rownames(rt),rt),
            file="03Cox/trianSet_riskScore.txt",sep="\t",quote=F,row.names=F)

lassoGene <-  cbind(Gene = lasso_fea,Coef = myCoefs[which(myCoefs != 0 )])
# lassoGene <-  cbind(Gene = lasso_fea,Coef = myCoefs[1:3])
write.table(lassoGene,file="03Cox/trianSet_lasso_coef.txt",row.names = F,quote = F,sep = "\t")


#GEO external data set validation
GEOSet <- read.table("01rawData/GSE70768_symbol_exprSet_Time.txt",sep = "\t",header = T,check.names = F,row.names = 1)
GEOSet[1:5,1:7]
geoSet <- GEOSet[,c("Disease_Free_Survival","Disease_Free_Status" ,lasso_fea)]
geoSet$riskScore <- apply(geoSet[,lasso_fea], 1, function(x) {x %*% myCoefs@x})

geoSet$risk=as.vector(ifelse(geoSet$riskScore>median(geoSet$riskScore),"High Risk","Low Risk"))

write.table( cbind(id=rownames(geoSet),geoSet),
             file="03Cox/GSE70768_riskScore.txt",
             sep="\t",row.names = F,
             quote=F)

