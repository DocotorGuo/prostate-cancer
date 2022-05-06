# Risk and clinical correlation analysis
library(survival)
library(survminer)
mycol <- c('#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")
barplot(1:10,col =mycol )

rt=read.table("09Independent/trianSet_riskScore.txt",sep="\t",header=T,row.names=1,check.names=F) 

rt$age <- ifelse(rt$age >55,"Age >=55","Age <55")
rt$age <- factor(rt$age ,levels = c("Age <55","Age >=55"))

my_comparisons <- list(c("Age <55", "Age >=55"))

p1 <- ggplot(rt,aes(age,riskScore,fill=age))+
  geom_boxplot(outlier.colour = NA,notch = T,size = 0.3)+
  geom_jitter(shape = 21,size=2,width = 0.2)+
  geom_violin(pDFSition = pDFSition_dodge(width = .75), 
              size = 0.3,alpha = 0.4,trim = T)+
  theme_classic()+
  theme(legend.pDFSition = 'none',
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title.x = element_blank())+
  scale_fill_manual(values=mycol[c(1,2)])+
  stat_compare_means(comparisons = my_comparisons,
                     label = 'p.signif')
  # stat_compare_means(label.y = max(plotData$StromalScore)+5.5)
ggsave("09Independent/age_violin.pdf",width = 4,height = 4)

rt=read.table("09Independent/trianSet_riskScore.txt",sep="\t",header=T,row.names=1,check.names=F) 
rt <- rt[!is.na(rt$pathologic_T ) ,]

rt$pathologic_T <- factor(rt$pathologic_T  ,ordered = T)

# my_comparisons <- list(c("Stage I-II", "Stage III-IV"))
my_comparisons <- split(t(combn(levels(rt$pathologic_T),2)),1:nrow(t(combn(levels(rt$pathologic_T),2))))

p2 <- ggplot(rt,aes(pathologic_T,riskScore,fill=pathologic_T))+
  geom_boxplot(outlier.colour = NA,notch = T,size = 0.3)+
  geom_jitter(shape = 21,size=2,width = 0.2)+
  geom_violin(pDFSition = pDFSition_dodge(width = .75), 
              size = 0.3,alpha = 0.4,trim = T)+
  theme_classic()+
  theme(legend.pDFSition = 'none',
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title.x = element_blank())+
  scale_fill_manual(values=mycol)+
  stat_compare_means(comparisons = my_comparisons,
                     label = 'p.signif')
ggsave("09Independent/pathologic_T_violin.pdf",width = 4,height = 4)

rt=read.table("09Independent/trianSet_riskScore.txt",sep="\t",header=T,row.names=1,check.names=F) 
rt <- rt[!is.na(rt$pathologic_N ) ,]

rt$pathologic_N <- factor(rt$pathologic_N)

my_comparisons <- split(t(combn(levels(rt$pathologic_N),2)),1:nrow(t(combn(levels(rt$pathologic_N),2))))

p3 <- ggplot(rt,aes(pathologic_N,riskScore,fill=pathologic_N))+
  geom_boxplot(outlier.colour = NA,notch = T,size = 0.3)+
  geom_jitter(shape = 21,size=2,width = 0.2)+
  geom_violin(pDFSition = pDFSition_dodge(width = .75), 
              size = 0.3,alpha = 0.4,trim = T)+
  theme_classic()+
  theme(legend.pDFSition = 'none',
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title.x = element_blank())+
  scale_fill_manual(values=mycol[c(1,2)])+
  stat_compare_means(comparisons = my_comparisons,
                     label = 'p.signif')
ggsave("09Independent/pathologic_N_violin.pdf",width = 4,height = 4)

rt=read.table("09Independent/trianSet_riskScore.txt",sep="\t",header=T,row.names=1,check.names=F) 
# rt <- rt[rt$pathologic_t != "TX",]
rt <- rt[!is.na(rt$psa_value ) ,]
rt$psa <- ifelse(rt$psa_value >=0.4,"PSA >=0.4","PSA <0.4")

rt$psa <- factor(rt$psa)
my_comparisons <- split(t(combn(levels(rt$psa),2)),1:nrow(t(combn(levels(rt$psa),2))))

p4 <- ggplot(rt,aes(psa,riskScore,fill=psa))+
  geom_boxplot(outlier.colour = NA,notch = T,size = 0.3)+
  geom_jitter(shape = 21,size=2,width = 0.2)+
  geom_violin(pDFSition = pDFSition_dodge(width = .75), 
              size = 0.3,alpha = 0.4,trim = T)+
  theme_classic()+
  theme(legend.pDFSition = 'none',
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title.x = element_blank())+
  scale_fill_manual(values=mycol)+
  stat_compare_means(comparisons = my_comparisons,
                     label = 'p.signif')
ggsave("09Independent/PSA_violin.pdf",width = 4,height = 4)


rt=read.table("09Independent/trianSet_riskScore.txt",sep="\t",header=T,row.names=1,check.names=F) 
# rt <- rt[rt$pathologic_m != "MX",]
rt <- rt[!is.na(rt$gleason_score ) ,]
rt$gleason_score <- ifelse(rt$gleason_score >=8,"gleason >=8","gleason <8")

rt$gleason_score <- factor(rt$gleason_score,ordered = T)

my_comparisons <- split(t(combn(levels(rt$gleason_score),2)),1:nrow(t(combn(levels(rt$gleason_score),2))))

p5 <- ggplot(rt,aes(gleason_score,riskScore,fill=gleason_score))+
  geom_boxplot(outlier.colour = NA,notch = T,size = 0.3)+
  geom_jitter(shape = 21,size=2,width = 0.2)+
  geom_violin(pDFSition = pDFSition_dodge(width = .75), 
              size = 0.3,alpha = 0.4,trim = T)+
  theme_classic()+
  theme(legend.pDFSition = 'none',
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title.x = element_blank())+
  scale_fill_manual(values=mycol)+
  stat_compare_means(comparisons = my_comparisons,
                     label = 'p.signif')
ggsave("09Independent/gleason_score_violin.pdf",width = 4,height = 4)



ggarrange(p1,p2,p3,p4,p5,labels = c("A", "B","C","D","E"),ncol = 3, nrow = 2)

ggsave("09Independent/all_clinical_violin.pdf",width = 11,height = 7)


# Clinical correlation heatmap analysis
rt=read.table("09Independent/trianSet_riskScore.txt",sep="\t",header=T,row.names=1,check.names=F) 
# rt$pathologic_stage <-  gsub("[ABCabc]","", rt$pathologic_stage)
rt <- na.omit(rt)

rt=rt[order(rt$risk,decreasing = T),]

rt$age <- ifelse(rt$age >=55,"Age >=55","Age <55")
rt$age <- factor(rt$age ,levels = c("Age >=55","Age <55"))

rt$pathologic_T <- factor(rt$pathologic_T,ordered = T)

rt$pathologic_N <- factor(rt$pathologic_N,ordered = T)
rt$PSA <- ifelse(rt$psa_value >=0.4,"PSA >=0.4","PSA <0.4")
rt$PSA <- factor(rt$PSA ,levels = c("PSA >=0.4","PSA <0.4"))

rt$gleason <- ifelse(rt$gleason_score >=8,"gleason >=8","gleason <8")
rt$gleason <- factor(rt$gleason ,levels = c("gleason >=8","gleason <8"))


rt1=rt[c(8:(ncol(rt)-4))]
rt1 =rt1[,c("RRM2","TK1","AK5","HAGHL","ISYNA1","ALDH1A2","CA14","ACSS3","PYGM","CD38")]
max(rt1) 
min(rt1)

rt2=t(rt1)
rownames(rt2) <- gsub("(.*?)\\|(.*?)\\|(.*?)","\\3\\",rownames(rt2))


annotation=data.frame(
  age = rt$age,
  pathologic_T = rt$pathologic_T,
  pathologic_N = rt$pathologic_N,
  PSA = rt$PSA,
  gleason = rt$gleason,
  
  riskScore = rt$risk)

rownames(annotation)=colnames(rt2)
ann_colors = list(
  age=c(`Age >=55` ="#f2994a",`Age <55`= "#f2c94c"),
  # genders = c(female="#FF9289",male="#96CA00"),
  # pathologic_stage=c(`Stge I`= "#cddc39",`Stge II`="#ffeb3b", `Stge III`="#ffc107", `Stge IV`="#ff9800"),
  `pathologic_T`= c(`T2` = "#ffaf7b",`T3` ="#d76d77" ,`T4`="#3a1c71"),
  `pathologic_n` = c(N0="#4caf50",`N1`="#00bcd4"),
  `PSA`  = c(`PSA >=0.4`="#764ba2",`PSA <0.4`="#667eea"),
  `gleason`  = c(`gleason >=8`="#fcb69f",`gleason <8`="#ffecd2"),
   riskScore=c(`Low Risk` = "#58a591",`High Risk` ="#f6685e"))

pdf(file="09Independent/clinical_heatmap.pdf",width = 8,height = 5)
pheatmap(rt2, 
         scale = "row",
         annotation_col=annotation, 
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         fontsize_row=11,
         show_colnames = F,
         fontsize_col=3,
         color = colorRampPalette(c("#028846", "white", "red"))(50) )
dev.off()



# Independent prognostic analysis
library(survival)
rt=read.table("09Independent/trianSet_riskScore.txt",header=T,check.names=F,row.names=1,sep = "\t")

rt <- na.omit(rt)
rt <- rt[,-(c(8:17,ncol(rt)))]
rt1 <- rt

rt$pathologic_T <- as.numeric(as.factor(rt$pathologic_T))
# rt$pathologic_m <- as.numeric(as.factor(rt$pathologic_m))
rt$pathologic_N <- as.numeric(as.factor(rt$pathologic_N))
colnames(rt)[6:7] <- c("PSA","gleason")

uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(Disease_Free_Survival, Disease_Free_Status) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(uniTab,file="09Independent/uniCox_independent.txt",sep="\t",row.names=F,quote=F)

#Factors with P value <0.05 were selected
(names<-uniTab$id[as.numeric(uniTab$pvalue)<0.05])

#Multivariate independent prognostic analysis
multiCox=coxph(Surv(Disease_Free_Survival, Disease_Free_Status) ~ ., data = rt[,c("Disease_Free_Survival", "Disease_Free_Status",names)])
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file="09Independent/multiCox_independent.txt",sep="\t",row.names=F,quote=F)



# Draw the forest map function
bioForest=function(coxFile=null,forestFile=null){
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

  pdf(file=forestFile, width = 6,height = 3.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(4,2))

  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.8-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.8-0.5*0.2,n+1,'Pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard Ratio',cex=text.cex,font=2,adj=1,)

  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)+1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard Ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=1.5)
  abline(v=1,col="black",lty=2,lwd=1.3)
  boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.8)
  axis(1)
  dev.off()
}


bioForest(coxFile="09Independent/uniCox_independent.txt",forestFile="09Independent/uniForest.pdf")


bioForest=function(coxFile=null,forestFile=null){
  #读取输入文件
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #输出图形
  pdf(file=forestFile, width = 6,height = 3)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(4,2))
  
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.8-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.8-0.5*0.2,n+1,'Pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard Ratio',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)+1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard Ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=1.5)
  abline(v=1,col="black",lty=2,lwd=1.3)
  boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=0.8)
  axis(1)
  dev.off()
}
bioForest(coxFile="09Independent/multiCox_independent.txt",forestFile="09Independent/multiForest.pdf")


##########  Line diagram and calibration curve
library(survival)
library(regplot)
library(rms)
library(nomogramEx)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
rt=read.table("09Independent/trianSet_riskScore.txt",sep="\t",header=T,row.names=1,check.names=F)  

rt <- na.omit(rt)
colnames(rt)[6:7] <- c("PSA","gleason")

dd <- datadist(rt)
options(datadist="dd")
options(na.action="na.delete")
summary(rt$Disease_Free_Survival)
coxpbc <- cph(formula = Surv(Disease_Free_Survival,Disease_Free_Status) ~ gleason+riskScore ,data=rt,x=T,y=T,surv = T,na.action=na.delete)  #,time.inc =2920

print(coxpbc)

surv <- Survival(coxpbc) 

surv1 <- function(x) surv(12,x)
surv3 <- function(x) surv(36,x)
surv5 <- function(x) surv(60,x)


x <- nomogram(coxpbc,fun = list(surv1,surv3,surv5),lp=T,
              funlabel = c('1-year survival Probability','2-year survival Probability','3-year survival Probability'),
              maxscale = 100,fun.at = c(0.95,0.9,0.8,0.7,0.5,0.3,0.1))

pdf("10nomogram/nomogram_classical.pdf",width = 9,height = 6)
plot(x, lplabel="Linear Predictor",
     xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.05,lwd=5,
     label.every = 1,col.grid = gray(c(0.8, 0.95)))
dev.off()

# Calibration Curve is drawn for verification
f1 <- cph(formula =  Surv(Disease_Free_Survival,Disease_Free_Status) ~ gleason+riskScore,data=rt,x=T,y=T,surv = T,na.action=na.delete,time.inc = 12) 
#m=50
cal1 <- calibrate(f1, cmethod="KM", method="boot",u=12,m=50,B=1000) 

f3 <- cph(formula =  Surv(Disease_Free_Survival,Disease_Free_Status) ~ gleason+riskScore ,data=rt,x=T,y=T,surv = T,na.action=na.delete,time.inc = 36) 
#m=50  
cal3 <- calibrate(f3, cmethod="KM", method="boot",u=36,m=50,B=1000) 

f5 <- cph(formula =  Surv(Disease_Free_Survival,Disease_Free_Status) ~ gleason+riskScore,data=rt,x=T,y=T,surv = T,na.action=na.delete,time.inc = 50) 
#m=50 
cal5 <- calibrate(f5, cmethod="KM", method="boot",u=50,m=50,B=1000) 

pdf("10nomogram/calibration_compare.pdf",width = 5,height = 5.5)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced DFS (%)",ylab = "Observed DFS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#F0027F"),
     xlim = c(0,1),ylim= c(0,1),col = c("#F0027F"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#F0027F"), pch = 16)

mtext("")

plot(cal5,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)


abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","3-year","5-year"), #图例文字
       col =c("#2166AC","#F0027F","#B2182B"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
dev.off()



pdf("10nomogram/calibration_1y.pdf",width = 5,height = 5.5)
plot(cal1,
     lwd = 2,#error bar的粗细
     lty = 1,#error bar的类型，可以是0-6
     errbar.col = c("#2166AC"),#error bar的颜色
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced DFS (%)",ylab = "Observed DFS (%)",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6) #字的大小
lines(cal1[,c('mean.predicted',"KM")], 
      type = 'b', #连线的类型，可以是"p","b","o"
      lwd = 2, #连线的粗细
      pch = 16, #点的形状，可以是0-20
      col = c("#2166AC")) #连线的颜色
mtext("")
box(lwd = 1) #边框粗细
abline(0,1,lty = 3, #对角线为虚线
       lwd = 2, #对角线的粗细
       col = c("#224444")#对角线的颜色
) 
dev.off()

pdf("10nomogram/calibration_3y.pdf",width = 5,height = 5.5)
plot(cal3,
     lwd = 2,
     lty = 1,
     errbar.col = c("#F0027F"),
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced DFS (%)",ylab = "Observed DFS (%)",
     col = c("#F0027F"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal3[,c('mean.predicted',"KM")],
      type= 'b',
      lwd = 2,
      col = c("#F0027F"),
      pch = 16)
mtext("")
box(lwd = 1)
abline(0,1,lty= 3,
       lwd = 2,
       col =c("#224444"))
dev.off()


pdf("10nomogram/calibration_5y.pdf",width = 5,height = 5.5)
plot(cal5,
     lwd = 2,
     lty = 1,
     errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced DFS (%)",ylab = "Observed DFS (%)",
     col = c("#B2182B"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal5[,c('mean.predicted',"KM")],
      type= 'b',
      lwd = 2,
      col = c("#B2182B"),
      pch = 16)
mtext("")
box(lwd = 1)
abline(0,1,lty= 3,
       lwd = 2,
       col =c("#224444"))
dev.off()


#Calculate cindex  
library(survival) 
library(survcomp) 

fit <- cph(formula = Surv(Disease_Free_Survival,Disease_Free_Status) ~ gleason+riskScore ,data=rt)  #,time.inc =2920
cindex <- concordance.index(predict(fit),surv.time = rt$Disease_Free_Survival,surv.event = rt$Disease_Free_Status,method = "noether")
cindex$c.index; cindex$lower; cindex$upper
# [1] 0.7560501
# [1] 0.7035654
# [1] 0.8085348

