rm(list=ls())
#survival analysis 
library(survival)
library(survminer)
inputfile<- "03Cox/GSE70768_riskScore.txt"
outputfile <- "GSE70768"
rt <- read.table(inputfile,header=T,sep="\t",check.names = F)

diff=survdiff(Surv(Disease_Free_Survival, Disease_Free_Status) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,length(diff$n) - 1)
HR = (diff$obs[2]/diff$exp[2])/(diff$obs[1]/diff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/diff$exp[2]+1/diff$exp[1]))


paste0(round(HR,3),"(" ,paste(round(low95,3), round(up95,3), sep = " - "),")" )

median(rt$Disease_Free_Survival)
pval <- ifelse(pValue < 0.001, "p < 0.001", 
               paste("p = ",round(pValue,3), sep = ""))
pval


fit <- survfit(Surv(Disease_Free_Survival, Disease_Free_Status) ~ risk, data = rt)
ggsurvplot(fit,
           pval = pval,
           linetype = "solid",  
           palette = c("#FF0033","#1B9E77"), #线的颜色对应高、低(自定义调色板) FF0033红色
           # surv.median.line = "hv",# 增加中位生存时间
           # conf.int = TRUE,# 增加置信区间
           title = outputfile,
           #ylab = "Progression-free survival (percentage)",
           xlab = " Time (Months)",
           # legend.title = "Survival Plot",
           legend = c(0.8,0.3),
           legend.labs = c("High Risk","Low Risk"),
           legend.title="",
           risk.table = T, # 添加风险表
           #surv.scale = "percent", # 生存曲线的比例转换。允许的值为"default"或"percent"。
           risk.table.title="",
           risk.table.height  = 0.25,
           #risk.table.fontsize = 0.2,
           surv.plot.height =0.75,# 生存图的高度，默认为0.75；
           ggtheme = theme_classic(),  #theme_bw(),
           # tables.theme = theme_cleantable()
)
dev.copy2pdf(file = paste0("08ROC/",outputfile,"_survival.pdf") , width = 4.3,height = 5)
dev.off()


library(survivalROC)
period_time <- 12

rt <- read.table(inputfile,header=T,sep="\t",check.names = F)

rocCol=c("#3cb346","#eeb401", "#ef1828","#942d8d")
aucText <- c()
pdf(file=paste0("08ROC/",outputfile,"_ROC.pdf"),width=4.5,height=4.5)
par(mar=c(4,4,2,1),mgp=c(2,0.5,0))
roc=survivalROC(Stime=rt$Disease_Free_Survival,
                status=rt$Disease_Free_Status,
                marker = rt$riskScore ,
                predict.time =period_time ,
                #span = 0.25*nrow(rt)^(-0.20),##span,NNE法的namda
                method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="#1B9E77",
     xlab="1-Specificity (False Positive)", ylab="Sensitivity (True Positive)",
     #main="ROC curve, Method = KM",
     main=outputfile,
     lwd = 1.5, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1,lty= 3)
j <- 0
for(i in c(1,3,5)){
        
        roc=survivalROC(Stime=rt$Disease_Free_Survival, status=rt$Disease_Free_Status, marker = rt$riskScore,
                        # span = 0.25*nrow(rt)^(-0.20),##span,NNE法的namda
                        predict.time =period_time*i, method="KM")
        
        j=j+1
        aucText=c(aucText,paste0(i," years"," (AUC=",sprintf("%.3f",roc$AUC),")"))
        lines(roc$FP, roc$TP, type="l",col=rocCol[j],xlim=c(0,1), ylim=c(0,1),lwd = 1.8, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
        
}

legend("bottomright", aucText,lty= 1,lwd=1.8,bty="n",col=rocCol)
dev.off()

## Training set risk curve

rt <- read.table(inputfile,header=T,sep="\t",check.names = F)
rt$Disease_Free_Status <- ifelse(rt$Disease_Free_Status==0,'NO Recurrence','Recurrence')
rt$risk <- factor(rt$risk,levels = c("Low Risk","High Risk"))
rt=rt[order(rt$riskScore),] 
rt$patient = 1:nrow(rt)

p1 <- ggplot(data=rt,mapping = aes(x=patient,y=riskScore))+
        geom_point(aes(color=risk),size=0.5)+
        theme_classic()+
        scale_color_manual(name="Risk Group",values = c("#028846","red"))+
        #x-axis
        theme(axis.line.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank())+
        ylab("Risk Score")+
        
        geom_vline(xintercept=sum(rt$risk=="Low Risk"),colour="black", linetype=2,size=0.5)+
        scale_x_continuous(expand = c(0,0))
p1 


p2 <- ggplot(data = rt,mapping = aes(x=patient,y=Disease_Free_Survival,col=Disease_Free_Status))+
        geom_point(aes(col=Disease_Free_Status),size=0.5)+ scale_color_manual(name="Status",values = c("#028846","red")) +
        theme_classic()+
        #x a-xis
        theme(axis.line.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.x=element_blank())+
        ylab("Survival time(Months)")+
       
        geom_vline(xintercept = sum(rt$risk=="High Risk") ,linetype=2,size=0.5)+
        scale_x_continuous(expand = c(0,0))
p2

middle = ggplot(rt, aes(
        x = patient,
        y = 1)) +
        geom_tile(aes(fill = risk))+theme_classic()+
        scale_fill_manual(name="Risk Group",  values = c("#028846","red"))+
        theme(
                legend.position = "none",
                axis.line = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                axis.title = element_blank(),
                plot.margin = unit(c(3,0,-3,3), "cm")
        )+theme(legend.title = element_text(size = 12),
                legend.text = element_text(size=12))+
        scale_x_continuous(expand = c(0,1))
middle

rt1=rt[c(4:(ncol(rt)-3))] 
max(rt1)
min(rt1)
pheatmap::pheatmap(t(rt1),scale = "row")
rt1=rt1 %>% scale() %>%as.data.frame()
rt1$id = 1:nrow(rt1)
rt1= tidyr::gather(rt1,variable,value,-id)
rt1$variable <- factor(rt1$variable ,levels =c("CD38","PYGM","ACSS3","CA14","ALDH1A2","ISYNA1","HAGHL","AK5","TK1","RRM2") )
p3 = ggplot(rt1, aes_string(x = 'id',y = 'variable',fill = 'value')) +
        geom_raster() + labs(fill="Expression")+
        theme(
                panel.grid = element_blank(),
                panel.background = element_blank(),
                axis.line = element_blank(),
                axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                axis.title = element_blank(),
                plot.background = element_blank() #the key to avoide legend overlap
        ) +
        scale_fill_gradient2(
                low = "blue",
                mid = "white",
                high = "red"
        ) +
        scale_x_continuous(expand = c(0,0))
p3


relative_heights=c(1,0.06,1)
library(aplot)
p1 %>% 
        insert_bottom(p2,height = relative_heights[1])%>% 
        insert_bottom(middle,height=relative_heights[2]) %>%
        insert_bottom(p3,height=relative_heights[3])

ggsave(file = paste0("08ROC/",outputfile,"_riskCurve.pdf") , width = 5,height = 5.5)
