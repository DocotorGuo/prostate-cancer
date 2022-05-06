#install.packages("survivalROC")

library(survivalROC)

rt=read.table("09Independent/trianSet_riskScore.txt",header=T,sep="\t",check.names=F,row.names=1)    #??ȡcox?ع??????ļ?
# rt$futime=rt$futime/365
rt <- rt[-c(8:17,ncol(rt))]
colnames(rt)[6:7] <- c("PSA","gleason")
rt <- na.omit(rt)
rt$pathologic_T <- as.numeric(as.factor(rt$pathologic_T))
rt$pathologic_N <- as.numeric(as.factor(rt$pathologic_N))

rocCol <- c('#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#8DD3C7',"#999999")

aucText=c()

#risk score ROC
pdf(file="09Independent/multiROC.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$Disease_Free_Survival, status=rt$Disease_Free_Status, marker = rt$riskScore, predict.time =12, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
  xlab="1-Specificity (False Positive)", ylab="Sensitivity (True Positive)",
  lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("risk score"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)


j=1
for(i in colnames(rt[,3:(ncol(rt)-1)])){
	roc=survivalROC(Stime=rt$Disease_Free_Survival, status=rt$Disease_Free_Status, marker = rt[,i], predict.time =12, method="KM")
	j=j+1
	aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
	lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 2)
}
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
dev.off()

