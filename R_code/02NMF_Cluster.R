rm(list = ls())
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(survival)
library(survminer)
library(Rtsne)
library(ggplot2)
library(ggsci)
library(pheatmap)
mycol <- c('#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")
# Differential metabolic genes were read for NFM cluster analysis
feExpr <- read.table("03Cox/TCGA_PRAD_diff_MGs_exprSet.txt",sep = "\t",header = T,row.names = 1,check.names = F)


feExpr_matrix <- as.matrix(log2(feExpr+1) )
max(feExpr_matrix);min(feExpr_matrix)

library(NMF)
library(doMPI) 

ranks <- 2:6
nmf.res <- nmf(feExpr_matrix,ranks,nrun=50,seed = 4)
pdf(file="03NMF/NMF.pdf", width = 7.5,height = 7)
plot(nmf.res)
dev.off()
#Run the NMF again with a better rank
nmf.rank3 <- nmf(feExpr_matrix,
                 rank = 2,
                 method = "brunet",
                 nrun=50,
                 seed = 4)

index <- extractFeatures(nmf.rank3,"max")
sig.order <- unlist(index)
NMF.Exp.rank3 <- feExpr_matrix[sig.order,]
NMF.Exp.rank3 <- na.omit(NMF.Exp.rank3) 
Cluster <- predict(nmf.rank3) 
table(Cluster)

# Save the clustering results
Cluster_res <- Cluster %>%
  as.data.frame() %>%
  rownames_to_column("id")# %>% arrange(Cluster)
colnames(Cluster_res) <-c("id","Cluster")
Cluster_res$Cluster <- ifelse(Cluster_res$Cluster==1,"Cluster1","Cluster2")
table(Cluster_res$Cluster)
write.table(Cluster_res,file = "03NMF/TCGA-PRAD_NMF_Cluster.txt",row.names = F,quote = F,sep = "\t")
pdf(file="03NMF/NMF_heatmap.pdf", width = 5,height = 5)
consensusmap(nmf.rank3,
             labRow = NA,
             labCol = NA,
             annCol = data.frame("cluster"=Cluster[colnames(NMF.Exp.rank3)]),
             annColors = list(cluster=c("1"=mycol[1],"2"=mycol[2])))
dev.off()


# Load Rtsne package to validate the performances of NMF
library(Rtsne)
set.seed(4)
# TSNE dimension reduction analysis was performed using Rtsne function  
feExpr1 <- as.matrix(t(feExpr_matrix))
tsne_out <- Rtsne(feExpr1,pca=FALSE,dims=2,
                  perplexity=15,theta=0.0) # Run TSNE
plot(tsne_out$Y,col=factor( Cluster_res$Cluster),asp=1)
tsne_plot  = data.frame(tSNE1  =tsne_out$Y[,1], tSNE2  = tsne_out$Y[,2],Cluster = Cluster_res$Cluster)
tsne_plot$Cluster <- factor(tsne_plot$Cluster ,labels = c("Cluster1","Cluster2"))
head(tsne_plot)
library(paletteer)
ggplot(tsne_plot,aes(x=tSNE1,y=tSNE2))+
    geom_point(aes(fill = Cluster),shape = 21,color = "black")+
    scale_fill_manual(values=mycol[1:2]) +
    # stat_ellipse(aes(color = Cluster,fill = Cluster),
    #              geom = "polygon",
    #              alpha = 0.3,
    #              linetype = 2)+
    # scale_color_paletteer_d("RColorBrewer::Set3")+
    # scale_fill_paletteer_d("RColorBrewer::Set3")+
    theme_classic()+
    theme(legend.position = "top")

ggsave("03NMF/TCGA-PRAD_Cluster_tSNE.pdf",width = 4.2,height = 4.5)

Cluster_res$id <- substr(Cluster_res$id,1,12)
Cluster_res <- Cluster_res %>%
    inner_join(clincal, "id")
head(Cluster_res)
write.table(Cluster_res,file = "03NMF/TCGA-PRAD_NMF_Cluster_survival.txt",row.names = F,quote = F,sep = "\t")

# The clustering results were analyzed for survival
library(survminer)
library(survival)
fit <- survfit(Surv(Disease_Free_Survival , Disease_Free_Status ) ~ Cluster, 
               data = Cluster_res)
ggsurvplot(fit,
           pval = TRUE,
           linetype = "solid",  
           palette = mycol[1:2],
           #surv.median.line = "hv", 
           title = "TCGA-PRAD",
           ylab = "Disease Free Survival",
           xlab = " Time (Months)",
           # legend.title = "Survival Plot",
           legend = c(0.8,0.30),
           legend.labs = c("Cluster1","Cluster2"),
           legend.title="",
           risk.table = T,
           risk.table.title="",
           tables.height = 0.2,
           ggtheme = theme_bw(),
           tables.theme = theme_cleantable()
)
dev.copy2pdf(file = "03NMF/TCGA-PRAD_Cluster_DFS_survival.pdf", width = 5,height = 5.5)
dev.off()
