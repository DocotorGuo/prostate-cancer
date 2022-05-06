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


feExpr <- read.table("03Cox/TCGA_PRAD_diff_MGs_exprSet.txt",sep = "\t",header = T,row.names = 1,check.names = F)
MGExpr <- read.table("01rawData/GSE70768_symbol_exprSet.txt",sep = "\t",header = T,row.names = 1,check.names = F)
MGExpr <- MGExpr[  rownames(MGExpr) %in% rownames(feExpr),]

feExpr_matrix <- as.matrix((MGExpr) )

max(feExpr_matrix);min(feExpr_matrix)


library(NMF)
library(doMPI)  

ranks <- 2:6
nmf.res <- nmf(feExpr_matrix,ranks,nrun=50,seed = 4)
pdf(file="03NMF/GSE70768_NMF.pdf", width = 7.5,height = 7)
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

# Save the NFM clustering results
Cluster_res <- Cluster %>%
  as.data.frame() %>%
  rownames_to_column("id")# %>% arrange(Cluster)
colnames(Cluster_res) <-c("id","Cluster")
Cluster_res$Cluster <- ifelse(Cluster_res$Cluster==1,"Cluster1","Cluster2")
table(Cluster_res$Cluster)
write.table(Cluster_res,file = "03NMF/GSE70768_NMF_Cluster.txt",row.names = F,quote = F,sep = "\t")
pdf(file="03NMF/GSE70768_NMF_heatmap.pdf", width = 5,height = 5)
consensusmap(nmf.rank3,
             labRow = NA,
             labCol = NA,
             annCol = data.frame("cluster"=Cluster[colnames(NMF.Exp.rank3)]),
             annColors = list(cluster=c("1"=mycol[1],"2"=mycol[2])))
dev.off()



# Load Rtsne package  
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

ggsave("03NMF/GSE70768_Cluster_tSNE.pdf",width = 4.2,height = 4.5)




feExpr <- read.table("03Cox/TCGA_PRAD_diff_MGs_exprSet.txt",sep = "\t",header = T,row.names = 1,check.names = F)
MGExpr <- read.table("01rawData/GSE70769_symbol_exprSet.txt",sep = "\t",header = T,row.names = 1,check.names = F)
MGExpr <- MGExpr[  rownames(MGExpr) %in% rownames(feExpr),]

feExpr_matrix <- as.matrix((MGExpr) )
# feExpr_matrix <- as.matrix(log2(feExpr+1) )
max(feExpr_matrix);min(feExpr_matrix)


library(NMF)
library(doMPI)  

ranks <- 2:6
nmf.res <- nmf(feExpr_matrix,ranks,nrun=50,seed = 4)
pdf(file="03NMF/GSE70769_NMF.pdf", width = 7.5,height = 7)
plot(nmf.res)
dev.off()
#选用较优rank再运行一次NMF
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

# Save the NFM clustering results
Cluster_res <- Cluster %>%
  as.data.frame() %>%
  rownames_to_column("id")# %>% arrange(Cluster)
colnames(Cluster_res) <-c("id","Cluster")
Cluster_res$Cluster <- ifelse(Cluster_res$Cluster==1,"Cluster1","Cluster2")
table(Cluster_res$Cluster)
write.table(Cluster_res,file = "03NMF/GSE70769_NMF_Cluster.txt",row.names = F,quote = F,sep = "\t")
pdf(file="03NMF/GSE70769_NMF_heatmap.pdf", width = 5,height = 5)
consensusmap(nmf.rank3,
             labRow = NA,
             labCol = NA,
             annCol = data.frame("cluster"=Cluster[colnames(NMF.Exp.rank3)]),
             annColors = list(cluster=c("1"=mycol[1],"2"=mycol[2])))
dev.off()



# Load Rtsne package  
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

ggsave("03NMF/GSE70769_Cluster_tSNE.pdf",width = 4.2,height = 4.5)
