
#biocLite("GSVA")
#biocLite("limma")
#biocLite("GSEABase")

library(GSVA)
library(limma)
library(clusterProfiler)
expFile="01rawData/TCGA-PRAD_TPM.txt"                                                 
pdFile = "03NMF/TCGA-PRAD_NMF_Cluster.txt"

#The prostate cancer expression files were read and the data were collated
testExpr <- read.table(expFile,sep="\t",header=T,row.names=1,check.names=F)
xx <- substr(colnames(testExpr),14,16)
table(xx)
TumorId <- colnames(testExpr)[substr(colnames(testExpr),14,15)=="01"]
testExpr<- testExpr[,TumorId]

# colnames(testExpr) <- substr(colnames(testExpr),1,12 )
pd<- read.table(pdFile,sep="\t",header=T,check.names=F)
pd <- pd[order(pd$Cluster,decreasing = F),]
table(pd$Cluster)

# Propose subtypes of samples and arrange them sequentially
lowRSName <- pd$id[which(pd$Cluster == "Cluster1")]
highRSName <-  pd$id[which(pd$Cluster == "Cluster2")]

group <- c(rep("Cluster1", 264), rep("Cluster2", 225))
group <- paste(group, collapse = " ")
group <- c(paste(c(489, 2, 1), collapse = " "), "# Cluster1 Cluster2", group)

write.table(file = "07GSEA/Cluster_group.cls", group, col.names = F, row.names = F, quote = F)


data <- testExpr[,c(lowRSName,highRSName)]
GSEA_df <- cbind(Name = rownames(data) ,DESCRIPTION = "na",data)
# GSEA_df <- cbind(gene_id = rownames(data) ,data)
write.table(GSEA_df,"07GSEA/Cluster_Exprset.txt",row.names = F,sep = "\t",quote = F)
