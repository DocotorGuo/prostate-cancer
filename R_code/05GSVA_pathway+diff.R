# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
# BiocManager::install("ComplexHeatmap")
 
# BiocManager::install("GSVA")
mycol <- c('#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")

library(GSVA)
library(GSEABase)
library(limma)
library(tidyverse)
library(clusterProfiler)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

exprFile <- "01rawData/TCGA-PRAD_TPM.txt" 
gmtFile <- "05GSVA/c2.cp.kegg.v7.4.symbols.gmt"
expr <- read.table(exprFile,sep="\t",header=T,check.names=F,row.names = 1)
xx <- substr(colnames(expr),14,16)
table(xx)
TumorId <- colnames(expr)[substr(colnames(expr),14,15)=="01"]
expr<- expr[,TumorId]

geneSet=getGmt(gmtFile, 
                collectionType=BroadCollection(category="c2"),
               geneIdType=SymbolIdentifier())

gsva_result <- gsva(expr=as.matrix(expr), gset.idx.list=geneSet, 
               min.sz=5, 
               max.sz=500, 
               kcdf="Gaussian", 
               method='ssgsea',abs.ranking=TRUE,parallel.sz=4)
gsvaOut <- gsva_result %>% as.data.frame() %>% rownames_to_column("id")
gsvaOut[1:3,1:5]
write.table(gsvaOut,file="05GSVA/GSVA_KEGG_result.txt",row.names = FALSE, sep = "\t", quote = FALSE)


# wilcox.test pathway difference
samples_input <- "03NMF/TCGA-PRAD_NMF_Cluster.txt"
samples_file <- read.table(samples_input,sep="\t",header=T,check.names=F)
gsva_result[1:3,1:4]
# colnames(gsva_result) <- substr(colnames(gsva_result),1,12)
gsva_result_df <- t(gsva_result ) %>% as.data.frame() %>% rownames_to_column("id")%>% inner_join(samples_file[,c("id","Cluster")], "id")

outTab <- data.frame()
for (i in 2:(ncol(gsva_result_df)-1)){
  rt <- data.frame(expression=gsva_result_df[,i],Cluster=gsva_result_df[,"Cluster"])
  # kruskalTest <- kruskal.test(get(pathway)~Cluster,gsva_result_df)
  wilcoxTest <- wilcox.test(expression~Cluster,rt)
  pvalue <- wilcoxTest$p.value
  Cluster1mean=mean(rt[rt$Cluster=="Cluster1","expression"])
  Cluster2mean=mean(rt[rt$Cluster=="Cluster2","expression"])
  diffMed=log2(Cluster1mean+1) -log2(Cluster2mean+1) 
  
  outTab <- rbind(outTab,
                  cbind(id=colnames(gsva_result_df)[i],
                        Cluster1Mean=Cluster1mean,
                        Cluster2Mean=Cluster2mean,
                        log2FC=diffMed,
                        p.Value=pvalue))
  }
outTab$p.adjust =p.adjust(outTab$p.Value, method = "BH")
write.table(outTab, file = "05GSVA/GSVA_KEGG_wilcoxTest.txt", row.names = FALSE, sep = "\t", quote = FALSE)



gmtFile <- "05GSVA/h.all.v7.4.symbols.gmt"
geneSet=getGmt(gmtFile, 
               collectionType=BroadCollection(category="h"),
               geneIdType=SymbolIdentifier())

gsva_result <- gsva(expr=as.matrix(expr), gset.idx.list=geneSet, 
                    min.sz=5, 
                    max.sz=500, 
                    kcdf="Gaussian", 
                    method='ssgsea',abs.ranking=TRUE,parallel.sz=4)
gsvaOut <- gsva_result %>% as.data.frame() %>% rownames_to_column("id")
gsvaOut[1:3,1:5]
write.table(gsvaOut,file="05GSVA/GSVA_HALLMARK_result.txt",row.names = FALSE, sep = "\t", quote = FALSE)



# wilcox.test pathway difference
samples_input <- "03NMF/TCGA-PRAD_NMF_Cluster.txt"
samples_file <- read.table(samples_input,sep="\t",header=T,check.names=F)
gsva_result[1:3,1:4]
# colnames(gsva_result) <- substr(colnames(gsva_result),1,12)
gsva_result_df <- t(gsva_result ) %>% as.data.frame() %>% rownames_to_column("id")%>% inner_join(samples_file[,c("id","Cluster")], "id")

outTab <- data.frame()
for (i in 2:(ncol(gsva_result_df)-1)){
  rt <- data.frame(expression=gsva_result_df[,i],Cluster=gsva_result_df[,"Cluster"])
  # kruskalTest <- kruskal.test(get(pathway)~Cluster,gsva_result_df)
  wilcoxTest <- wilcox.test(expression~Cluster,rt)
  pvalue <- wilcoxTest$p.value
  Cluster1mean=mean(rt[rt$Cluster=="Cluster1","expression"])
  Cluster2mean=mean(rt[rt$Cluster=="Cluster2","expression"])
  diffMed=log2(Cluster1mean+1) -log2(Cluster2mean+1) 
  
  outTab <- rbind(outTab,
                  cbind(id=colnames(gsva_result_df)[i],
                        Cluster1Mean=Cluster1mean,
                        Cluster2Mean=Cluster2mean,
                        log2FC=diffMed,
                        p.Value=pvalue))
}
outTab$p.adjust =p.adjust(outTab$p.Value, method = "BH")
write.table(outTab, file = "05GSVA/GSVA_HALLMARK_wilcoxTest.txt", row.names = FALSE, sep = "\t", quote = FALSE)


#Customized gene sets were analyzed for GSVA enrichment
gmtFile <- "05GSVA/GSVA_geneSet.txt"

gmt <- read.table(gmtFile,sep="\t",check.names=F,row.names = 1)
geneSetlist <- list()
geneSet <- as.data.frame(t(gmt))
for (i in colnames(geneSet)) {
  geneSetlist[[i]] <-  geneSet[,i][which(geneSet[,i]!="")]

}
geneSet <- geneSetlist

gsva_result <- gsva(expr=as.matrix(expr), gset.idx.list=geneSet, 
                    min.sz=1, 
                    max.sz=500, 
                    kcdf="Gaussian", 
                    method='ssgsea',abs.ranking=TRUE,parallel.sz=4)
gsvaOut <- gsva_result %>% as.data.frame() %>% rownames_to_column("id")
gsvaOut[1:3,1:5]
write.table(gsvaOut,file="05GSVA/GSVA_geneSet_result.txt",row.names = FALSE, sep = "\t", quote = FALSE)


#  wilcox.test pathway difference
samples_input <- "03NMF/TCGA-PRAD_NMF_Cluster.txt"
samples_file <- read.table(samples_input,sep="\t",header=T,check.names=F)
gsva_result[1:3,1:4]
# colnames(gsva_result) <- substr(colnames(gsva_result),1,12)
gsva_result_df <- t(gsva_result ) %>% as.data.frame() %>% rownames_to_column("id")%>% inner_join(samples_file[,c("id","Cluster")], "id")

outTab <- data.frame()
for (i in 2:(ncol(gsva_result_df)-1)){
  rt <- data.frame(expression=gsva_result_df[,i],Cluster=gsva_result_df[,"Cluster"])
  # kruskalTest <- kruskal.test(get(pathway)~Cluster,gsva_result_df)
  wilcoxTest <- wilcox.test(expression~Cluster,rt)
  pvalue <- wilcoxTest$p.value
  Cluster1mean=mean(rt[rt$Cluster=="Cluster1","expression"])
  Cluster2mean=mean(rt[rt$Cluster=="Cluster2","expression"])
  diffMed=log2(Cluster1mean+1) -log2(Cluster2mean+1) 
  
  outTab <- rbind(outTab,
                  cbind(id=colnames(gsva_result_df)[i],
                        Cluster1Mean=Cluster1mean,
                        Cluster2Mean=Cluster2mean,
                        logFC=diffMed,
                        p.Value=pvalue))
}
outTab$p.adjust =p.adjust(outTab$p.Value, method = "BH")
write.table(outTab, file = "05GSVA/GSVA_geneSet_wilcoxTest.txt", row.names = FALSE, sep = "\t", quote = FALSE)

#Heatmaps of GSVA score differences of the custom gene set pathway were compared
exprSet_m6A <-  read.table("05GSVA/GSVA_geneSet_result.txt",header = T,sep = "\t",check.names = F,row.names = 1)

outDiff.mRNA <-  read.table("05GSVA/GSVA_geneSet_wilcoxTest.txt",header = T,sep = "\t",check.names = F,row.names = 1)

geneNum=30
# outDiff.mRNA=outDiff.mRNA[outDiff.mRNA$p.adjust < 0.05,]
outDiff.mRNA=outDiff.mRNA[order(as.numeric(as.vector(outDiff.mRNA$logFC))),]
diffGeneName=rownames( outDiff.mRNA )
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(geneNum*2)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp= exprSet_m6A[hmGene,]

max(hmExp)
min(hmExp)
library(pheatmap)
Cluster_res <- read.table("03NMF/TCGA-PRAD_NMF_Cluster.txt",header = T,sep = "\t",check.names = F)
Cluster1 <- Cluster_res$id[Cluster_res$Cluster =="Cluster1"]
Cluster2 <- Cluster_res$id[Cluster_res$Cluster =="Cluster2"]

hmExp <-(hmExp[,c(Cluster1,Cluster2)])


annotation_col = data.frame(
  Type = factor(rep(c("Cluster1","Cluster2"), c(length(Cluster1) ,length(Cluster2))),levels = c("Cluster1","Cluster2"))
)
rownames(annotation_col) = colnames(hmExp)
head(annotation_col)
# 自定注释信息的颜色列表 "#E69F00","#56B4E9"
ann_colors = list(
  Type = c(Cluster1 = mycol[1], Cluster2 = mycol[2] )
  )
head(ann_colors)
pdf(file="05GSVA/GSVA_geneSet_heatmap.pdf",height=5,width=9)
pheatmap(hmExp,
  scale = "row",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_colnames = F, 
  cluster_cols =F,
  border=FALSE,
  # cellwidth = 4,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()


#Heatmaps of GSVA score differences of the HALLMARK pathway were compared
exprSet_m6A <-  read.table("05GSVA/GSVA_HALLMARK_result.txt",header = T,sep = "\t",check.names = F,row.names = 1)
outDiff.mRNA <-  read.table("05GSVA/GSVA_HALLMARK_wilcoxTest.txt",header = T,sep = "\t",check.names = F,row.names = 1)
outDiff.mRNA=outDiff.mRNA[outDiff.mRNA$p.adjust < 0.05,]
geneNum=10
# outDiff.mRNA=outDiff.mRNA[outDiff.mRNA$p.adjust < 0.05,]
outDiff.mRNA=outDiff.mRNA[order(as.numeric(as.vector(outDiff.mRNA$log2FC))),]
diffGeneName=rownames( outDiff.mRNA )
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(geneNum*2)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp= exprSet_m6A[hmGene,]


library(pheatmap)
Cluster_res <- read.table("03NMF/TCGA-PRAD_NMF_Cluster.txt",header = T,sep = "\t",check.names = F)
Cluster1 <- Cluster_res$id[Cluster_res$Cluster =="Cluster1"]
Cluster2 <- Cluster_res$id[Cluster_res$Cluster =="Cluster2"]

hmExp <-(hmExp[,c(Cluster1,Cluster2)])

annotation_col = data.frame(
  Type = factor(rep(c("Cluster1","Cluster2"), c(length(Cluster1) ,length(Cluster2))),levels = c("Cluster1","Cluster2"))
)
rownames(annotation_col) = colnames(hmExp)
head(annotation_col)

ann_colors = list(
  Type = c(Cluster1 = mycol[1], Cluster2 = mycol[2] ) 
)
head(ann_colors)

pdf(file="05GSVA/GSVA_HALLMARK_heatmap.pdf",height=5,width=9)
pheatmap(hmExp,
         scale = "row",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = F, 
         cluster_cols =F,
         border=FALSE,
         # cellwidth = 4,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()


scores =read.table("05GSVA/GSVA_geneSet_result.txt",header=T,sep="\t",check.names=F,row.names = 1)
boxdata <- as.data.frame(t(scores)) %>% rownames_to_column("id")
# boxdata$id <- substr(boxdata$id ,1,12)
risk <- read.table("03NMF/TCGA-PRAD_NMF_Cluster.txt",header = T,sep = "\t",check.names = F)

boxdata <- merge(boxdata,risk[,c("id","Cluster")],by="id")

write.table(boxdata, "05GSVA/GSVA_geneSet_result_risk.txt", quote = F, row.names = F,sep = "\t")


boxdata1 <- boxdata[,-1]
library(tidyr)
library(ggpubr)
mycol <- c('#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")

boxdata1 <- gather(boxdata1,cellType, value,colnames(boxdata1)[-ncol(boxdata1)])
boxdata1$Cluster <- factor(boxdata1$Cluster,levels = c("Cluster1","Cluster2"))

p2 <- ggboxplot(boxdata1, x = "cellType", y = "value",
                color = "Cluster", palette = mycol[c(1:2)],
                add = "jitter",
                add.params = list(alpha=0.6),)+ rotate_x_text(45) +xlab("")+ylab("ssGSEA score")
p2 + stat_compare_means(aes(group = Cluster),method="wilcox.test", label = "p.signif")
# p2 + stat_compare_means(comparisons = split(t(combn(levels(boxdata1$Cluster),2)),1:nrow(t(combn(levels(boxdata1$Cluster),2)))))
ggsave("05GSVA/GSVA_geneSet_ssGSEA_score_boxplot.pdf",height=5.5,width=10)
