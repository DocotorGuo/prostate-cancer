#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
rm(list = ls())  ## 一键清空~
options(stringsAsFactors = F)
expFile="01rawData/TCGA-PRAD_TPM.txt"                                                 
pdFile = "03NMF/TCGA-PRAD_NMF_Cluster.txt"

#read gene expression files and clean the data
testExpr <- read.table(expFile,sep="\t",header=T,row.names=1,check.names=F)
xx <- substr(colnames(testExpr),14,16)
table(xx)
TumorId <- colnames(testExpr)[substr(colnames(testExpr),14,15)=="01"]
testExpr<- testExpr[,TumorId]

# colnames(testExpr) <- substr(colnames(testExpr),1,12 )
pd<- read.table(pdFile,sep="\t",header=T,row.names = 1,check.names=F)

# Drug prediction R package of pRRophetic  
library(pRRophetic)
library(ggplot2)
library(ggpubr)
library(cowplot)
mycol <- c( '#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")

#Drug name
GCP.drug <- read.table("06immunotherapy/pRRophetic/drug_sig.txt") 
GCP.drug <- GCP.drug$V1

exprData <- as.matrix(testExpr[,rownames(pd)])
# Customize as many box colors as possible
jco <- c('#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")

### Drug sensitivity prediction、
GCPinfo <- GCP.IC50 <- GCP.expr <- cvOut <- predictedPtype <- predictedBoxdat <- list() # 初始化列表
plotp <- list()

for (drug in GCP.drug) {
  set.seed(4) #set the seeds so that the results can be repeated
  cat(drug," starts!\n")  
  # 预测IC50值，采用默认参数，详细可参考??pRRopheticPredict参数列表
  predictedPtype[[drug]] <- pRRopheticPredict(testMatrix = exprData ,
                                              drug = drug,
                                              tissueType = "allSolidTumors",
                                              batchCorrect = "eb",
                                              selection = 1, # 1Means if there are duplicate genes, take mean treatment
                                              dataset = "cgp2014")
  if(!all(names(predictedPtype[[drug]])==rownames(pd))) {stop("Name mismatched!\n")} # If the names do not match, an error is reported and exit
  
  predictedBoxdat[[drug]] <- data.frame("est.ic50"=predictedPtype[[drug]],
                                        "Cluster"=ifelse(pd$Cluster == "Cluster1","Cluster1","Cluster2"), 
                                        row.names = names(predictedPtype[[drug]])) 
  predictedBoxdat[[drug]]$Cluster <- factor(predictedBoxdat[[drug]]$Cluster,levels = c("Cluster1","Cluster2"),ordered = T) 
  
  # plot
  p <- ggplot(data = predictedBoxdat[[drug]], aes(x=Cluster, y=est.ic50))
  p <- p + geom_boxplot(aes(fill = Cluster),outlier.colour = NA,notch = T,size = 0.3)+
    geom_jitter(aes(fill = Cluster),shape = 21,size=2,width = 0.2)+
    geom_violin(aes(fill = Cluster),position = position_dodge(width = .75), 
                size = 0.3,alpha = 0.4,trim = T)+
    scale_fill_manual(values=mycol[c(3,2)]) + #自定义box的配色
    theme_classic()+
    stat_compare_means(comparisons = list(c("Cluster1","Cluster2")),method="wilcox.test", label = "p.format")+
    theme(legend.position="none") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),plot.title = element_text(size = 12, hjust = 0.5)) + 
    xlab("") + ylab("Estimated IC50") + 
    ggtitle(drug) # title
  
  plotp[[drug]] <- p 
  cat(drug," has been finished!\n") 
}


# Suitable for displaying two drugs
p1 <- plot_grid(plotp[[1]],plotp[[2]],plotp[[3]],plotp[[4]],labels = c("A","B","C","D"),nrow = 1) 
ggsave("06immunotherapy/pRRophetic/TOP4 boxplot of predicted IC50.pdf", width = 8, height = 3.6)

# Suitable for displaying a variety of drugs
p2 <- plot_grid(plotlist=plotp, ncol=5)
ggsave("06immunotherapy/pRRophetic/boxplot of predicted IC50_multiple.pdf", width = 10, height = 26)

# The differences between groups were examined
p <- vector()
for (drug in GCP.drug) {
  tmp <- wilcox.test(as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$Cluster %in% "Cluster1"),"est.ic50"]),
                     as.numeric(predictedBoxdat[[drug]][which(predictedBoxdat[[drug]]$Cluster %in% "Cluster2"),"est.ic50"]),alternative = "less")$p.value
  p <- append(p,tmp)
}
names(p) <- GCP.drug
print(p) #
#保存到文件
write.table(p,"06immunotherapy/pRRophetic/drug_sig_pvalue.txt", quote = F, sep = "\t")



#tide 通过TIDE比较不同亚型免疫应答情况
library(pheatmap)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
TIDE <- testExpr
# 这里为了得到比较好的结果，采用two direction median centered
TIDE <- sweep(TIDE,2, apply(TIDE,2,median,na.rm=T))
TIDE <- sweep(TIDE,1, apply(TIDE,1,median,na.rm=T))
write.table(TIDE,"06immunotherapy/TIDE_input.self_subtract",sep = "\t",row.names = T,col.names = NA,quote = F)



# See the TIDE tutorial in the folder to get the output file tide_output.csv  
TIDE.res <- read.csv("06immunotherapy/SubMap/TIDE_output.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)

pd$TIDE <- TIDE.res[rownames(pd),"Responder"]
print(table(pd$TIDE,pd$Cluster))

print(fisher.test(table(pd$TIDE,pd$Cluster))) 

# submap预测亚型的免疫治疗响应性
# 自定义函数用来产生submap需要的数据格式
generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct, 
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
  
  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

# 创建submap需要的数据格式 (SKCM)
skcm.immunotherapy.logNC <- read.table("06immunotherapy/SubMap/skcm.immunotherapy.47samples.log2CountsNorm.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F) #原文提供的log2转化的标准化count值
rownames(skcm.immunotherapy.logNC) <- toupper(rownames(skcm.immunotherapy.logNC)) # 基因大写，因为我使用的数据是把基因名都大写的
skcm.immunotherapy.info <- read.table("06immunotherapy/SubMap/skcm.immunotherapy.47sampleInfo.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

skcm.immunotherapy.info <- skcm.immunotherapy.info[order(skcm.immunotherapy.info$label),]
skcm.immunotherapy.info$rank <- rep(c(1,2,3,4),times=as.character(table(skcm.immunotherapy.info$label))) #1: CTLA4_noR 2: CTLA4_R 3:PD1_noR 4:PD1_R

tmp <- testExpr
# colnames(tmp) <- substr(colnames(tmp),1,12 )
GENELIST <- intersect(rownames(tmp),rownames(skcm.immunotherapy.logNC)) # 取交集后的基因列表

sam_info <- skcm.immunotherapy.info
in_gct <- skcm.immunotherapy.logNC[GENELIST,rownames(skcm.immunotherapy.info)]

# The file name that generates the output data
gct_file <- "skcm.immunotherapy.for.SubMap.gct"
cls_file <- "skcm.immunotherapy.for.SubMap.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

# Propose subtypes of samples and arrange them sequentially
samples.C1 <- rownames(pd)[which(pd$Cluster == "Cluster1")]
samples.C2 <- rownames(pd)[which(pd$Cluster == "Cluster2")]

sam_info <- data.frame("Type"=c(samples.C1,samples.C2),row.names = c(samples.C1,samples.C2))
sam_info$rank <- rep(c(1,2),times=c(length(samples.C1),length(samples.C2))) #1: C1,即HPV16-IMM 2: C2,即HPV16-KRT

# The file name that generates the output data
gct_file <- "TCGA-PRAD.Immune2.for.SubMap.gct"
cls_file <- "TCGA-PRAD.Immune2.for.SubMap.cls"

in_gct <- log2(tmp[GENELIST,rownames(sam_info)] + 1) # Produces a form similar to the sample data, the normalized count value converted by log2  
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")


library(GenePattern)
gctfileA <- '06immunotherapy/SubMap/TCGA-PRAD.Immune2.for.SubMap.gct'
# a <- read.gct(gctfileA)
clsfileA <- '06immunotherapy/SubMap/TCGA-PRAD.Immune2.for.SubMap.cls'
# b <- read.cls(clsfileA)

gctfileB <- '06immunotherapy/SubMap/skcm.immunotherapy.for.SubMap.gct'
# a <- read.gct(gctfileB)
clsfileB <- '06immunotherapy/SubMap/skcm.immunotherapy.for.SubMap.cls'
# b <- read.cls(clsfileB)

source('D:/GEO/script/submap.R')
submap.main(  gctfileA,
              gctfileB,
              clsfileA,
              clsfileB,
              output.filename="06immunotherapy/SubMap/",
              ntag=100,
              nperm=50,
              nperm.fisher=1000,
              weighted.score.type=1,
              null.dist="pool",
              p.corr="Bonferroni",
              clust.row=1,
              clust.col=1,
              nom.p.mat="T",
              create.legend="T",
              # rnd.seed=47365321
              rnd.seed=5321) 

library(pheatmap)
heatmap.YlGnPe <- c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F")
cherry    <- "#700353"
lightgrey <- "#dcddde"

# Put the values in the submap result submap_submapresult. TXT file in the appropriate places  
# Input the nominal p values in the file and the corrected P values in the heat map  
tmp <- matrix(c(0.987,0.412,0.997,0.001,0.040,0.388,0.086,0.976, 
                1,1,1,0.008,0.320,1,0.691,1), # Bonferroni校正p值
              nrow = 4,byrow = T,dimnames = list(c("Cluster1_p","Cluster2_p","Cluster1_b","Cluster2_b"),c("CTAL4-noR","CTLA4-R","PD1-noR","PD1-R")))

library(pheatmap)
submap.data.A <- matrix(c( 0.799200799, 0.93206793, 0.2537463, 0.7962038,
                           0.004995005, 0.04095904, 0.8841159, 0.2547453),nrow = 2,byrow = T)
submap.data.B <- matrix(c(1.00000000, 1.0000000,  1,  1,
                          0.03996004, 0.3276723,  1,  1),nrow = 2,byrow = T)
submap.data <- rbind(submap.data.A,submap.data.B)
row.names(submap.data) <- c('Cluster1_p','Cluster2_p','Cluster1_b','Cluster2_b')
colnames(submap.data) <- c("CTAL4-noR","CTLA4-R","PD1-noR","PD1-R")
annotation_row <- data.frame(pvalue=rep(c('Nominal p value','Bonferroni corrected'),c(2,2)))
rownames(annotation_row) <- row.names(submap.data)
pdf(file = '06immunotherapy/SubMap/subclass mapping.pdf',width= 12,height= 8)
pheatmap(submap.data,
         show_colnames = T,
         color = heatmap.YlGnPe[5:1],
         display_numbers = matrix(ifelse(submap.data < 0.05,'p < .05',''),nrow(submap.data)),number_format = "%.3f",
         cluster_rows = F,cluster_cols = F,
         annotation_row=annotation_row,
         annotation_colors = list(pvalue=c("Nominal p value"=lightgrey,"Bonferroni corrected"=cherry)),
         cellwidth=30,cellheight=30,main = "",fontsize_number = 9,
         #fontsize_row = 20,fontsize_col = 25,
         gaps_row=2,
         fontsize=12)
dev.off()


# Draw a boxplot of the TIde response
ann <-read.table("03NMF/TCGA-PRAD_NMF_Cluster.txt",header = T,sep = "\t",check.names = F)
head(ann)
print(table(ann$Cluster))

TIDE.res <- read.csv("06immunotherapy/SubMap/TIDE_output.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)

TIDE.res$id <- rownames(TIDE.res)


# To test whether immunotherapy response was correlated with subtype, P <0.05 indicated correlation  
ann$TIDE <- TIDE.res[ann$id,"Responder"]
print(table(ann$TIDE,ann$Cluster))
print(fisher.test(table(ann$TIDE,ann$Cluster))) 



boxdata <- merge(TIDE.res,ann,by="id")
boxdata$Cluster <- factor(boxdata$Cluster,levels=c("Cluster1","Cluster2"))


library(ggpubr)
ggplot(boxdata,aes(Cluster,TIDE.x,fill=Cluster))+
  geom_boxplot(outlier.colour = NA,notch = T,size = 0.3)+
  geom_jitter(shape = 21,size=2,width = 0.2)+
  geom_violin(position = position_dodge(width = .75), 
              size = 0.3,alpha = 0.4,trim = T)+
  theme_classic()+ ylab('TIDE')+ggtitle('TIDE prediction')+
  theme(legend.position = 'none',
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title.x = element_blank())+
  scale_fill_manual(values=mycol[c(3,2)])+
  stat_compare_means(comparisons = list(c("Cluster1","Cluster2")),
                     label = 'p.signif')+
  stat_compare_means(label.y = max(boxdata$TIDE.x)+0)
ggsave("06immunotherapy/SubMap/TIDE_boxplot.pdf",width = 3.5,height = 3.5)



#Wilcoxon rank-sum test was drawn to compare the expression of immunoassay sites of different subtypes  （PDCD1, CD274, PDCD1LG2, TIGIT, HAVCR2, ID
expr <-  read.table("01rawData/TCGA-PRAD_TPM.txt",header = T,sep = "\t",check.names = F,row.names = 1)
xx <- substr(colnames(expr),14,16)
table(xx)
TumorId <- colnames(expr)[substr(colnames(expr),14,15)=="01"]
expr<- expr[,TumorId]

genes <-  read.table("06immunotherapy/HLA_genes.txt",header = T,sep = "\t",check.names = F)

Cluster_res <- read.table("03NMF/TCGA-PRAD_NMF_Cluster.txt",header = T,sep = "\t",check.names = F)
Cluster1 <- Cluster_res$id[Cluster_res$Cluster =="Cluster1"]
Cluster2 <- Cluster_res$id[Cluster_res$Cluster =="Cluster2"]

scores <- expr[rownames(expr)%in%genes$gene, ]
scores <-as.matrix(scores[,c(Cluster1,Cluster2)])


boxdata <- as.data.frame(t( (scores) )) %>% rownames_to_column("id")
boxdata <- as.data.frame(t( log2(scores+1) )) %>% rownames_to_column("id")
Cluster <- read.table("03NMF/TCGA-PRAD_NMF_Cluster.txt",header = T,sep = "\t",check.names = F)


boxdata <- merge(boxdata,Cluster[,c("id","Cluster")],by="id")

write.table(boxdata, "06immunotherapy/HLA_genes_type.txt", quote = F, row.names = F,sep = "\t")


boxdata1 <- boxdata[,-1]
library(tidyr)
library(ggpubr)
mycol <- c('#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")

boxdata1 <- gather(boxdata1,cellType, value,colnames(boxdata1)[-ncol(boxdata1)])
boxdata1$Cluster <- factor(boxdata1$Cluster,levels = c("Cluster1","Cluster2"))

p2 <- ggboxplot(boxdata1, x = "cellType", y = "value",
                color = "Cluster", palette = mycol,
                add = "jitter",
                add.params = list(alpha=0.6),)+ rotate_x_text(45) +xlab("")+ylab("Gene expression log2(TPM+1)")
p2 + stat_compare_means(aes(group = Cluster),method="wilcox.test", label = "p.signif")
# p2 + stat_compare_means(comparisons = split(t(combn(levels(boxdata1$Cluster),2)),1:nrow(t(combn(levels(boxdata1$Cluster),2)))))
ggsave("06immunotherapy/HLA_genes_type_boxplot.pdf",height=4.5,width=7.5)


genes <-  read.table("06immunotherapy/immune-checkpoint_genes.txt",header = T,sep = "\t",check.names = F)

Cluster_res <- read.table("03NMF/TCGA-PRAD_NMF_Cluster.txt",header = T,sep = "\t",check.names = F)
Cluster1 <- Cluster_res$id[Cluster_res$Cluster =="Cluster1"]
Cluster2 <- Cluster_res$id[Cluster_res$Cluster =="Cluster2"]

scores <- expr[rownames(expr)%in%genes$gene, ]
scores <-as.matrix(scores[,c(Cluster1,Cluster2)])

# scores =read.table("GSVA/GSVA_geneSet_result.txt",header=T,sep="\t",check.names=F,row.names = 1)
boxdata <- as.data.frame(t( (scores) )) %>% rownames_to_column("id")
boxdata <- as.data.frame(t( log2(scores+1) )) %>% rownames_to_column("id")
Cluster <- read.table("03NMF/TCGA-PRAD_NMF_Cluster.txt",header = T,sep = "\t",check.names = F)


# boxdata$id <- substr(rownames(boxdata) ,1,12)
# Cluster$Cluster <- ifelse(Cluster$Cluster =="Cluster2","Cluster2","Cluster1 + Cluster3")
boxdata <- merge(boxdata,Cluster[,c("id","Cluster")],by="id")

write.table(boxdata, "06immunotherapy/immune-checkpoint_genes_type.txt", quote = F, row.names = F,sep = "\t")


boxdata1 <- boxdata[,-1]
library(tidyr)
library(ggpubr)
mycol <- c('#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")

boxdata1 <- gather(boxdata1,cellType, value,colnames(boxdata1)[-ncol(boxdata1)])
boxdata1$Cluster <- factor(boxdata1$Cluster,levels = c("Cluster1","Cluster2"))

p2 <- ggboxplot(boxdata1, x = "cellType", y = "value",
                color = "Cluster", palette = mycol,
                add = "jitter",
                add.params = list(alpha=0.6),)+ rotate_x_text(45) +xlab("")+ylab("Gene expression log2(TPM+1)")
p2 + stat_compare_means(aes(group = Cluster),method="wilcox.test", label = "p.signif")
# p2 + stat_compare_means(comparisons = split(t(combn(levels(boxdata1$Cluster),2)),1:nrow(t(combn(levels(boxdata1$Cluster),2)))))
ggsave("06immunotherapy/immunecheck_genes_type_boxplot.pdf",height=4.5,width=7)
