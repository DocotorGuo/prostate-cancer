rm(list=ls())
## Load the R package
# BiocManager::install("paletteer")
source("D:/GEO/script/Filter_TCGA_Replicate_Samples.R")
library(limma)
library(tidyverse)

#Extracting clinical information
clinical <- read.table("01rawData/prad_tcga_clinical_data.tsv",check.names = F,sep="\t",quote = "",header = T)
clinical$id <- clinical$`Patient ID`
# Disease Free (Months)
clinical$Disease_Free_Survival <- clinical$`Disease Free (Months)`
clinical$Disease_Free_Status <- ifelse(clinical$`Disease Free Status` =="0:DiseaseFree",0,1)

clinical$age <- clinical$`Diagnosis Age`
clinical$pathologic_T <-  gsub("[ABCabc]","", clinical$`American Joint Committee on Cancer Tumor Stage Code`)
linical$pathologic_N <-  gsub("[ABCabc]","", clinical$`Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code`)
clinical$psa_value <- clinical$`Psa most recent results`
clinical$gleason_score <- clinical$`Radical Prostatectomy Gleason Score for Prostate Cancer`

clincal <- clinical[c("id","Disease_Free_Survival","Disease_Free_Status","age","pathologic_T","pathologic_N","psa_value","gleason_score")]
# Samples with clinical information of NA were removed
clincal <- clincal[!is.na(clincal$Disease_Free_Status ), ]
clincal$Disease_Free_Survival <- as.numeric(clincal$Disease_Free_Survival)

clincal <- unique(clincal)
# Load the prostate cancer expression data collected from TCGA official website 
load("D:/TCGA/HTSeq-FPKM/TCGA-PRAD.Rdata")
tcgaReplicateFilterId <- tcgaReplicateFilter(colnames(HTSeq_FPKM)[-1],analyte_target="RNA")
HTSeq_FPKM <- HTSeq_FPKM[,c("gene_id",tcgaReplicateFilterId)]
xx <- substr(colnames(HTSeq_FPKM)[-1],14,16)
table(xx)
HTSeq_FPKM[1:3,1:3]
colnames(HTSeq_FPKM) <- substr(colnames(HTSeq_FPKM),1,16)

NormalId <- colnames(HTSeq_FPKM)[-1][substr(colnames(HTSeq_FPKM)[-1],14,15) =="11"]
TumorId <- colnames(HTSeq_FPKM)[-1][substr(colnames(HTSeq_FPKM)[-1],14,15)=="01"]

coTumorId <- TumorId[substr(TumorId,1,12) %in% clincal$id]

# 489 tumor cases
str(clincal)
clincal <- clincal[clincal$id %in%  substr(coTumorId,1,12),]
write.table(clincal,file = "01rawData/TCGA_PRAD_clinical.txt",row.names = F,quote = F,sep = "\t")
#There were 489 tumor cases after the duplication was deleted

HTSeq_FPKM <- HTSeq_FPKM[,c("gene_id",NormalId,coTumorId)] %>% column_to_rownames("gene_id")
HTSeq_FPKM[1:3,1:5]
# fpkm To Tpm
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
HTSeq_TPM <- apply(HTSeq_FPKM,2,fpkmToTpm)
HTSeq_TPM <- HTSeq_TPM %>%as.data.frame()%>% rownames_to_column("gene_id")
HTSeq_TPM[1:3,1:5]
colSums(HTSeq_TPM[-1])
# Load gene annotation
load("D:/TCGA/gencode.v22.annotation.Rdata")
gtf_df[1:3,1:3]
gtf_df <- gtf_df[gtf_df$"gene_type" == "protein_coding",]
gtf_df$"gene_id" <- substr(gtf_df$"gene_id", 1,15)



exprSet <- HTSeq_TPM %>%
  dplyr::inner_join(gtf_df,by="gene_id") %>%

  dplyr::select(-c(gene_id,gene_type))%>%
 
  dplyr::select(gene_name,everything()) %>%
  
  dplyr::group_by(gene_name) %>%
  dplyr::summarise_all(mean)
exprSet[1:3,1:4]
write.table(exprSet,file = "01rawData/TCGA-PRAD_TPM.txt",row.names = F,quote = F,sep = "\t")
save(exprSet,file = "01rawData/TCGA-PRAD_TPM.Rdata")
load("01rawData/TCGA-PRAD_TPM.Rdata")
datExpr <- exprSet %>% column_to_rownames("gene_name")

group_list=c(rep("Normal",length(NormalId)),
             rep("Tumor",length(coTumorId)))  

group_list=factor(group_list)

group_list <- relevel(group_list, ref="Normal")
table(group_list)
# group_list
# Normal  Tumor 
# 52    489 

exprSet <- datExpr[,c(NormalId,coTumorId)]

# Low expression genes were filtered out
exprSet=exprSet[rowMeans(exprSet)>0.001,]
exprSet1 <- log2(exprSet+1)

library(limma)

# Differential gene analysis was performed using limma
design=model.matrix(~ group_list)
fit=lmFit(exprSet1,design)
fit=eBayes(fit)
res=topTable(fit,adjust='fdr',coef="group_listTumor",number=Inf)
allDiff <- na.omit(res)

logFCCutoff <- 1
pvalueCutoff <- 0.05

outDiff=allDiff[(abs(allDiff$logFC)>logFCCutoff & allDiff$adj.P.Val<pvalueCutoff),]
outDiff <- outDiff %>% rownames_to_column(var = "id")
write.table(outDiff,file="02DEGs/TCGA-PRAD_limma_diff.txt",row.names=F,quote=F,sep = "\t")

allDiff <- allDiff  %>% rownames_to_column(var = "id")
write.table(allDiff, '02DEGs/TCGA-PRAD_limma_alldiff.txt', sep = '\t', row.names=F, quote = FALSE)

#Heat maps of differential genes were created using Pheatmap  
library(pheatmap)
geneNum=50
outDiff=outDiff[order(as.numeric(as.vector(outDiff$logFC))),]
# diffGeneName=rownames( outDiff )
diffGeneName=outDiff$id 
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(geneNum*2)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp= exprSet[hmGene,]

hmExp = log2(hmExp+1)


Type= factor(c(rep("Normal",52),rep("Tumor",489)),levels = c("Normal","Tumor"))
names(Type)=colnames(hmExp)
Type=as.data.frame(Type)

# loc <- order(Type,colSums(hmExp),decreasing = T)
pdf(file="02DEGs/mRNA_heatmap.pdf",height=8,width=9)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8,
         border=FALSE)
dev.off()

#Draw a volcano map using ggplot2
library(ggplot2)
library(ggrepel)
res <- allDiff

Significant=ifelse((res$adj.P.Val< 0.05 & abs(res$logFC)> 1), ifelse(res$logFC > 1,"Up","Down"), "Not")

p = ggplot(res, aes(logFC, -log10(adj.P.Val)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("#004BFB", "#BFBFBF", "#F91F10"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-1,1), colour="black", linetype="dashed")+
  geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+xlab("log2FC")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+theme_bw()


#保存为图片
pdf("02DEGs/mRNA_vol.pdf",width=5.5,height=5)
print(p)
dev.off()

#Draw a Venn diagram using ggvenn  
# GSEA metabolism-related gene
library(clusterProfiler)
gmtFile <- "01rawData/c2.cp.kegg.v7.4.symbols.gmt"
geneSet=read.gmt(gmtFile)
metabolism <- geneSet[grepl("METABOLISM",geneSet$term),] 
write.table(metabolism,file = "01rawData/METABOLISM_geneSet.txt",row.names = F,quote = F,sep = "\t")
outDiff  <- read.table("02DEGs/TCGA-PRAD_limma_diff.txt",check.names = F,sep="\t",header = T)

x <- list(DEGs=outDiff$id,MGs=unique(metabolism$gene) )
library(ggvenn)
mycol <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#ADD1E5",
           '#CCEBC5', '#FCCDE5', '#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")
ggvenn(x,c("DEGs","MGs"),
       stroke_size = 0.3,
       show_percentage = F,
       fill_color =mycol[4:8])
ggsave("02DEGs/venn.pdf",width = 4.2,height = 4.5)

intersectGenes<- intersect(outDiff$id,unique(metabolism$gene))
write.table(intersectGenes,file = "02DEGs/venn_intersectGenes.txt",row.names = F,quote = F,sep = "\t",col.names = F)


# Draw heat maps of differential metabolic genes in normal and tumor samples
load("01rawData/TCGA-PRAD_TPM.Rdata")
exprSet0 <- exprSet %>% column_to_rownames("gene_name")
HYPOXIA <- read.table("02DEGs/venn_intersectGenes.txt",check.names = F,sep="\t",header = F)
hmExp <- exprSet0[rownames(exprSet0) %in% HYPOXIA$V1,]
Type= factor(c(rep("Normal",52),rep("Tumor",489)),levels = c("Normal","Tumor"))
names(Type)=colnames(hmExp)
Type=as.data.frame(Type)

# loc <- order(Type,colSums(hmExp),decreasing = T)
pdf(file="02DEGs/MGs_mRNA_heatmap.pdf",height=6,width=9)
pheatmap(log(hmExp+1) ,
         # hmExp,
         annotation=Type,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8,
         border=FALSE)
dev.off()

#Make input files for univariate COX analysis  
load("01rawData/TCGA-PRAD_TPM.Rdata")
exprSet <- exprSet %>% column_to_rownames("gene_name")

HYPOXIA <- read.table("02DEGs/venn_intersectGenes.txt",check.names = F,sep="\t",header = F)

clincal <- read.table("01rawData/TCGA_PRAD_clinical.txt",sep = "\t",header = T,check.names = F)
exprSet0 <- exprSet[ HYPOXIA$V1   ,c(coTumorId)]
exprSet0_out <- exprSet0 %>% rownames_to_column("id")
write.table(exprSet0_out,file = "03Cox/TCGA_PRAD_diff_MGs_exprSet.txt",row.names = F,quote = F,sep = "\t")

exprSet0 <- as.data.frame( t(exprSet0))
# exprSet$Id <- rownames(exprSet)
exprSet0$id <- substr(rownames(exprSet0),1,12)
exprSet_Time <- merge(clincal,exprSet0,by="id")


exprSet_Time[1:3,1:4]

write.table(exprSet_Time,file = "03Cox/TCGA_PRAD_diff_MGs_exprSetTime.txt",row.names = F,quote = F,sep = "\t")


# Heatmap of differential metabolic genes in NFM clustering
exprSet <- read.table("03Cox/TCGA_PRAD_diff_MGs_exprSet.txt",check.names = F,sep="\t",header = T)
exprSet0 <- exprSet %>% column_to_rownames("id")
Cluster_res <- read.table("03NMF/TCGA-PRAD_NMF_Cluster.txt",check.names = F,sep="\t",header = T)
Cluster1 <- Cluster_res$id[Cluster_res$Cluster =="Cluster1"]
Cluster2 <- Cluster_res$id[Cluster_res$Cluster =="Cluster2"]
hmExp <- exprSet0[,c(Cluster1,Cluster2)]
Type= factor(c(rep("Cluster1",264),rep("Cluster2",225)),levels = c("Cluster1","Cluster2"))
names(Type)=colnames(hmExp)
Type=as.data.frame(Type)

# loc <- order(Type,colSums(hmExp),decreasing = T)
pdf(file="03NMF/MGs_Cluster_heatmap.pdf",height=6,width=9)
pheatmap(log(hmExp+1) ,
         # hmExp,
         annotation=Type,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8,
         border=FALSE)
dev.off()
