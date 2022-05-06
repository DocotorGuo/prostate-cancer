# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
# BiocManager::install("ComplexHeatmap")
 
# BiocManager::install("GSVA")

library(GSVA)
library(ComplexHeatmap)
library(dplyr)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

# read marker genes of immune cells
immunity <- read.csv("04immune/immunitygene.csv", header = T)

# Remove celltypes that are not immune cells  
immunity <- immunity[!immunity$CellType %in% c("Blood vessels", "Normal mucosa", "SW480 cancer cells", "Lymph vessels"), ] %>% tidyr::separate_rows(Gene.Symbol,sep= " /// ")


immunity <- immunity %>% 
  split(., .$CellType) %>% 
  lapply(., function(x)(x$Symbol))
immunity <- lapply(immunity, unique)

# The infiltration level was quantified using ssGSEA
mycol <- c('#B3DE69', '#FDB462','#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7',"#999999")

library("limma")                                        
expFile="01rawData/TCGA-PRAD_TPM.txt" 
expr=read.table(expFile,sep="\t",header=T,check.names=F,row.names = 1)
xx <- substr(colnames(expr),14,16)
table(xx)
TumorId <- colnames(expr)[substr(colnames(expr),14,15)=="01"]
expr<- expr[,TumorId]

expr[1:3,1:5]
class(expr)
expr1 <- as.matrix(expr)

# ssGSEA, immune infiltration level was quantified using gsva function in GSVA package 
gsva1 <- as.data.frame(t(gsva(expr1, immunity, method = "ssgsea")))

ssgseaOut <- gsva1
ssgseaOut_df=cbind(id=rownames(ssgseaOut),ssgseaOut)
ssgseaOut_df[1:3,1:6]
write.table(ssgseaOut_df,file="04immune/TCGA-PRAD_ssGSEA_score.txt",sep="\t",quote=F,row.names=F)


gsva <- t(ssgseaOut)

library(pheatmap)
Cluster_res <- read.table("03NMF/TCGA-PRAD_NMF_Cluster.txt",header = T,sep = "\t",check.names = F)
Cluster1 <- Cluster_res$id[Cluster_res$Cluster =="Cluster1"]
Cluster2 <- Cluster_res$id[Cluster_res$Cluster =="Cluster2"]

gsva <- gsva[,c(Cluster1,Cluster2)]

estimate <- read.table("04immune/TCGA-PRAD_estimateScores.txt",header = T,sep = "\t",check.names = F,row.names = 1)
estimate <- estimate[c(Cluster1,Cluster2),]
annotation_col = data.frame(
  Type = factor(rep(c("Cluster1","Cluster2"), c(length(Cluster1) ,length(Cluster2))),levels = c("Cluster1","Cluster2")),
  StromalScore =estimate$StromalScore,
  ImmuneScore =estimate$ImmuneScore,
  ESTIMATEScore =estimate$ESTIMATEScore
  
)
rownames(annotation_col) = colnames(gsva)
head(annotation_col)
# A color list of custom comment information "#E69F00","#56B4E9"
ann_colors = list(
  # Type = c(`Cluster1` ="#028846", `Cluster2` = "red" )
  Type = c(Cluster1 = mycol[1], Cluster2 = mycol[2] ) )
head(ann_colors)

gsva <- gsva[apply(gsva, 1, function(x) sd(x)!=0),]
loc <- order(annotation_col$Type,colSums(gsva),decreasing = T)

pdf(file="04immune/TCGA-PRAD_ssGSEA_score_heatmap.pdf",height=8,width=10)
pheatmap(gsva,
         # gsva[,loc],
         scale = "row", 
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = F, 
         cluster_cols =F,
         border=FALSE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

dev.off()


# Boxplots
boxdata <-  as.data.frame(t(gsva))
risk <- read.table("03NMF/TCGA-PRAD_NMF_Cluster.txt",header = T,sep = "\t",check.names = F)
boxdata$id <- rownames(boxdata)
boxdata <- merge(boxdata,risk[,c("id","Cluster")],by="id")


write.table(boxdata, "04immune/TCGA-PRAD_ssGSEA_score_Type.txt", quote = F, row.names = F,sep = "\t")


boxdata1 <- boxdata[,-1]
library(tidyr)
library(ggpubr)
boxdata1 <- gather(boxdata1,cellType, value,colnames(boxdata1)[-ncol(boxdata1)])
boxdata1$Cluster <- factor(boxdata1$Cluster,levels = c("Cluster1","Cluster2"))

p2 <- ggboxplot(boxdata1, x = "cellType", y = "value",
               color = "Cluster", palette = mycol[c(1,2)],
               add = "jitter",
               add.params = list(alpha=0.6),)+ rotate_x_text(45) +xlab("")+ylab("ssGSEA score")
p2 + stat_compare_means(aes(group = Cluster),method="wilcox.test", label = "p.signif")

ggsave("04immune/TCGA-PRAD_ssGSEA_score_boxplot.pdf",height=5,width=10)