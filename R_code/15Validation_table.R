rm(list=ls())

#install.packages("table1")
library(table1)
Cluster_res <- read.table("03NMF/GSE70768_NMF_Cluster.txt",header = T,sep = "\t",check.names = F)
Cluster_res <- Cluster_res[,c(1,ncol(Cluster_res))]
# Cluster_res$id <- substr(Cluster_res$id ,1,12)
clincal <- read.table("01rawData/GSE70768_clinical.txt",sep = "\t",header = T,check.names = F)
colnames(clincal) <-  c("id","age","Clinical_stage","Pathology_stage","psa_at_diag","total_follow_up","tumour_gleason")
dat <- merge(clincal,Cluster_res,by="id")
dat <- dat[,-c(6)]
dat$tumour_gleason  <-  as.numeric( substr(dat$tumour_gleason,0,1))

dat$psa_value<- ifelse(is.na(dat$psa_at_diag) ,NA,ifelse(dat$psa_at_diag >=4,">=4","<4"))
dat$gleason_score <- ifelse(dat$tumour_gleason >=8,">=8","<8")
dat$Clinical_stage <-  ( substr(dat$Clinical_stage,1,2))
dat$Pathology_stage <-  ( substr(dat$Pathology_stage,1,2))

label(dat$Clinical_stage)   <- "Clinical_stage"
# label(dat$`expression-based subtype`)   <- "Expression-based Subtype"
label(dat$Pathology_stage)    <- "Pathology_stage"
label(dat$gleason_score) <- "Gleason_score"
label(dat$Cluster)   <- "Cluster"


head(dat)

pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}
# table_output <- table1(~ `clinical stage` + morphology  + `% intratumoral tils`+`% stromal tils`| risk, data=dat, overall = F,
#        topclass="Rtable1-zebra")
table_output <- table1(~  Clinical_stage   + Pathology_stage +psa_value +gleason_score  | Cluster, data=dat, overall = F,
       topclass="Rtable1-zebra",extra.col=list(`P-value`=pvalue))
# table1(~  age + gender + `who class`| Cluster, data=dat, overall = F,
#        topclass="Rtable1-zebra",extra.col=list(`P-value`=pvalue))
table_output <- print(table_output)
write.csv(table_output,file = "01rawData/GSE70768_table_output.csv")



Cluster_res <- read.table("03NMF/GSE70769_NMF_Cluster.txt",header = T,sep = "\t",check.names = F)
Cluster_res <- Cluster_res[,c(1,ncol(Cluster_res))]
# Cluster_res$id <- substr(Cluster_res$id ,1,12)
clincal <- read.table("01rawData/GSE70769_clinical.txt",sep = "\t",header = T,check.names = F)
colnames(clincal) <-  c("id","Clinical_stage","Pathology_stage","psa_at_diag","total_follow_up","tumour_gleason")
dat <- merge(clincal,Cluster_res,by="id")
dat <- dat[,-c(5)]
dat$tumour_gleason  <-  as.numeric( substr(dat$tumour_gleason,0,1))

dat$psa_value<- ifelse(is.na(dat$psa_at_diag) ,NA,ifelse(dat$psa_at_diag >=4,">=4","<4"))
dat$gleason_score <- ifelse(dat$tumour_gleason >=8,">=8","<8")
dat$Clinical_stage <-  ( substr(dat$Clinical_stage,1,2))
dat$Pathology_stage <-  ( substr(dat$Pathology_stage,1,2))

label(dat$Clinical_stage)   <- "Clinical_stage"
# label(dat$`expression-based subtype`)   <- "Expression-based Subtype"
label(dat$Pathology_stage)    <- "Pathology_stage"
label(dat$gleason_score) <- "Gleason_score"
label(dat$Cluster)   <- "Cluster"


head(dat)

pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}
# table_output <- table1(~ `clinical stage` + morphology  + `% intratumoral tils`+`% stromal tils`| risk, data=dat, overall = F,
#        topclass="Rtable1-zebra")
table_output <- table1(~  Clinical_stage   + Pathology_stage +psa_value +gleason_score  | Cluster, data=dat, overall = F,
                       topclass="Rtable1-zebra",extra.col=list(`P-value`=pvalue))
# table1(~  age + gender + `who class`| Cluster, data=dat, overall = F,
#        topclass="Rtable1-zebra",extra.col=list(`P-value`=pvalue))
table_output <- print(table_output)
write.csv(table_output,file = "01rawData/GSE70769_table_output.csv")
