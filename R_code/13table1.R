rm(list=ls())

#install.packages("table1")
library(table1)
Cluster_res <- read.table("03NMF/TCGA-PRAD_NMF_Cluster.txt",header = T,sep = "\t",check.names = F)
Cluster_res <- Cluster_res[,c(1,ncol(Cluster_res))]
Cluster_res$id <- substr(Cluster_res$id ,1,12)
clincal <- read.table("01rawData/TCGA_PRAD_clinical.txt",sep = "\t",header = T,check.names = F)
dat <- merge(clincal,Cluster_res,by="id")
dat <- dat[,-c(1:3)]

dat$age <- ifelse(dat$age >=60,">=60","<60")
dat$psa_value<- ifelse(is.na(dat$psa_value) ,NA,ifelse(dat$psa_value >=4,">=4","<4"))
dat$gleason_score <- ifelse(dat$gleason_score >=8,">=8","<8")

label(dat$pathologic_T)   <- "pathologic_T"
# label(dat$`expression-based subtype`)   <- "Expression-based Subtype"
label(dat$pathologic_N)    <- "pathologic_N"
label(dat$gleason_score) <- "gleason_score"
label(dat$Cluster)   <- "Cluster"

head(dat)

# The independent Chi-square test is used for categorical variables and the T-test is used for continuous variables (other tests may be used if desired)  
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
table_output <- table1(~ age  + pathologic_T   + pathologic_N +psa_value +gleason_score  | Cluster, data=dat, overall = F,
       topclass="Rtable1-zebra",extra.col=list(`P-value`=pvalue))
# table1(~  age + gender + `who class`| Cluster, data=dat, overall = F,
#        topclass="Rtable1-zebra",extra.col=list(`P-value`=pvalue))
table_output <- print(table_output)
write.csv(table_output,file = "01rawData/table_output.csv")
