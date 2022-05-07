# Characterization of metabolism-associated molecular patterns in prostate cancer

    c2.cp.kegg.v7.4.symbols.gmt which downloaded from Molecular Signatures Database (MSigDB, https://www.gsea-msigdb.org/gsea/msigdb/). 
    In the present study, The transcriptomic expression profiles and corresponding clinical data of 489 prostate cancer tissues and 52 non-tumor tissues were obtained from the TCGA-PARD dataset in The Cancer Genome Atlas (TCGA, https://portal.gdc.cancer.gov/) and cBioPortal database (https://www.cbioportal.org/).GSE70768 and GSE70769 datasets were used as a testing set.

    the differentially expressed genes (DEGs) between prostate cancer tissues and non-tumor tissues were screened using the Limma package in R with the criteria of absolute log2 (fold change, FC) > 1 and adjusted P-value < 0.05, and the results visualized using ggplot package in R. Then, the differentially expressed metabolism-associated genes (MAGs) were obtained by intersecting the DEGs and MAGs. And the differentially expressed MAGs were visualized using pheatmap package in R. 
    
    Based on the differentially expressed MAGs, the 489 prostate cancer samples from the TCGA-PARD dataset were classed into different molecular subclusters using unsupervised non-negative matrix factorization (NMF) clustering via the NMF R package.Then, t-distributed stochastic neighbor embedding (t-SNE) was performed to verify the classification performance using the mRNA expression data of DE-MAGs. Kaplan-Meier (KM) DFS curves were drawn using the survival R package 
    
    Estimate R package was performed to evaluate the ESTIMATE score, immune score, and stromal score of each sample.Besides, a single-sample Gene set expression analysis (ssGSEA) was conducted to estimate the immune infiltration by calculating the enrich score of each gene in a special immune cell marker gene set based on the mRNA expression data. SsGSEA was performed using the GSVA package in R.GSEA was performed to explore the potential molecular mechanism between molecular subclusters.
    
    
