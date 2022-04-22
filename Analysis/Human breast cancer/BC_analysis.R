###############################################################
#              Human breast cancer data analysis              #
###############################################################
rm(list = ls())
load("./Output/BC.RData")
source("./R/RunAll.R")
source("./R/Cauchy.R")

se.genes <- rownames(spark_BC@counts)[which(spark_BC@res_mtest[,"adjusted_pvalue"] < 0.05)]

##-----------------------------------------------------------------------------------
## Case 7: Multiple covariates assisted analysis (GWAS data, TCGA BRCA and TCGA THCA)
##-----------------------------------------------------------------------------------
aux.genes <- intersect(rownames(TCGA.br.result), rownames(TCGA.th.result))
aux.genes <- intersect(aux.genes, rownames(GWAS.result))
aux.Pvals <- rbind(TCGA.br.result[aux.genes,"pvalue"], 
             GWAS.result[aux.genes,"pvalue"], 
             TCGA.th.result[aux.genes,"pvalue"])
aux.pvalue7 <- Cauchy(aux.Pvals)
names(aux.pvalue7) <- aux.genes

target.genes <- intersect(aux.genes, rownames(spark_BC@counts))

pval_spark <- spark_BC@res_mtest[target.genes,"combined_pvalue"]
covar7 <- aux.pvalue7[target.genes]

Rej.res7 <- RunAll(pval_spark, covar7, alpha = 0.05)

genes.by <- target.genes[Rej.res7$rej.by]
genes.bh <- target.genes[Rej.res7$rej.bh]
genes.sabha7 <- target.genes[Rej.res7$rej.sabha]
genes.osem7 <- target.genes[Rej.res7$rej.osem]

##-----------------------------------------------------------------------------------

##-------------------------------------------------------------
## Case 1: Single Covariate-assisted analysis (TCGA-BRCA)
##-------------------------------------------------------------
covar1 <- TCGA.br.result[target.genes, "pvalue"]

Rej.res1 <- RunAll(pval_spark, covar1, alpha = 0.05)

genes.sabha1 <- target.genes[Rej.res1$rej.sabha]
genes.osem1 <- target.genes[Rej.res1$rej.osem]

##-------------------------------------------------------------
## Case 2: Single Covariate-assisted analysis (TCGA-THCA)
##-------------------------------------------------------------
covar2 <- TCGA.th.result[target.genes, "pvalue"]

Rej.res2 <- RunAll(pval_spark, covar2, alpha = 0.05)

genes.sabha2 <- target.genes[Rej.res2$rej.sabha]
genes.osem2 <- target.genes[Rej.res2$rej.osem]

##-------------------------------------------------------------
## Case 3: Single Covariate-assisted analysis (GWAS)
##-------------------------------------------------------------
covar3 <- GWAS.result[target.genes, "pvalue"]

Rej.res3 <- RunAll(pval_spark, covar3, alpha = 0.05)

genes.sabha3 <- target.genes[Rej.res3$rej.sabha]
genes.osem3 <- target.genes[Rej.res3$rej.osem]

##--------------------------------------------------------------------------------------------
## Case 4: Multiple covariates assisted analysis (GWAS data and TCGA-BRCA RNA sequencing data)
##--------------------------------------------------------------------------------------------
aux.genes <- intersect(rownames(TCGA.br.result), rownames(GWAS.result))
aux.Pvals <- rbind(TCGA.br.result[aux.genes,"pvalue"], 
                   GWAS.result[aux.genes,"pvalue"])
aux.pvalue4 <- Cauchy(aux.Pvals)
names(aux.pvalue4) <- aux.genes

covar4 <- aux.pvalue4[target.genes]

Rej.res4 <- RunAll(pval_spark, covar4, alpha = 0.05)

genes.sabha4 <- target.genes[Rej.res4$rej.sabha]
genes.osem4 <- target.genes[Rej.res4$rej.osem]

##----------------------------------------------------------------------------------
## Case 5: Multiple covariates assisted analysis (TCGA-BRCA data and TCGA-THCA data)
##----------------------------------------------------------------------------------
aux.genes <- intersect(rownames(TCGA.br.result), rownames(TCGA.th.result))
aux.genes <- aux.genes[complete.cases(aux.genes)]
aux.Pvals <- rbind(TCGA.br.result[aux.genes,"pvalue"],
                   TCGA.th.result[aux.genes,"pvalue"])
aux.pvalue5 <- Cauchy(aux.Pvals)
names(aux.pvalue5) <- aux.genes

covar5 <- aux.pvalue5[target.genes]

Rej.res5 <- RunAll(pval_spark, covar5, alpha = 0.05)

genes.sabha5 <- target.genes[Rej.res5$rej.sabha]
genes.osem5 <- target.genes[Rej.res5$rej.osem]

##-------------------------------------------------------------------------------
## Case 6: Multiple covariates assisted analysis (GWAS data and TCGA-THCA data)
##-------------------------------------------------------------------------------
aux.genes <- intersect(rownames(TCGA.th.result), rownames(GWAS.result))
aux.Pvals <- rbind(GWAS.result[aux.genes,"pvalue"], 
                   TCGA.th.result[aux.genes,"pvalue"])
aux.pvalue6 <- Cauchy(aux.Pvals)
names(aux.pvalue6) <- aux.genes

covar6 <- aux.pvalue6[target.genes]

Rej.res6 <- RunAll(pval_spark, covar6, alpha = 0.05)

genes.sabha6 <- target.genes[Rej.res6$rej.sabha]
genes.osem6 <- target.genes[Rej.res6$rej.osem]

##-----------------------------------------------------------------------------------

save.image("./Output/BC_result.RData")


#--------------------------------------------------------------------------------
# Plot the results
#--------------------------------------------------------------------------------
load("./Output/BC_result.Rdata")
library(UpSetR)

listInput.all <- list(BY = genes.by7, BH = genes.bh7, SABHA.Case7 = genes.sabha7, 
                      OrderShapeEM.Case1 = genes.osem1, OrderShapeEM.Case2 = genes.osem2,
                      OrderShapeEM.Case3 = genes.osem3, OrderShapeEM.Case4 = genes.osem4,
                      OrderShapeEM.Case5 = genes.osem5, OrderShapeEM.Case6 = genes.osem6,
                      OrderShapeEM.Case7 = genes.osem7)
pdf(file="Output/plots/BC_result.pdf",onefile = FALSE,width=9,height=6)
upset(fromList(listInput.all), nsets = 10, sets.x.label = "Number of SE genes", 
      matrix.color = "darkblue", main.bar.color = "darkblue", 
      sets.bar.color = "darkblue", order.by = "freq")
dev.off()