
################################################################
#                  Human heart data analysis                   #
################################################################

rm(list = ls())
load("./Output/HH.Rdata")
source("./R/RunAll.R")

alpha <- 0.05
# se.genes <- rownames(spark_HH@counts)[which(spark_HH@res_mtest[,"adjusted_pvalue"] < alpha)]

target.genes <- intersect(sc_HH$gene, rownames(spark_HH@counts))
# Primary pvalues
pval_spark <- spark_HH@res_mtest[target.genes,"combined_pvalue"]

# Auxiliary covariates
covar <- sc_HH[target.genes,"p_val"]
Rej.res <- RunAll(pval_spark, covar, alpha = alpha)

genes.by <- target.genes[Rej.res$rej.by]
genes.bh <- target.genes[Rej.res$rej.bh]
genes.sabha <- target.genes[Rej.res$rej.sabha]
genes.osem <- target.genes[Rej.res$rej.osem]

save(genes.by, genes.bh, genes.sabha, genes.osem, file = "./Output/HH_result.RData")

#--------------------------------------------------------------------------------
# Plot the results
#--------------------------------------------------------------------------------
load("./Output/HH_result.Rdata")
library(UpSetR)
listInput <- list(BY = genes.by, BH = genes.bh, SABHA = genes.sabha, OrderShapeEM = genes.osem)
pdf(file="Output/plots/HH_result.pdf",onefile = FALSE,width=9,height=6)
upset(fromList(listInput), sets.x.label = "Number of SE genes", 
      matrix.color = "darkblue", main.bar.color = "darkblue", 
      sets.bar.color = "darkblue", order.by = "freq")
dev.off()