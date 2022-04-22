
################################################################
#              Mouse olfactory bulb data analysis              #
################################################################

rm(list = ls())
load("./Output/MOB.Rdata")
source("./R/RunAll.R")

alpha <- 0.05
se.genes <- rownames(spark_MOB@counts)[which(spark_MOB@res_mtest[,"adjusted_pvalue"] < alpha)]

target.genes <- intersect(sc_MOB$gene, rownames(spark_MOB@counts))
# Primary pvalues
pval_spark <- spark_MOB@res_mtest[target.genes,"combined_pvalue"]

# Auxiliary covariates
covar <- sc_MOB[target.genes,"p_val"]
Rej.res <- RunAll(pval_spark, covar, alpha = alpha)

genes.by <- target.genes[Rej.res$rej.by]
genes.bh <- target.genes[Rej.res$rej.bh]
genes.sabha <- target.genes[Rej.res$rej.sabha]
genes.osem <- target.genes[Rej.res$rej.osem]



save(se.genes, genes.by, genes.bh, genes.sabha, genes.osem, file = "./Output/MOB_result.RData")

#--------------------------------------------------------------------------------
# Plot the results
#--------------------------------------------------------------------------------
load("./Output/MOB_result.Rdata")
library(UpSetR)
listInput <- list(BY = genes.by, BH = genes.bh, SABHA = genes.sabha, OrderShapeEM = genes.osem)
pdf(file="Output/plots/MOB_result.pdf",onefile = FALSE,width=9,height=6)
upset(fromList(listInput), sets.x.label = "Number of SE genes", 
      matrix.color = "darkblue", main.bar.color = "darkblue", 
      sets.bar.color = "darkblue", order.by = "freq")
dev.off()