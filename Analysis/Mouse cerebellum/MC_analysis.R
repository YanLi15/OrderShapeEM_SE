
################################################################
#              Mouse olfactory bulb data analysis              #
################################################################

rm(list = ls())
load("./Output/MC.Rdata")
source("./R/RunAll.R")

se.genes <- rownames(spark_MC@counts)[which(spark_MC@res_mtest[,"adjusted_pvalue"] < 0.05)]

target.genes <- intersect(rownames(sc_MC), rownames(spark_MC@counts))

# Primary pvalues
pval_spark <- spark_MC@res_mtest[target.genes,"combined_pvalue"]
# Auxiliary covariates
covar <- sc_MC[target.genes,"p_val"]
Rej.res <- RunAll(pval_spark, covar, alpha = 0.05)

genes.by <- target.genes[Rej.res$rej.by]
genes.bh <- target.genes[Rej.res$rej.bh]
genes.sabha <- target.genes[Rej.res$rej.sabha]
genes.osem <- target.genes[Rej.res$rej.osem]

save(se.genes, genes.by, genes.bh, genes.sabha, genes.osem, file = "./Output/MC_result.RData")

# write.csv(genes.by, file = "./Output/by.csv")
# write.csv(genes.bh, file = "./Output/bh.csv")
# write.csv(genes.sabha, file = "./Output/sabha.csv")
# write.csv(genes.osem, file = "./Output/osem.csv")

#--------------------------------------------------------------------------------
# Plot the results
#--------------------------------------------------------------------------------
load("./Output/MC_result.Rdata")
library(UpSetR)
listInput <- list(BY = genes.by, BH = genes.bh, SABHA = genes.sabha, OrderShapeEM = genes.osem)
pdf(file="Output/plots/MC_result.pdf",onefile = FALSE,width=9,height=6)
upset(fromList(listInput), sets.x.label = "Number of SE genes", 
      matrix.color = "darkblue", main.bar.color = "darkblue", 
      sets.bar.color = "darkblue", order.by = "freq")
dev.off()