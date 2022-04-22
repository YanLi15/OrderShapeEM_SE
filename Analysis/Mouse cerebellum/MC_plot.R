
##------------------------------------------------------------------------
## Summarize Pattern using Hclust based on the results of OrderShapeEM
##------------------------------------------------------------------------

rm(list = ls())

source("./R/PlotFuncs.R")
library(amap)

load("./Output/MC.RData")
load("./Output/MC_result.RData")
## SPARK method (all genes)
# anscombe variance stabilizing transformation: NB
vst_count <- var_stabilize(spark_MC@counts) # R function in funcs.R

# the SE genes based on an FDR cutoff of 0.05
sig_vst_count <- vst_count[genes.osem, ]
sig_vst_res <- t(apply(sig_vst_count, 1, LMReg, T = log(spark_MC@lib_size)))

hc <- hcluster(sig_vst_res, method = "euc", link = "ward", nbproc = 1, 
               doubleprecision = TRUE)
numC <- 2
memb <- cutree(hc, k = numC)

# The mean residuals of the three patterns for each location
cent <- NULL
for (k in 1:numC) {
  cent <- cbind(cent, colMeans(sig_vst_res[memb == k, , drop = FALSE]))
}

position_cord <- spark_MC@location
rownames(position_cord) <- rownames(cent)

# The relative residuals
rel_cent <- t(apply(cent, 1, relative_func))
pd <- setNames(cbind.data.frame(position_cord, rel_cent), c("x", "y", paste0("Pattern ", c("II", "I"))))
# pd2 <- setNames(cbind.data.frame(position_cord, rel_cent), c("x", "y", paste0("Pattern ", c("III", "II", "I"))))
MBP <- lapply(1:numC, function(x) {
  pattern_plot2(pd, x, xy = T, main = T, titlesize = 1.5)
})

grid.arrange(grobs = MBP[numC:1], nrow = numC)
