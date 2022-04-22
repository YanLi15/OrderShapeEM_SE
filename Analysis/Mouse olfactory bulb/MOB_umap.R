library(umap)
library(amap)
source("./R/PlotFuncs.R")
load("./Output/MOB.RData")
load("./Output/MOB_result.RData")

# anscombe variance stabilizing transformation: NB
vst_count <- var_stabilize(spark_MOB@counts) # R function in funcs.R

# the SE genes based on an FDR cutoff of 0.05
sig_vst_count <- vst_count[genes.osem, ]
sig_vst_res <- t(apply(sig_vst_count, 1, LMReg, T = log(spark_MOB@lib_size)))

hc <- hcluster(sig_vst_res, method = "euc", link = "ward", nbproc = 1, 
               doubleprecision = TRUE)
numC <- 3
memb <- cutree(hc, k = numC)

counts <- spark_MOB@counts[genes.osem,]
scaled_counts <- t(scale(t(counts)))
umap_results <- umap::umap(scaled_counts)
labels <- memb
labels[which(labels == 1)] <- "I"
labels[which(labels == 2)] <- "II"
labels[which(labels == 3)] <- "III"

df <- data.frame(umap_results$layout[,1], umap_results$layout[,2], labels)
colnames(df) <- c("UMAP1","UMAP2", "Cluster")
ggplot(df, aes(x = UMAP1, y = UMAP2, color = Cluster), color = labels) + geom_point()


