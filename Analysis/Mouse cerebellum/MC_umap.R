library(umap)
library(amap)
source("./R/PlotFuncs.R")
load("./Output/MC.RData")
load("./Output/MC_result.RData")

# anscombe variance stabilizing transformation: NB
vst_count <- var_stabilize(spark_MC@counts) # R function in funcs.R

# the SE genes based on an FDR cutoff of 0.05
sig_vst_count <- vst_count[genes.osem, ]
sig_vst_res <- t(apply(sig_vst_count, 1, LMReg, T = log(spark_MC@lib_size)))

hc <- hcluster(sig_vst_res, method = "euc", link = "ward", nbproc = 1, 
               doubleprecision = TRUE)
numC <- 2
memb <- cutree(hc, k = numC)

counts <- spark_MC@counts[genes.osem,]
# normalized_counts <- t(apply(counts, 1, min_max_norm))
# scaled_counts <- t(apply(counts, 1, stand))
scaled_counts <- t(scale(t(as.matrix(counts))))
umap_results <- umap::umap(scaled_counts)
labels <- memb
labels[which(labels == 1)] <- "I"
labels[which(labels == 2)] <- "II"

df <- data.frame(umap_results$layout[,1], umap_results$layout[,2], labels)
colnames(df) <- c("UMAP1","UMAP2", "Cluster")
ggplot(df, aes(x = UMAP1, y = UMAP2, color = Cluster), color = labels) + geom_point()


