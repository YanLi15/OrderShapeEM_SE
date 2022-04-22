##-------------------------------------------------------------
## Spatial Distribution of Representative Genes 
##-------------------------------------------------------------
source("./R/PlotFuncs.R")
library(amap)
load("./Output/MOB.RData")
load("./Output/MOB_result.RData")

genes.overlap    <- intersect(genes.osem, genes.by) # overlap between OrderShapeEM and BY
genes.osem.only  <- genes.osem[!genes.osem%in%genes.overlap]
# genes.bh.only <- genes.bh[!genes.bh%in%genes.overlap]

osem.only.pval <- spark_MOB@res_mtest[genes.osem.only,"combined_pvalue"]
names(osem.only.pval) <- genes.osem.only
osem.only.pval <- sort(osem.only.pval)

gene_plot <- c("Elmod1", "Srsf3", "Col3a1", "Scg2", "Arhgef7", "Gsta4")
vst_ct <- var_stabilize(spark_MOB@counts) # R function in funcs.R
sig_vst_ct <- vst_ct[gene_plot, ]
rel_vst_ct <- apply(sig_vst_ct, 1, relative_func)


pltdat <- cbind.data.frame(spark_MOB@location[,1:2],rel_vst_ct)
genetitle <- c(expression("Elmod1 "),
               expression("Srsf3 "),
               expression("Col3a1 "),
               expression("Scg2 "),
               expression("Arhgef7"),
               expression("Gsta4"))

pp <- lapply(1:(ncol(pltdat)-2),
             function(x){pattern_plot2(pltdat,x,main=T,titlesize=1.5,title=genetitle[x])})
grid.arrange(grobs=pp, ncol=3)
