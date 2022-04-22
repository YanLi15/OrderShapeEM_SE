rm(list = ls())
library(spdep)
load("./Output/MOB.RData")
load("./Output/MOB_result.RData")

counts = spark_MOB@counts
coords = as.matrix(spark_MOB@location)
# sp = SpatialPoints(spark_MOB@location)
k1 = knn2nb(knearneigh(coords))
all.linked <- max(unlist(nbdists(k1, coords)))
# Define neighboring polygons
col.nb <- dnearneigh(coords, 0, all.linked, 
                     row.names=row.names(coords), 
                     longlat = FALSE)
# Assign weights to the neighbors
# col.W <- nb2listw(col.nb, style = "W")
col.S <- nb2listw(col.nb, style = "S")
# Computing the Moran’s I statistic
# MI <- apply(counts, 1, moran, listw = col.W, n = length(col.nb), S0 = Szero(col.W))
MI.S <- apply(counts, 1, moran, listw = col.S, n = length(col.nb), S0 = Szero(col.S))
# Performing a hypothesis test
# MI.test <- apply(counts, 1, moran.test, listw = col.W)
f <- function(x){x = x[[1]]}
MI.S <- lapply(MI.S, f)
MI <- NULL
for(i in 1:length(MI.S)){
  MI <- c(MI, MI.S[[i]])
}
names(MI) <- row.names(spark_MOB@counts)

target.genes <- intersect(sc_MOB$gene, rownames(spark_MOB@counts))
# Moran's I statistics of the SE genes identified by OrderShapeEM and SPARK(BY)
mi.by <- MI[genes.by]
mi.osem <- MI[genes.osem]
mi.all <- MI[target.genes]


factor <- factor(rep(c("OrderShapeEM", "SPARK","All"), 
                     times = c(length(genes.osem), length(genes.by), length(mi.all))))
dataset <- data.frame(value = c(mi.osem, mi.by, mi.all), group = factor)
boxplot(value ~ group, dataset, col = c("antiquewhite2", "coral", "cornflowerblue"), 
        xlab = "Methods", ylab = "Moran's I statistics", outline = FALSE, cex.lab = 1.2)