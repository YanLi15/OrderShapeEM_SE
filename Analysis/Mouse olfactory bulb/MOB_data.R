###############################################################
#             Mouse olfactory bulb data analysis              #
###############################################################
library(SPARK)
library(dplyr)
setwd("E:/postdoc/research/spatial genetic structure/manuscript/OrderShapeEM_SE")
##-------------------------------------------------------------
## Load the spatial MOB data and apply SPARK
##-------------------------------------------------------------
# read the raw counts (spatial data)
counts_MOB <- read.table("./RealData/Rep11_MOB_count_matrix-1.tsv", check.names = F)
# extract the coordinates from the raw data
info_MOB <- cbind.data.frame(x = as.numeric(sapply(strsplit(rownames(counts_MOB), split = "x"), "[", 1)), 
                             y = as.numeric(sapply(strsplit(rownames(counts_MOB), split = "x"), "[", 2)))
rownames(info_MOB) <- rownames(counts_MOB)
# filter genes and cells/spots (excluding the gene that are lowly expressed)
spark_MOB <- CreateSPARKObject(counts = t(counts_MOB), location = info_MOB, 
                               percentage = 0.1, min_total_counts = 10)
# total counts for each cell/spot
spark_MOB@lib_size <- apply(spark_MOB@counts, 2, sum)

# Estimating Parameter Under Null with SPARK
spark_MOB <- spark.vc(spark_MOB, covariates = NULL, lib_size = spark_MOB@lib_size, 
                      num_core = 10, verbose = T, fit.maxiter = 500)
spark_MOB <- spark.test(spark_MOB, check_positive = T, verbose = T)

##----------------------------------------------------------------------------------
## Load the single-cell RNA data from MOB and make differential analysis with Seurat
##----------------------------------------------------------------------------------
library(loomR)
library(Seurat)
library(dplyr)

# Download l1_olfactory.loom from http://mousebrain.org/tissues.html
# Connect to the loom file in read/write mode
scRNA_MOB <- connect(filename = "./Covariates/l1_olfactory.loom", mode = "r+")
# Access the full data matrix, n rows(cells) * m columns(genes)
scRNA_MOB_count <- scRNA_MOB[["matrix"]][,]

# Access all gene names
MOB_genes <- scRNA_MOB$row.attrs$Gene[]
MOB_cellID <- scRNA_MOB$col.attrs$CellID[]

rm(scRNA_MOB)
colnames(scRNA_MOB_count) <- MOB_genes
rownames(scRNA_MOB_count) <- MOB_cellID
scRNA_MOB_count <- t(scRNA_MOB_count)

scRNA_MOB_seurat <- CreateSeuratObject(counts = scRNA_MOB_count, min.cells = 1, min.features = 50)

# QC
scRNA_MOB_seurat[["percent.mt"]] <- PercentageFeatureSet(scRNA_MOB_seurat, pattern = "^MT-")
scRNA_MOB_seurat <- subset(scRNA_MOB_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Nomalize the data
scRNA_MOB_seurat <- NormalizeData(scRNA_MOB_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly expressed genes
scRNA_MOB_seurat <- FindVariableFeatures(scRNA_MOB_seurat, selection.method = "vst", nfeatures = 2000)

# Scaling the data
sc.genes <- rownames(scRNA_MOB_seurat)
scRNA_MOB_seurat <- ScaleData(scRNA_MOB_seurat, features = sc.genes)

# PCA linear dimensional reduction
scRNA_MOB_seurat <- RunPCA(scRNA_MOB_seurat, features = VariableFeatures(object = scRNA_MOB_seurat))

# Determine the dimensionality of the dataset
scRNA_MOB_seurat <- JackStraw(scRNA_MOB_seurat, num.replicate = 100)
scRNA_MOB_seurat <- ScoreJackStraw(scRNA_MOB_seurat, dims = 1:20)
# plot
JackStrawPlot(scRNA_MOB_seurat, dims = 1:20)

# Cluster the cells
scRNA_MOB_seurat <- FindNeighbors(scRNA_MOB_seurat, dims = 1:20)
scRNA_MOB_seurat <- FindClusters(scRNA_MOB_seurat, resolution = 0.5)

# TSNE clustering
scRNA_MOB_seurat <- RunTSNE(scRNA_MOB_seurat, dims = 1:10)
DimPlot(scRNA_MOB_seurat, reduction = "tsne")
DimPlot(scRNA_MOB_seurat, reduction = "tsne", label = TRUE)

# compute p-values for every cluster compared to all remaining cells
scRNA_MOB <- FindAllMarkers(scRNA_MOB_seurat, return.thresh = 1)
sc_MOB <- scRNA_MOB %>% group_by(gene) %>% filter(p_val == min(p_val)) %>% dplyr::slice(1)
sc_MOB <- as.data.frame(sc_MOB)
rownames(sc_MOB) <- sc_MOB$gene

save(spark_MOB, scRNA_MOB, sc_MOB, file = "./Output/MOB.RData")