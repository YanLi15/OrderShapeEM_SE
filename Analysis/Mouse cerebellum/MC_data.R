rm(list=ls())
library(data.table)
library(Matrix)
library(SPARK)

##------------------------------------------------------------------------------------
## Load the spatial transcriptomic data and make SE analysis with the original SPARK
##------------------------------------------------------------------------------------
load("./RealData/SlideseqV2_ROI.rds")

# remove mt genes
mt_idx   <- grep("mt-",rownames(sp_count))
sp_count <- sp_count[-mt_idx,]

# remove cells without any expression
removed_cell <- which(as.vector(sp_sums_Rcpp(sp_count))==0)
if(length(removed_cell)>0){
  sp_count    <- sp_count[,-removed_cell]
  location    <- location[-removed_cell,]
}

rownames(location) <- colnames(sp_count)
location = as.data.frame(location)
spark_MC <- CreateSPARKObject(counts = sp_count, location = location, 
                           percentage = 0.1, min_total_counts = 10)
spark_MC@lib_size <- apply(spark_MC@counts, 2, sum)
spark_MC <- spark.vc(spark_MC, covariates = NULL, lib_size = spark_MC@lib_size, 
                  num_core = 5, verbose = T, fit.maxiter = 500)
spark_MC  <- spark.test(spark_MC, check_positive = T, verbose = T)

##-----------------------------------------------------------------------------------
## Load the single-cell RNA data from MC and make differential analysis with Seurat
##-----------------------------------------------------------------------------------
library(loomR)
library(Seurat)
library(dplyr)

# Download l1_cerebellum.loom from http://mousebrain.org/tissues.html
# Connect to the loom file in read/write mode
scRNA_MC <- connect(filename = "./Covariates/l1_cerebellum.loom", mode = "r+", skip.validate = TRUE)
# Access the full data matrix, n rows(cells) * m columns(genes)
scRNA_MC_count <- scRNA_MC[["matrix"]][,]

# Access all gene names
MC_genes <- scRNA_MC$row.attrs$Gene[]
MC_cellID <- scRNA_MC$col.attrs$CellID[]

rm(scRNA_MC)
colnames(scRNA_MC_count) <- MC_genes
rownames(scRNA_MC_count) <- MC_cellID
scRNA_MC_count <- t(scRNA_MC_count)

scRNA_MC_seurat <- CreateSeuratObject(counts = scRNA_MC_count, min.cells = 1, min.features = 50)

# QC
scRNA_MC_seurat[["percent.mt"]] <- PercentageFeatureSet(scRNA_MC_seurat, pattern = "^MT-")
scRNA_MC_seurat <- subset(scRNA_MC_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Nomalize the data
scRNA_MC_seurat <- NormalizeData(scRNA_MC_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly expressed genes
scRNA_MC_seurat <- FindVariableFeatures(scRNA_MC_seurat, selection.method = "vst", nfeatures = 2000)

# Scaling the data
sc.genes <- rownames(scRNA_MC_seurat)
scRNA_MC_seurat <- ScaleData(scRNA_MC_seurat, features = sc.genes)

# PCA linear dimensional reduction
scRNA_MC_seurat <- RunPCA(scRNA_MC_seurat, features = VariableFeatures(object = scRNA_MC_seurat))

# Determine the dimensionality of the dataset
scRNA_MC_seurat <- JackStraw(scRNA_MC_seurat, num.replicate = 100)
scRNA_MC_seurat <- ScoreJackStraw(scRNA_MC_seurat, dims = 1:20)
# plot
JackStrawPlot(scRNA_MC_seurat, dims = 1:20)

# Cluster the cells
scRNA_MC_seurat <- FindNeighbors(scRNA_MC_seurat, dims = 1:20)
scRNA_MC_seurat <- FindClusters(scRNA_MC_seurat, resolution = 0.5)

# TSNE clustering
scRNA_MC_seurat <- RunTSNE(scRNA_MC_seurat, dims = 1:10)
DimPlot(scRNA_MC_seurat, reduction = "tsne")
DimPlot(scRNA_MC_seurat, reduction = "tsne", label = TRUE)

# find markers for every cluster compared to all remaining cells
scRNA_MC <- FindAllMarkers(scRNA_MC_seurat, return.thresh = 1)
sc_MC <- scRNA_MC %>% group_by(gene) %>% filter(p_val == min(p_val)) %>% dplyr::slice(1)
sc_MC <- as.data.frame(sc_MC)
rownames(sc_MC) <- sc_MC$gene

save(spark_MC, scRNA_MC, sc_MC, file = "./Output/MC.RData")