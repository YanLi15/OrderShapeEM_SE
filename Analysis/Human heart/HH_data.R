rm(list=ls())

##--------------------------------------------------------------
## scRNA-seq data analysis (Auxiliary covariate)
##--------------------------------------------------------------
require(Seurat)
require(data.table)
require(dplyr)
mat <- fread("./Covariates/exprMatrix.tsv.gz")
meta <- read.table("./Covariates/meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)
scRNA_HH_seurat <- CreateSeuratObject(counts = mat, project = "adultPancreas", meta.data = meta)

# QC and selecting cells
scRNA_HH_seurat[["percent.mt"]] <- PercentageFeatureSet(scRNA_HH_seurat, pattern = "^MT-")
# Normalizing the data
scRNA_HH_seurat <- NormalizeData(scRNA_HH_seurat)
# Identify highly expressed genes
scRNA_HH_seurat <- FindVariableFeatures(scRNA_HH_seurat, selection.method = "vst", nfeatures = 2000)
# Scaling the data
sc.genes <- rownames(scRNA_HH_seurat)
scRNA_HH_seurat <- ScaleData(scRNA_HH_seurat)

# PCA linear dimensional reduction
scRNA_HH_seurat <- RunPCA(scRNA_HH_seurat, features = VariableFeatures(object = scRNA_HH_seurat))

# Determine the dimensionality of the dataset
scRNA_HH_seurat <- JackStraw(scRNA_HH_seurat, num.replicate = 100)
scRNA_HH_seurat <- ScoreJackStraw(scRNA_HH_seurat, dims = 1:20)

# Cluster the cells
scRNA_HH_seurat <- FindNeighbors(scRNA_HH_seurat, dims = 1:20)
scRNA_HH_seurat <- FindClusters(scRNA_HH_seurat, resolution = 0.5)

# find markers for every cluster compared to all remaining cells
scRNA_HH <- FindAllMarkers(scRNA_HH_seurat, return.thresh = 1)
sc_HH <- scRNA_HH %>% group_by(gene) %>% filter(p_val == min(p_val)) %>% dplyr::slice(1)
sc_HH <- as.data.frame(sc_HH)
rownames(sc_HH) <- sc_HH$gene

##--------------------------------------------------------------
## scRNA-seq data analysis
##--------------------------------------------------------------
library(SPARK)
source('./R/spark.test.R')

sp_count <- Matrix::readMM("./RealData/V1_Human_Heart_filtered_feature_bc_matrix/matrix.mtx")
genes <- read.delim("./RealData/V1_Human_Heart_filtered_feature_bc_matrix/features.tsv", header=FALSE)
barcodes <- read.table("./RealData/V1_Human_Heart_filtered_feature_bc_matrix/barcodes.tsv", quote="\"", comment.char="")
location <- read.csv("./RealData/V1_Human_Heart_filtered_feature_bc_matrix/tissue_positions_list.csv", header=FALSE)

barcodes <- as.vector(t(barcodes))
rownames(sp_count) <- genes$V2
colnames(sp_count) <- barcodes
rownames(location) <- location$V1
location <- location[barcodes,5:6]
colnames(location) <- c("x","y")

# remove mt genes
mt_idx  <- grep("MT-",rownames(sp_count))
if(length(mt_idx)!=0){
  sp_count  <- sp_count[-mt_idx, ]
}
# remove genes that are not expressed in any spot
sp_count <- sp_count[-which(rowSums(as.matrix(sp_count)) == 0),] # 20904 * 4247

target.genes <- intersect(sc_HH$gene, rownames(sp_count))
counts <- sp_count[target.genes,]
## filter genes and cells/spots (excluding the gene that are lowly expressed) and creat the spark object
spark_HH <- CreateSPARKObject(counts     = counts, 
                              location   = location,
                              percentage = 0, 
                              min_total_counts = 0)
## total counts for each cell/spot
spark_HH@lib_size <- apply(spark_HH@counts, 2, sum)

## Estimating Parameter Under Null
spark_HH <- spark.vc(spark_HH, 
                     covariates = NULL, 
                     lib_size = spark_HH@lib_size, 
                     num_core = 10,
                     verbose = F)

spark_HH <- spark.test(spark_HH, 
                       check_positive = T, 
                       verbose = F)
# save(spark_HH, file = "./Output/HH.RData")

save(spark_HH, scRNA_HH, sc_HH, file = "./Output/HH.RData")