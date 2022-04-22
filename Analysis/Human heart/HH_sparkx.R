rm(list = ls())
load("./Output/HH.Rdata")
source("./R/RunAll.R")
library(SPARK)

alpha = 0.05
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

sparkx_HH <- sparkx(counts, location, numCores=1, option="mixture")

# Primary pvalues
pval_sparkx <- sparkx_HH$res_mtest[target.genes,"combinedPval"]
# Auxiliary covariates
covar <- sc_HH[target.genes,"p_val"]

Rej.res <- RunAll(pval_sparkx, covar, alpha = alpha)

sparkx.by <- target.genes[Rej.res$rej.by]
sparkx.bh <- target.genes[Rej.res$rej.bh]
sparkx.sabha <- target.genes[Rej.res$rej.sabha]
sparkx.osem <- target.genes[Rej.res$rej.osem]

save(sparkx.by, sparkx.bh, sparkx.sabha, sparkx.osem, file = "./Output/HH_sparkx_result.RData")


#----------------------------------------------
# Validate
#----------------------------------------------

single.cell.genes <- t(read.table("./ValidateData/HH/Single-cell maker genes.txt"))

length(intersect(sparkx.by, single.cell.genes)) # 104
length(intersect(sparkx.bh, single.cell.genes)) # 148
length(intersect(sparkx.sabha, single.cell.genes)) # 160
length(intersect(sparkx.osem, single.cell.genes)) # 214


## Genes related to the olfactory bulb in Harmonizome database
library(rjson)
library(stringr)
# abs(Standardized Value) = -log10(p-value)
Harmonizome_Heart <- fromJSON(paste(readLines("./ValidateData/HH/Harmonizome_Heart.json")))
Harmonizome_Heart <- Harmonizome_Heart[["associations"]]

Harmonizome.genes <- NULL
Harmonizome.pvalue <- NULL

for (i in 1:length(Harmonizome_Heart)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_Heart[[i]][["gene"]][["symbol"]])
  Harmonizome.pvalue <- c(Harmonizome.pvalue, 10^(-abs(Harmonizome_Heart[[i]][["standardizedValue"]])))
}

names(Harmonizome.pvalue) <- Harmonizome.genes
Harmonizome.genes <- Harmonizome.genes[!duplicated(Harmonizome.genes)] # 6723 genes

length(intersect(sparkx.by, Harmonizome.genes)) # 194
length(intersect(sparkx.bh, Harmonizome.genes)) # 309
length(intersect(sparkx.sabha, Harmonizome.genes)) # 353
length(intersect(sparkx.osem, Harmonizome.genes)) # 465