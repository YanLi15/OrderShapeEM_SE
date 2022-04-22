###########################################################
#                Validate the BC results                  #
#               Take Case I as an example                 #
###########################################################
rm(list = ls())
load("./Output/BC_result.RData")

genes.overlap    <- intersect(genes.by, genes.bh) # the same with BH and SABHA
genes.bh.only    <- genes.bh7[!genes.bh%in%genes.overlap]
genes.sabha7.only <- genes.sabha7[!genes.sabha7%in%genes.overlap]
genes.osem1.only  <- genes.osem1[!genes.osem1%in%genes.overlap]
genes.osem2.only  <- genes.osem2[!genes.osem2%in%genes.overlap]
genes.osem3.only  <- genes.osem3[!genes.osem3%in%genes.overlap]
genes.osem4.only  <- genes.osem4[!genes.osem4%in%genes.overlap]
genes.osem5.only  <- genes.osem5[!genes.osem5%in%genes.overlap]
genes.osem6.only  <- genes.osem6[!genes.osem6%in%genes.overlap]
genes.osem7.only  <- genes.osem7[!genes.osem7%in%genes.overlap]
##----------------------------------------------------
## The highlighted marker genes in the original study 
##----------------------------------------------------
cancer.genes <- c("IGFBP5", "SPARC", "VIM", "FN1", "POSTN", "MUCL1", "PIP", 
                  "SCGB2A2", "GAS6", "KRT17", "PEG10", "AREG", "MMP14", "DCN")
length(intersect(genes.overlap, cancer.genes))
length(intersect(genes.bh.only, cancer.genes))
length(intersect(genes.sabha7.only, cancer.genes))
length(intersect(genes.osem1.only, cancer.genes))
length(intersect(genes.osem2.only, cancer.genes))
length(intersect(genes.osem3.only, cancer.genes))
length(intersect(genes.osem4.only, cancer.genes))
length(intersect(genes.osem5.only, cancer.genes))
length(intersect(genes.osem6.only, cancer.genes))
length(intersect(genes.osem7.only, cancer.genes))

##------------------------------------------------------------
## The breast cancer related genes in the CancerMine database 
##------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)

cancerMine <- read.table("./ValidateData/BC/cancermine_collated.csv", sep = ",", header = T)
cancerMine <- cancerMine[which(cancerMine["cancer_normalized"] == "breast cancer"), "gene_entrez_id"]
cancerMine <- as.character(cancerMine)
cancerMine.genes <- mapIds(org.Hs.eg.db,
                           keys=cancerMine,
                           column="SYMBOL",  # Target gene type
                           keytype="ENTREZID")
# 1642 genes

length(intersect(genes.overlap, cancerMine.genes))
length(intersect(genes.bh.only, cancerMine.genes))
length(intersect(genes.sabha7.only, cancerMine.genes))
length(intersect(genes.osem1.only, cancerMine.genes))
length(intersect(genes.osem2.only, cancerMine.genes))
length(intersect(genes.osem3.only, cancerMine.genes))
length(intersect(genes.osem4.only, cancerMine.genes))
length(intersect(genes.osem5.only, cancerMine.genes))
length(intersect(genes.osem6.only, cancerMine.genes))
length(intersect(genes.osem7.only, cancerMine.genes))

##-------------------------------------------------------------
## Genes related to breast cancer in Harmonizome database
##-------------------------------------------------------------
library(rjson)
# abs(Standardized Value) = -log10(p-value)
Harmonizome_DISEASES <- fromJSON(paste(readLines("./ValidateData/BC/Harmonizome_DISEASES.json")))
Harmonizome_GAD      <- fromJSON(paste(readLines("./ValidateData/BC/Harmonizome_GAD.json")))
Harmonizome_GWAS     <- fromJSON(paste(readLines("./ValidateData/BC/Harmonizome_GWAS.json")))
Harmonizome_OMIM     <- fromJSON(paste(readLines("./ValidateData/BC/Harmonizome_OMIM.json")))
Harmonizome_PSP      <- fromJSON(paste(readLines("./ValidateData/BC/Harmonizome_PSP.json")))
Harmonizome_GWASdb   <- fromJSON(paste(readLines("./ValidateData/BC/Harmonizome_GWASdb.json")))

Harmonizome_DISEASES <- Harmonizome_DISEASES[["associations"]]
Harmonizome_GAD      <- Harmonizome_GAD[["associations"]]
Harmonizome_GWAS     <- Harmonizome_GWAS[["associations"]]
Harmonizome_OMIM     <- Harmonizome_OMIM[["associations"]]
Harmonizome_PSP      <- Harmonizome_PSP[["associations"]]
Harmonizome_GWASdb     <- Harmonizome_GWASdb[["associations"]]

Harmonizome.genes <- NULL
for (i in 1:length(Harmonizome_DISEASES)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_DISEASES[[i]][["gene"]][["symbol"]])
}

for (i in 1:length(Harmonizome_GAD)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_GAD[[i]][["gene"]][["symbol"]])
}

for (i in 1:length(Harmonizome_GWAS)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_GWAS[[i]][["gene"]][["symbol"]])
}

for (i in 1:length(Harmonizome_OMIM)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_OMIM[[i]][["gene"]][["symbol"]])
}

for (i in 1:length(Harmonizome_PSP)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_PSP[[i]][["gene"]][["symbol"]])
}

for (i in 1:length(Harmonizome_GWASdb)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_GWASdb[[i]][["gene"]][["symbol"]])
}

Harmonizome.genes <- Harmonizome.genes[!duplicated(Harmonizome.genes)]

# 3505 genes

length(intersect(genes.overlap, Harmonizome.genes))
length(intersect(genes.bh.only, Harmonizome.genes))
length(intersect(genes.sabha7.only, Harmonizome.genes))
length(intersect(genes.osem1.only, Harmonizome.genes))
length(intersect(genes.osem2.only, Harmonizome.genes))
length(intersect(genes.osem3.only, Harmonizome.genes))
length(intersect(genes.osem4.only, Harmonizome.genes))
length(intersect(genes.osem5.only, Harmonizome.genes))
length(intersect(genes.osem6.only, Harmonizome.genes))
length(intersect(genes.osem7.only, Harmonizome.genes))
