##########################################################
#              Validate the MOB results                  #
##########################################################
rm(list = ls())
load("./Output/MOB_result.RData")

length(intersect(genes.by, genes.osem))

genes.overlap    <- intersect(genes.by, genes.osem)  # 739 genes (the same for BH and SABHA)
genes.bh.only    <- genes.bh[!genes.bh%in%genes.overlap]
genes.sabha.only <- genes.sabha[!genes.sabha%in%genes.overlap]
genes.osem.only  <- genes.osem[!genes.osem%in%genes.overlap]

##--------------------------------------------------------------------
## The highlighted marker genes in the original study 
##--------------------------------------------------------------------
marker.genes <- c("Doc2g", "Slc17a7", "Reln", "Cdhr1", "Sv2b", "Shisa3", "Plcxd2", "Nmb", "Uchl1", "Rcan2")
length(intersect(genes.by, marker.genes)) # 7
length(intersect(genes.bh, marker.genes)) # 8
length(intersect(genes.sabha, marker.genes)) # 8
length(intersect(genes.osem, marker.genes)) # 8

##--------------------------------------------------------------------
## Single-cell-specific marker genes in the olfactory bulb (Tepe list) 
##--------------------------------------------------------------------
single.cell.genes <- read.table("./ValidateData/MOB/single-cell gene list.txt")
single.cell.genes <- t(single.cell.genes)

length(intersect(genes.by, single.cell.genes)) # 458
length(intersect(genes.bh, single.cell.genes)) # 595
length(intersect(genes.sabha, single.cell.genes)) # 659
length(intersect(genes.osem, single.cell.genes)) # 730

length(intersect(genes.overlap, single.cell.genes)) # 458
length(intersect(genes.bh.only, single.cell.genes)) # 137
length(intersect(genes.sabha.only, single.cell.genes)) # 201
length(intersect(genes.osem.only, single.cell.genes)) # 272

##-------------------------------------------------------------
## Genes related to the olfactory bulb in Harmonizome database
##-------------------------------------------------------------
library(rjson)
library(stringr)
# abs(Standardized Value) = -log10(p-value)
Harmonizome_Glomerular <- fromJSON(paste(readLines("./ValidateData/MOB/Harmonizome_Glomerular.json")))
Harmonizome_Granule    <- fromJSON(paste(readLines("./ValidateData/MOB/Harmonizome_Granule.json")))
Harmonizome_Mitral     <- fromJSON(paste(readLines("./ValidateData/MOB/Harmonizome_Mitral.json")))

Harmonizome_Glomerular <- Harmonizome_Glomerular[["associations"]]
Harmonizome_Granule    <- Harmonizome_Granule[["associations"]]
Harmonizome_Mitral     <- Harmonizome_Mitral[["associations"]]

Harmonizome.genes <- NULL
Harmonizome.pvalue <- NULL

for (i in 1:length(Harmonizome_Glomerular)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_Glomerular[[i]][["gene"]][["symbol"]])
  Harmonizome.pvalue <- c(Harmonizome.pvalue, 10^(-abs(Harmonizome_Glomerular[[i]][["standardizedValue"]])))
}

for (i in 1:length(Harmonizome_Granule)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_Granule[[i]][["gene"]][["symbol"]])
  Harmonizome.pvalue <- c(Harmonizome.pvalue, 10^(-abs(Harmonizome_Granule[[i]][["standardizedValue"]])))
}

for (i in 1:length(Harmonizome_Mitral)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_Mitral[[i]][["gene"]][["symbol"]])
  Harmonizome.pvalue <- c(Harmonizome.pvalue, 10^(-abs(Harmonizome_Mitral[[i]][["standardizedValue"]])))
}

names(Harmonizome.pvalue) <- Harmonizome.genes
Harmonizome.genes <- Harmonizome.genes[!duplicated(Harmonizome.genes)]
Harmonizome.genes <- tolower(Harmonizome.genes)
Harmonizome.genes <- str_to_title(Harmonizome.genes)

length(intersect(genes.by, Harmonizome.genes)) # 231
length(intersect(genes.bh, Harmonizome.genes)) # 339
length(intersect(genes.sabha, Harmonizome.genes)) # 392
length(intersect(genes.osem, Harmonizome.genes)) # 415

length(intersect(genes.overlap, Harmonizome.genes)) # 231
length(intersect(genes.bh.only, Harmonizome.genes)) # 108
length(intersect(genes.sabha.only, Harmonizome.genes)) # 161
length(intersect(genes.osem.only, Harmonizome.genes)) # 184

