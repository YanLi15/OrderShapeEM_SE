rm(list = ls())
load("./Output/HH_result.RData")

##---------------------------------------------------------------------
## Single-cell-specific marker genes in the human heart (Cui)
##---------------------------------------------------------------------
single.cell.genes <- t(read.table("./ValidateData/HH/Single-cell maker genes.txt"))

length(intersect(genes.by, single.cell.genes)) # 121
length(intersect(genes.bh, single.cell.genes)) # 139
length(intersect(genes.sabha, single.cell.genes)) # 140
length(intersect(genes.osem, single.cell.genes)) # 199
##-------------------------------------------------------------
## Genes related to heart cell types in Harmonizome database
##-------------------------------------------------------------
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
# Harmonizome.genes <- tolower(Harmonizome.genes)
# Harmonizome.genes <- str_to_title(Harmonizome.genes) # 6723 genes

length(intersect(genes.by, Harmonizome.genes)) # 122
length(intersect(genes.bh, Harmonizome.genes)) # 153
length(intersect(genes.sabha, Harmonizome.genes)) # 157
length(intersect(genes.osem, Harmonizome.genes)) # 261