##########################################################
#              Validate the MOB results                  #
##########################################################
rm(list = ls())
load("./Output/MC_result.RData")

length(intersect(genes.by, genes.osem))

genes.overlap    <- intersect(genes.by, genes.osem) # the same with BH and SABHA
genes.bh.only    <- genes.bh[!genes.bh%in%genes.overlap]
genes.sabha.only <- genes.sabha[!genes.sabha%in%genes.overlap]
genes.osem.only  <- genes.osem[!genes.osem%in%genes.overlap]

##---------------------------------------------------------------------------------
## spatially non-random genes in the Purkinje layer identified from Slide-seq data
##---------------------------------------------------------------------------------
SlideSeq.genes <- read.table("./ValidateData/MC/Slide-seq.txt")
SlideSeq.genes <- t(SlideSeq.genes) # 669 genes

length(intersect(genes.by, SlideSeq.genes)) # 104
length(intersect(genes.bh, SlideSeq.genes)) # 116
length(intersect(genes.sabha, SlideSeq.genes)) # 140
length(intersect(genes.osem, SlideSeq.genes)) # 138

##-------------------------------------------------------------
## Genes related to the cerebellum in Harmonizome database
##-------------------------------------------------------------
library(rjson)
library(stringr)
# abs(Standardized Value) = -log10(p-value)
Harmonizome_Cerebellum <- fromJSON(paste(readLines("./ValidateData/MC/Harmonizome_Cerebellum.json")))
Harmonizome_Cerebellum <- Harmonizome_Cerebellum[["associations"]]

Harmonizome.genes <- NULL
Harmonizome.pvalue <- NULL

for (i in 1:length(Harmonizome_Cerebellum)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_Cerebellum[[i]][["gene"]][["symbol"]])
  Harmonizome.pvalue <- c(Harmonizome.pvalue, 10^(-abs(Harmonizome_Cerebellum[[i]][["standardizedValue"]])))
}

names(Harmonizome.pvalue) <- Harmonizome.genes
Harmonizome.genes <- Harmonizome.genes[!duplicated(Harmonizome.genes)]
Harmonizome.genes <- tolower(Harmonizome.genes)
Harmonizome.genes <- str_to_title(Harmonizome.genes) # 1867 genes

length(intersect(genes.by, Harmonizome.genes)) # 66
length(intersect(genes.bh, Harmonizome.genes)) # 73
length(intersect(genes.sabha, Harmonizome.genes)) # 91
length(intersect(genes.osem, Harmonizome.genes)) # 85

length(intersect(genes.overlap, Harmonizome.genes))
length(intersect(genes.bh.only, Harmonizome.genes))
length(intersect(genes.sabha.only, Harmonizome.genes))
length(intersect(genes.osem.only, Harmonizome.genes))

