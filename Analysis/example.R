
rm(list = ls())
##------------------------------------------------------------------
## Load the spatial human breast cancer data
##------------------------------------------------------------------
# library(dplyr)

load("./RealData/Layer2_BC_Count.rds")

## extract the coordinates from the rawdata, i.e., location or coordinates
info_BC <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",1)),
                            y=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",2)),
                            total_counts=apply(rawcount,2,sum))
rownames(info_BC) <- colnames(rawcount)

##----------------------------------------------------------------
## Load the TCGA breast cancer data and make differential analysis
##----------------------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(DESeq2)

query.br <- GDCquery(project = "TCGA-BRCA",
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification", 
                     workflow.type = "HTSeq - Counts")

# All 1222 barcode
samples.br <- getResults(query.br,cols=c("cases"))
# Extract the 1102 TP (cancer) samples from the sampleDown
dataSmTP.br <- TCGAquery_SampleTypes(barcode = samples.br,
                                     typesample = "TP")
# Extract the 113 normal tissues samples
dataSmNT.br <- TCGAquery_SampleTypes(barcode = samples.br,
                                     typesample = "NT")
queryDown.br <- GDCquery(project = "TCGA-BRCA", 
                         data.category = "Transcriptome Profiling",
                         data.type = "Gene Expression Quantification", 
                         workflow.type = "HTSeq - Counts", 
                         barcode = c(dataSmTP.br, dataSmNT.br))

# Download data
GDCdownload(query = queryDown.br, method = "api")

dataPrep.br <- GDCprepare(query = queryDown.br, save = FALSE)

# Remove outliers
brca.count <- TCGAanalyze_Preprocessing(object = dataPrep.br,
                                        cor.cut = 0.6,
                                        datatype = "HTSeq - Counts")

# Set grouping information and build dds objects
group.br <- factor(c(rep("tumor", length(dataSmTP.br)), rep("normal", length(dataSmNT.br))), levels = c("tumor", "normal"))
colData.br <- data.frame(row.names = colnames(brca.count), group = group.br)
dds.br <- DESeqDataSetFromMatrix(countData = brca.count,
                                 colData = colData.br,
                                 design = ~ group)
# dispersion estimates and differential analysis
dds.br <- DESeq(dds.br)
res.br <- results(dds.br)
TCGA.br.result <- cbind(stat = res.br@listData$stat, pvalue = res.br@listData$pvalue)
# Obtain the gene names
genes.br <-mapIds(org.Hs.eg.db, 
                  keys=rownames(brca.count),
                  column="SYMBOL",  
                  keytype="ENSEMBL" ) 
# rownames(brca.count) <- genes
row.names(TCGA.br.result) <- genes.br
TCGA.br.result <- TCGA.br.result[complete.cases(TCGA.br.result),]

##-----------------------------------------------------------------
## Load the TCGA Thyroid cancer data and make differential analysis
##-----------------------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
# options(connectionObserver = NULL) # run if library(org.Hs.eg.db) failed
library(org.Hs.eg.db)
library(DESeq2)

query.th <- GDCquery(project = "TCGA-THCA",
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification", 
                     workflow.type = "HTSeq - Counts")

# All 568 barcode
samples.th <- getResults(query.th,cols=c("cases"))
# Extract the 502 TP (cancer) samples from the sampleDown
dataSmTP.th <- TCGAquery_SampleTypes(barcode = samples.th,
                                     typesample = "TP")
# Extract the 58 normal tissues samples
dataSmNT.th <- TCGAquery_SampleTypes(barcode = samples.th,
                                     typesample = "NT")
queryDown.th <- GDCquery(project = "TCGA-THCA", 
                         data.category = "Transcriptome Profiling",
                         data.type = "Gene Expression Quantification", 
                         workflow.type = "HTSeq - Counts", 
                         barcode = c(dataSmTP.th, dataSmNT.th))

# Download data
GDCdownload(query = queryDown.th, method = "api")

dataPrep.th <- GDCprepare(query = queryDown.th, save = FALSE)

# Remove outliers
thca.count <- TCGAanalyze_Preprocessing(object = dataPrep.th,
                                        cor.cut = 0.6,
                                        datatype = "HTSeq - Counts")

# Set grouping information and build dds objects
# group <- factor(c(rep("tumor",100),rep("normal",100)), levels = c("tumor","normal"))
group.th <- factor(c(rep("tumor", length(dataSmTP.th)), rep("normal", length(dataSmNT.th))), levels = c("tumor", "normal"))
colData.th <- data.frame(row.names=colnames(thca.count), group=group.th)
dds.th <- DESeqDataSetFromMatrix(countData = thca.count,
                                 colData = colData.th,
                                 design = ~ group)
# dispersion estimates and differential analysis
dds.th <- DESeq(dds.th)
res.th <- results(dds.th)
TCGA.th.result <- cbind(stat = res.th@listData$stat, padj = res.th$padj, pvalue = res.th@listData$pvalue, log2FoldChange = res.th$log2FoldChange)
# Obtain the gene names
genes.th <-mapIds(org.Hs.eg.db, 
                  keys=rownames(thca.count),
                  column="SYMBOL", 
                  keytype="ENSEMBL" ) 
# rownames(brca.count) <- genes
row.names(TCGA.th.result) <- genes.th
TCGA.th.result <- TCGA.th.result[complete.cases(TCGA.th.result),]

##-------------------------------------------------------------
## Gene-based association analysis based on GWAS data
##-------------------------------------------------------------
library(readr)
library(sumFREGAT)
library(stringr)

# Download GWAS summary data from http://bcac.ccge.medschl.cam.ac.uk/bcacdata/oncoarray/oncoarray-and-combined-summary-result/gwas-summary-results-breast-cancer-risk-2017/
oncoarray_bcac_public_release_oct17 <- read_table2("./Covariates/oncoarray_bcac_public_release_oct17.txt")
bcac.data <- oncoarray_bcac_public_release_oct17[,c("phase3_1kg_id",
                                                    "bcac_gwas_all_P1df",
                                                    "chr",
                                                    "position_b37",
                                                    "a1",
                                                    "bcac_gwas_all_beta",
                                                    "bcac_gwas_all_eaf_controls",
                                                    "a0")]
colnames(bcac.data) <- c("ID","P","CHROM","POS","EA","BETA","EAF","REF")
bcac.data <- as.data.frame(bcac.data)
bcac.data <- bcac.data[-which(bcac.data[,"P"]=="NULL"),]
bcac.data$ID <- str_split_fixed(bcac.data$ID,':',4)[,1]
bcac.data$P <- as.double(bcac.data$P)
prep.score.files(bcac.data[1:4], reference = "GWAS/ref1KG.MAC5.EUR_AF.RData",output.file.prefix = "GWAS/bcac")

GWAS.result <- ACAT(score.file = "GWAS/bcac.vcf.gz", gene.file = "hg38", genes = "all", anno.type = "",
                    beta.par = c(1, 1), weights.function = NULL, user.weights = FALSE,
                    gen.var.weights = "none", write.file = FALSE, quiet = FALSE)
row.names(GWAS.result) <- GWAS.result$gene
GWAS.result <- GWAS.result[-which(is.na(GWAS.result$pvalue)),]

##-------------------------------------------------------------------------------
## SE analysis with the proposed method incorporating multi-omics auxiliary data
##-------------------------------------------------------------------------------
library(SPARK)
source("./R/Cauchy.R")
source('./R/OrderShapeEM.R')
## filter genes and cells/spots (excluding the gene that are lowly expressed) and creat the spark object
spark_BC <- CreateSPARKObject(counts=rawcount, 
                              location=info_BC[,1:2],
                              percentage = 0.1, 
                              min_total_counts = 10)

## total counts for each cell/spot
spark_BC@lib_size <- apply(spark_BC@counts, 2, sum)

t1 <- proc.time()
## Estimating Parameter Under Null
spark_BC <- spark.vc(spark_BC, 
                     covariates = NULL, 
                     lib_size = spark_BC@lib_size, 
                     num_core = 5,
                     verbose = F)
spark_BC <- spark.test(spark_BC, 
                       check_positive = T, 
                       verbose = F)

## Combine the three auxiliary datasets using the Cauchy combination rule
aux.genes <- intersect(rownames(TCGA.br.result), rownames(TCGA.th.result))
aux.genes <- intersect(aux.genes, rownames(GWAS.result))
aux.Pvals <- rbind(TCGA.br.result[aux.genes,"pvalue"], 
                   GWAS.result[aux.genes,"pvalue"], 
                   TCGA.th.result[aux.genes,"pvalue"])
aux.pvalue <- Cauchy(aux.Pvals)
names(aux.pvalue) <- aux.genes

## Apply OrderShapeEM to conduct SE analysis incorporating 
target.genes <- intersect(aux.genes, rownames(spark_BC@counts))
pval_spark <- spark_BC@res_mtest[target.genes,"combined_pvalue"]
covar <- aux.pvalue[target.genes]

orderfdr.obj <- OrderShapeEM(pval_spark, covar, OrderShapeEM.control(trace = FALSE))
t1 <- proc.time()

sum(orderfdr.obj$fdr < 0.05)
rej.osem <- which(orderfdr.obj$fdr < 0.05)

