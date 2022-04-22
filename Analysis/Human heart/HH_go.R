##########################################################
#        Enrichment analysis of the HH results          #
##########################################################
rm(list = ls())
load("./Output/HH_result.RData")
load("./Output/HH.Rdata")
library(clusterProfiler)
library(org.Hs.eg.db)
# options(connectionObserver = NULL) # run if library(org.Hs.eg.db) failed

genes.all <- bitr(rownames(spark_HH@counts),fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
genes.spark <- bitr(genes.by,fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)

go.spark <- enrichGO(gene          = genes.spark$ENTREZID,
                     universe      = genes.all$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "All" ,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.99,
                     qvalueCutoff  = 0.99,
                     readable      = TRUE,
                     pool=TRUE)
sum(go.spark$p.adjust<0.05)  # 384

genes.osem <- bitr(genes.osem, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
go.osem <- enrichGO(gene          = genes.osem$ENTREZID,
                    universe      = genes.all$ENTREZID,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "All" ,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.99,
                    qvalueCutoff  = 0.99,
                    readable      = TRUE,
                    pool=TRUE)
sum(go.osem$p.adjust<0.05) # 499

go.overlap   <- intersect(go.spark@result[go.spark@result$p.adjust<0.05,]$ID,                              go.osem@result[go.osem@result$p.adjust<0.05,]$ID)
go.osem.only <- go.osem@result[go.osem@result$p.adjust<0.05,]
go.osem.only <- go.osem.only[!go.osem.only$ID%in%go.overlap,]
write.csv(go.osem.only,file = "./Output/HH_go_osem_only.csv",quote = FALSE)
