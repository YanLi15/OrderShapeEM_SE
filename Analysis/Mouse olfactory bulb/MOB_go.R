##########################################################
#        Enrichment analysis of the MOB results          #
##########################################################
rm(list = ls())
load("./Output/MOB_result.RData")
load("./Output/MOB.Rdata")
library(clusterProfiler)
library(org.Mm.eg.db)
# options(connectionObserver = NULL) # run if library(org.Mm.eg.db) failed

genes.all <- bitr(rownames(spark_MOB@counts),fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
genes.spark <- bitr(genes.by,fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)

go.spark <- enrichGO(gene          = genes.spark$ENTREZID,
                     universe      = genes.all$ENTREZID,
                     OrgDb         = org.Mm.eg.db,
                     ont           = "All" ,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.99,
                     qvalueCutoff  = 0.99,
                     readable      = TRUE,
                     pool=TRUE)
sum(go.spark$p.adjust<0.05)  # 1281

genes.osem <- bitr(genes.osem, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
go.osem <- enrichGO(gene          = genes.osem$ENTREZID,
                    universe      = genes.all$ENTREZID,
                    OrgDb         = org.Mm.eg.db,
                    ont           = "All" ,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.99,
                    qvalueCutoff  = 0.99,
                    readable      = TRUE,
                    pool=TRUE)
sum(go.osem$p.adjust<0.05) # 1557

go.overlap   <- intersect(go.spark@result[go.spark@result$p.adjust<0.05,]$ID,                              go.osem@result[go.osem@result$p.adjust<0.05,]$ID)
go.osem.only <- go.osem@result[go.osem@result$p.adjust<0.05,]
go.osem.only <- go.osem.only[!go.osem.only$ID%in%go.overlap,]
write.csv(go.osem.only,file = "./Output/MOB_go_osem_only.csv",quote = FALSE)
