##########################################################
#        nEnrichment analysis of the BC results          #
##########################################################
rm(list=ls())
load("./Output/BC_result.RData")
load("./Output/BC.Rdata")
library(clusterProfiler)
library(org.Hs.eg.db)

genes.all <- bitr(rownames(spark_BC@counts),fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
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
sum(go.spark$p.adjust<0.05) # 553

genes.osem2 <- bitr(genes.osem2, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
go.osem2 <- enrichGO(gene          = genes.osem2$ENTREZID,
                     universe      = genes.all$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "All" ,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.99,
                     qvalueCutoff  = 0.99,
                     readable      = TRUE,
                     pool=TRUE)
sum(go.osem2$p.adjust<0.05) # 650

go.overlap   <- intersect(go.spark@result[go.spark@result$p.adjust<0.05,]$ID, 
                          go.osem2@result[go.osem2@result$p.adjust<0.05,]$ID)
go.osem.only <- go.osem2@result[go.osem2@result$p.adjust<0.05,]
go.osem.only <- go.osem.only[!go.osem.only$ID%in%go.overlap,]
write.csv(go.osem.only,file = "./Output/BC_go_osem_only.csv",quote = FALSE)
