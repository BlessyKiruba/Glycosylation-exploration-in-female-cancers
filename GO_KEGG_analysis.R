library(clusterProfiler)
library(org.Hs.eg.db) 
library(enrichplot)
library(ggplot2)
library(DOSE)
library(dplyr)

survival_genes <- c("ZNF90","SGPP2","SPON1","MUC16","CP","SLC34A2","LMO7", "ING3", "ST6GAL1","CCNB2","CDC20","UBE2C","DAPL1","C3","HIST3H2A","CLDN10","THSD4","SORL1","HOXC9","SLC16A3","AGR2", "EPS8L1","TNFAIP2","TMC5", "CDKN3","GALNT6","GPX3", "EPHA2","KIFC3","DEFB1","B4GALT1","MYH14","PPL","CEACAM1","TMPRSS4","A4GALT","WFDC2","APOL1","MUC20","SERPINA1","LTF",  "CFH",  "TES",  "TFF3", "PTTG1","PLCG2","PROM2","CXCL1","SLC40A1","IL1R1","VCAN", "DERL3","FBLN1","PIGR", "NFIX", "TGFA","IL1RN","PROM1","BMPR1B","ST6GALNAC1","FABP5","SLPI","SLC44A4")

gene_df <- bitr(
  survival_genes,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

ego_bp <- enrichGO(
  gene          = unique(gene_df$ENTREZID),
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",              
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1,
  readable      = TRUE
)
write.csv(as.data.frame(ego_bp),
          file = "/path/to/GO and KEGG analysis results/63survivalgenes_GO_BP.csv",
          row.names = FALSE)
tiff('/path/to/GO and KEGG analysis results/63Survival_genes_GOBP.tiff',width = 12, height = 10,units = 'in', res=600)
barplot(ego_bp, showCategory = 20, drop = TRUE) +
  ggtitle("63 Survival genes - GO Biological Process")
dev.off()


survival_genes_go <- enrichGO(
  gene   = survival_genes,
  OrgDb  = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable     = TRUE
)

write.csv(as.data.frame(survival_genes_go), "/path/to/GO and KEGG analysis results/63survivalgenesGO.csv", row.names = FALSE)

tiff('/path/to/GO and KEGG analysis results/63Survival_genes_GO.tiff',width = 12, height = 8,units = 'in', res=600)
barplot(survival_genes_go, showCategory = 15, title = "63 Survival genes - GO Enrichment")
dev.off()


ov_up_entrez <- bitr(ov_up_genes, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

ov_up_kegg <- enrichKEGG(
  gene         = ov_up_entrez$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.5,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.5
)

write.csv(as.data.frame(ov_up_kegg), "/path/to/GO and KEGG analysis results/63survivalgenesKEGG.csv", row.names = FALSE)

tiff('/path/to/GO and KEGG analysis results/63Survival_genes_KEGG.tiff',width = 12, height = 8,units = 'in', res=600)
barplot(ov_up_kegg, showCategory = 15, title = "63 Survival genes - KEGG")
dev.off()

commongenes_from_rf_and_gb <- c("MUC20","CP","SLC16A3","TFF3","DAPL1")

commongenes_from_rf_and_gb_go <- enrichGO(
  gene         = commongenes_from_rf_and_gb,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.5,
  qvalueCutoff = 0.5,
  readable     = TRUE
)

write.csv(as.data.frame(commongenes_from_rf_and_gb_go), "/path/to/GO and KEGG analysis results/commongenes_from_rf&gb_GO.csv", row.names = FALSE)

tiff('/path/to/GO and KEGG analysis results/commongenes_from_rf&gb_GO.tiff',width = 12, height = 8,units = 'in', res=600)
barplot(commongenes_from_rf_and_gb_go, showCategory = 15, title = "Common Genes from RF & GB - GO Enrichment")
dev.off()

commongenes_from_rf_and_gb_go_entrez <- bitr(commongenes_from_rf_and_gb, 
                       fromType = "SYMBOL", 
                       toType = "ENTREZID", 
                       OrgDb = org.Hs.eg.db)
print(paste("Number of genes to analyze:", length(ov_down_entrez$ENTREZID)))

commongenes_from_rf_and_gb_kegg <- enrichKEGG(
  gene         = commongenes_from_rf_and_gb_go_entrez $ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.5,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.5
)

write.csv(as.data.frame(commongenes_from_rf_and_gb_kegg), "/path/to/GO and KEGG analysis results/commongenes_from_rf&gb_KEGG.csv", row.names = FALSE)

tiff('/path/to/GO and KEGG analysis results/commongenes_from_rf&gb_KEGG.tiff',width = 12, height = 8,units = 'in', res=600)
barplot(commongenes_from_rf_and_gb_kegg, showCategory = 15, title = "Common Genes from RF & GB - KEGG")
dev.off()

glycoprotein_pathway_genes <- c("SLC35D1", "UST", "VCAN", "XYLT1", "XYLT2", "PARP2", "AGER", "HEXB", "HPSE", "HYAL4", "IDS", "IDUA", "LIPC", "SGSH", "EGFLAM", "GAL3ST3", "GCNT1", "GXYLT1", "GOLPH3", "SIRT3", "SIRT1", "NDNF", "SPOCK2", "SPOCK3", "ENTPD5", "ERP44", "UGGT1", "HS3ST4", "LECT1", "BMPR1B", "DCN", "GCNT3", "TMEM165", "TET1", "IGF1", "MUC16", "MAN2A1", "ART4", "ATP7A", "SLC35C2", "VANGL2", "PSEN1", "VEGFB", "PORCN", "SIRT6", "ALG1", "B3GALNT1", "C1GALT1", "PCSK6", "B3GALT1", "B3GALT2", "CHST8", "FUT3", "FUT8", "RPN2", "SIRT4", "SERP1", "MGAT2", "FUT7", "FUT1", "FUT10", "GCNT4", "FUT2", "GYLTL1B", "COL11A1", "FUT4", "FUT5", "FUT6", "FUT9", "ST3GAL1", "ST3GAL2", "COL2A1", "ST3GAL3", "ST3GAL4", "ST3GAL5", "ST3GAL6", "ST6GALNAC1", "TNIP1", "ST6GALNAC2", "ST6GALNAC3", "ST6GALNAC4", "ST6GALNAC5", "ST6GALNAC6", "GOLGA2", "ST8SIA2", "ARSB", "GORASP1", "ST8SIA3", "ST8SIA4", "NEU4", "ST8SIA6", "A4GALT", "A4GNT", "MUC13", "ABO", "MUC2", "MUC4", "MUC5AC", "ADAMTS5", "DPAGT1", "MUC6", "ADAMTSL1", "ADAMTSL4", "ALG10", "ALG10B", "ALG12", "SPON1", "ALG1L", "ALG1L2", "ALG2", "PARP3", "ALG3", "ALG6", "ALG8", "ALG9", "ART1", "IMPAD1", "ART3", "ART5", "B3GALNT2", "B3GALT4", "B3GALT5", "KCNE1", "B3GALTL", "B3GNT5", "B3GNT6", "B3GNT9", "C1GALT1C1", "C3ORF39", "C3ORF64", "CANT1", "CCDC126", "NUS1", "BMP2", "COG3", "COG7", "CSPG4", "SGK196", "DERL3", "DHDDS", "DOLPP1", "DPM3", "DPY19L1", "DPY19L2", "DPY19L2P2", "DPY19L3", "DPY19L4", "ASGR2", "EDEM1", "EDEM2", "EDEM3", "FKRP", "IL1B", "FKTN", "FOXL1", "FUT11", "GAL3ST4", "GALNT1", "GALNT10", "GALNT11", "GALNT12", "GALNT13", "GALNT14", "GALNT2", "GALNT3", "GALNT4", "VCP", "GALNT6", "GALNT7", "MUC20", "GALNT8", "GALNT9", "GALNTL1", "GALNTL2", "GALNTL4", "GALNTL5", "GALNTL6", "GBGT1", "DAD1", "GCNT6", "GCNT7", "GGTA1P", "GNPTAB", "GNPTG", "ADAMTS7", "THBS1", "B4GALNT1", "C20ORF173", "DAG1", "GXYLT2", "HS6ST3", "SRD5A3", "ISPD", "KIAA2018", "PPARD", "LMAN1", "LMF1", "LOC100288842", "MAGT1", "MAN1A1", "MAN1A2", "MAN1B1", "MAN1C1", "MAN2A2", "ALG5", "MCFD2", "MGAT3", "MGAT4A", "TIPARP", "MGAT4B", "MGAT4C", "MGAT5", "MGAT5B", "MOGS", "MUC1", "ADAMTS12", "MUC12", "ADAMTS13", "MUC15", "MUC17", "MUC19", "MUC21", "MUC3A", "MUC3B", "MUC5B", "MUC7", "MUCL1", "NAGPA", "ARF4", "DPM1", "NUDT14", "DPM2", "OAS2", "OSTC", "PGM3", "PARP1", "PARP10", "PARP16", "PARP4", "MGAT1", "PHLDA1", "B4GALNT2", "CHST4", "CHST7", "SULF1", "CSGALNACT1", "SULF2", "EXTL2", "LARGE", "LRP2", "MGEA5", "OGT", "ST6GAL1", "ACPL2", "B3GALT6", "B3GAT1", "PMM2", "B3GAT2", "POGLUT1", "B3GAT3", "POMGNT1", "B3GNT1", "POMT1", "B3GNT2", "POMT2", "B3GNT3", "B3GNT4", "B3GNT7", "GAL3ST1", "B3GNT8", "PRKCSH", "B4GALT2", "B4GALT3", "B4GALT4", "B4GALT5", "B4GALT6", "B4GALT7", "CFP", "BCAN", "BMPR2", "RPN1", "BGN", "CHPF", "CHPF2", "SDF2", "CHST11", "SDF2L1", "CHST12", "SERP2", "CHST13", "SIRT5", "CHST14", "ST6GAL2", "CHST15", "ST8SIA1", "CHST3", "ST8SIA5", "CHST9", "STT3A", "CHSY1", "STT3B", "CHSY3", "SYVN1", "CSGALNACT2", "TET2", "CSPG5", "TET3", "CYTL1", "DSE", "DSEL", "TMEM115", "EXT1", "TMEM5", "B4GALT1", "EXT2", "TNKS", "EXTL1", "TNKS2", "EXTL3", "TRAK1", "TRAK2", "GALNT5", "TUSC3", "GCNT2", "GLCE", "GPC1", "UBE2G2", "NPC1", "UBE2J1", "UGGT2", "HIF1A", "SLC34A1", "HEXA", "WBSCR17", "ADAMTS9", "HS3ST3B1", "HS3ST5", "DDOST", "HS6ST1", "SIRT2", "CELA1", "HS6ST2", "HYAL1", "FBXO2", "FBXO6", "NCAN", "NDST1", "NDST2", "POFUT1", "NDST3", "POFUT2", "NDST4", "NGLY1")

glycoprotein_pathway_genes_go <- enrichGO(
  gene         = glycoprotein_pathway_genes,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable     = TRUE
)

glycoprotein_pathway_genes_entrez <- bitr(glycoprotein_pathway_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
glycoprotein_pathway_genes_kegg <- enrichKEGG(
  gene         = glycoprotein_pathway_genes_entrez$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)
tiff('/path/to/GO and KEGG analysis results/Glycoprotein_pathway_genes_GO.tiff',width = 12, height = 8,units = 'in', res=600)
barplot(glycoprotein_pathway_genes_go, showCategory = 15, title = "Glycoprotein Pathway Genes - GO Enrichment")
dev.off()
tiff('/path/to/GO and KEGG analysis results/Glycoprotein_pathway_genes_KEGG.tiff',width = 12, height = 8,units = 'in', res=600)
barplot(glycoprotein_pathway_genes_kegg, showCategory = 15, title = "Glycoprotein Pathway Genes - KEGG Enrichment")
dev.off()
write.csv(as.data.frame(glycoprotein_pathway_genes_go), "/path/to/GO and KEGG analysis results/Glycoprotein_pathway_genes_GO.csv", row.names = FALSE)
write.csv(as.data.frame(glycoprotein_pathway_geneskegg), "/path/to/GO and KEGG analysis results/Glycoprotein_pathway_genes_KEGG.csv", row.names = FALSE)

