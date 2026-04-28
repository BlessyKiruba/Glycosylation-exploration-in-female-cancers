glycogenes_msig <- c("SLC35D1","UST",'VCAN','XYLT1','XYLT2','PARP2','AGER','HEXB','HPSE','HYAL4',
'IDS','IDUA','LIPC','SGSH','EGFLAM','GAL3ST3','GCNT1','GXYLT1','GOLPH3','SIRT3',
'SIRT1','NDNF','SPOCK2','SPOCK3','ENTPD5','ERP44','UGGT1','HS3ST4','LECT1','BMPR1B',
'DCN','GCNT3','TMEM165','TET1','IGF1','MUC16','MAN2A1','ART4','ATP7A','SLC35C2','VANGL2',
'PSEN1','VEGFB','PORCN','SIRT6','ALG1','B3GALNT1','C1GALT1','PCSK6','B3GALT1','B3GALT2',
'CHST8','FUT3','FUT8','RPN2','SIRT4','SERP1','MGAT2','FUT7','FUT1','FUT10','GCNT4','FUT2',
'GYLTL1B','COL11A1','FUT4','FUT5','FUT6','FUT9','ST3GAL1','ST3GAL2','COL2A1','ST3GAL3',
'ST3GAL4','ST3GAL5','ST3GAL6','ST6GALNAC1','TNIP1','ST6GALNAC2','ST6GALNAC3',
'ST6GALNAC4','ST6GALNAC5','ST6GALNAC6','GOLGA2','ST8SIA2','ARSB','GORASP1',
'ST8SIA3','ST8SIA4','NEU4','ST8SIA6','A4GALT','A4GNT','MUC13','ABO','MUC2','MUC4',
'MUC5AC','ADAMTS5','DPAGT1','MUC6','ADAMTSL1','ADAMTSL4','ALG10','ALG10B','ALG12',
'SPON1','ALG1L','ALG1L2','ALG2','PARP3','ALG3','ALG6','ALG8','ALG9','ART1','IMPAD1',
'ART3','ART5','B3GALNT2','B3GALT4','B3GALT5','KCNE1','B3GALTL','B3GNT5','B3GNT6',
'B3GNT9','C1GALT1C1','C3ORF39','C3ORF64','CANT1','CCDC126','NUS1','BMP2','COG3',
'COG7','CSPG4','SGK196','DERL3','DHDDS','DOLPP1','DPM3','DPY19L1','DPY19L2',
'DPY19L2P2','DPY19L3','DPY19L4','ASGR2','EDEM1','EDEM2','EDEM3','FKRP','IL1B',
'FKTN','FOXL1','FUT11','GAL3ST4','GALNT1','GALNT10','GALNT11','GALNT12','GALNT13',
'GALNT14','GALNT2','GALNT3','GALNT4','VCP','GALNT6','GALNT7','MUC20','GALNT8','GALNT9',
'GALNTL1','GALNTL2','GALNTL4','GALNTL5','GALNTL6','GBGT1','DAD1','GCNT6','GCNT7',
'GGTA1P','GNPTAB','GNPTG','ADAMTS7','THBS1','B4GALNT1','C20ORF173','DAG1','GXYLT2',
'HS6ST3','SRD5A3','ISPD','KIAA2018','PPARD','LMAN1','LMF1','LOC100288842',
'MAGT1','MAN1A1','MAN1A2','MAN1B1','MAN1C1','MAN2A2','ALG5','MCFD2','MGAT3',
'MGAT4A','TIPARP','MGAT4B','MGAT4C','MGAT5','MGAT5B','MOGS','MUC1','ADAMTS12',
'MUC12','ADAMTS13','MUC15','MUC17','MUC19','MUC21','MUC3A','MUC3B','MUC5B','MUC7',
'MUCL1','NAGPA','ARF4','DPM1','NUDT14','DPM2','OAS2','OSTC','PGM3','PARP1','PARP10',
'PARP16','PARP4','MGAT1','PHLDA1','B4GALNT2','CHST4','CHST7','SULF1','CSGALNACT1',
'SULF2','EXTL2','LARGE','LRP2','MGEA5','OGT','ST6GAL1','ACPL2','B3GALT6','B3GAT1','PMM2',
'B3GAT2','POGLUT1','B3GAT3','POMGNT1','B3GNT1','POMT1','B3GNT2','POMT2','B3GNT3','B3GNT4',
'B3GNT7','GAL3ST1','B3GNT8','PRKCSH','B4GALT2','B4GALT3','B4GALT4','B4GALT5','B4GALT6',
'B4GALT7','CFP','BCAN','BMPR2','RPN1','BGN','CHPF','CHPF2','SDF2','CHST11','SDF2L1','CHST12',
'SERP2','CHST13','SIRT5','CHST14','ST6GAL2','CHST15','ST8SIA1','CHST3','ST8SIA5','CHST9',
'STT3A','CHSY1','STT3B','CHSY3','SYVN1','CSGALNACT2','TET2','CSPG5','TET3','CYTL1','DSE',
'DSEL','TMEM115','EXT1','TMEM5','B4GALT1','EXT2','TNKS','EXTL1','TNKS2','EXTL3','TRAK1',
'TRAK2','GALNT5','TUSC3','GCNT2','GLCE','GPC1','UBE2G2','NPC1','UBE2J1','UGGT2','HIF1A',
'SLC34A1','HEXA','WBSCR17','ADAMTS9','HS3ST3B1','HS3ST5','DDOST',"HS6ST1",'SIRT2','CELA1',
'HS6ST2','HYAL1','FBXO2','FBXO6','NCAN','N',
'DST1','NDST2','POFUT1','NDST3','POFUT2','NDST4','NGLY1')




library(AUCell)
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
glycogenes <- glycogenes_msig

glycogenes <- intersect(glycogenes, rownames(epi))
length(glycogenes)

epi_epi <- subset(
  epi,
  subset = cellType == "Epithelial cells"
)

exprMat <- GetAssayData(
  epi_epi,
  assay = "RNA",
  slot = "data"   # log-normalized counts
)

exprMat <- as.matrix(exprMat)

cells_rankings <- AUCell_buildRankings(
  exprMat,
  nCores = 8,
  plotStats = TRUE
)


geneSets <- list(GlycoSignature = glycogenes)

auc_results <- AUCell_calcAUC(
  geneSets,
  cells_rankings,
  nCores = 8
)

epi_epi$Glyco_AUCell <- as.numeric(
  getAUC(auc_results)["GlycoSignature", ]
)


VlnPlot(
  epi_epi,
  features = "Glyco_AUCell",
  group.by = "CancerType",
  pt.size = 0
) +
  theme_classic()



library(ggplot2)

epi_epi@meta.data %>%
  ggplot(aes(x = CancerType, y = Glyco_AUCell, fill = CancerType)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  labs(
    y = "Glycogene AUCell score",
    x = ""
  ) +
  theme(legend.position = "none")



VlnPlot(
  epi_epi,
  features = "Glyco_AUCell",
  group.by = "CancerType",
  split.by = "Type",
  pt.size = 0.001,
  cols = pub.colors
)

#----------------
# Median
#----------------

median_auc <- median(epi_epi$Glyco_AUCell)

epi_epi$GT_AUC_Group <- ifelse(
  epi_epi$Glyco_AUCell >= median_auc,
  "High_GT_AUC",
  "Low_GT_AUC"
)

epi_epi$GT_AUC_Group <- factor(
  epi_epi$GT_AUC_Group,
  levels = c("Low_GT_AUC", "High_GT_AUC")
)

table(epi_epi$GT_AUC_Group)

Idents(epi_epi) <- epi_epi$GT_AUC_Group

deg_GT_AUC <- FindMarkers(
  epi_epi,
  ident.1 = "High_GT_AUC",
  ident.2 = "Low_GT_AUC",
  assay = "RNA",
  slot = "data",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"
)

deg_GT_AUC_sig <- deg_GT_AUC %>%
  filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) > 1
  )

nrow(deg_GT_AUC_sig)

GT_AUC_genes <- rownames(deg_GT_AUC_sig)

write.table(
  GT_AUC_genes,
  file = "GT_AUC_High_vs_Low_genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


VlnPlot(
  epi_epi,
  features = "Glyco_AUCell",
  group.by = "GT_AUC_Group",
  pt.size = 0
) +
  theme_classic()

top_genes <- head(
  rownames(deg_GT_AUC_sig[order(deg_GT_AUC_sig$avg_log2FC, decreasing = TRUE), ]),
  10
)

DoHeatmap(
  epi_epi,
  features = top_genes,
  group.by = "GT_AUC_Group"
)