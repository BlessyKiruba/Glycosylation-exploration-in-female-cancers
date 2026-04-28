############################################
#   CELLCHAT PIPELINE FOR MULTIPLE CANCERS
############################################

# Load packages
library(Seurat)
library(CellChat)
library(patchwork)
library(dplyr)
library(future)

############################################
# 1. INPUT 
############################################

# Path to your merged.data Seurat object
merged.data_file <- "PATH/TO/merged.data_seurat.rds"

# Column specifying cancer type
condition_column <- "Type"        # e.g. "Breast Cancer", "Cervical Cancer"

# Column specifying cell type annotation
celltype_column <- "cellType"

# Cancer groups to analyze
cancer_groups <- c("Breast Cancer", 
                   "Cervical Cancer",
                   "Endometrial Cancer",
                   "Ovarian Cancer")

# Minimum cells per cell type per group
min_cells <- 100

############################################
# 2. Load merged.data object and split by cancer
############################################

merged.data <- readRDS(merged.data_file)

# Check distribution
table(merged.data@meta.data[[condition_column]], merged.data@meta.data[[celltype_column]])

# Split by cancer type
seurat.list <- lapply(cancer_groups, function(x){
  subset(merged.data, subset = !!as.name(condition_column) == x)
})
names(seurat.list) <- cancer_groups

############################################
# 3. Remove cell types with too few cells
############################################

for(i in names(seurat.list)) {
  obj <- seurat.list[[i]]
  counts <- table(obj@meta.data[[celltype_column]])
  
  valid <- names(counts[counts > min_cells])
  obj <- subset(obj, subset = !!as.name(celltype_column) %in% valid)
  
  seurat.list[[i]] <- obj
}

############################################
# 4. Convert to CellChat objects
############################################

cellchat.list <- list()

for(i in names(seurat.list)) {
  message("Processing: ", i)
  
  data.input <- GetAssayData(seurat.list[[i]], assay = "RNA", slot = "data")
  meta <- seurat.list[[i]]@meta.data
  
  cellchat <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = celltype_column
  )
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  
  future::plan("multisession", workers = 8)
  options(future.globals.maxSize = 40 * 1024^3)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  cellchat.list[[i]] <- cellchat
}

############################################
# 5. Merge and compare all cancers
############################################

object.list <- cellchat.list

cellchat.merge <- mergeCellChat(
  object.list,
  add.names = names(cellchat.list)
)
cellchat.merge <- mergeCellChat(object.list, add.names = names(cellchat.list))

saveRDS(cellchat.merge, "cellchat_merge_4cancers.rds")

############################################
# 6. Visualizations (examples)
############################################

# Compare overall signaling strength
pdf("compare_net_strength.pdf", width = 10, height = 6)
compareInteractions(cellchat.merge, show.legend = TRUE)
dev.off()

# Compare number of interactions
pdf("compare_edge_count.pdf", width = 10, height = 6)
compareInteractions(cellchat.merge, show.legend = TRUE, group = 2)
dev.off()

# Heatmap of selected pathways
# Heatmap of selected pathways
pathways.show <- c("MIF", "TNF", "TGFb")

# 1. Compute centrality for merged object
cellchat.merge <- netAnalysis_computeCentrality(cellchat.merge)

# 2. Plot the selected pathways
pdf("compare_pathways.pdf", width = 10, height = 6)
netAnalysis_signalingRole_network(cellchat.merge, signaling = pathways.show)
dev.off()


############################################
# END OF SCRIPT
############################################


gg1 <- compareInteractions(cellchat.merge, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat.merge, show.legend = F, group = c(1,2), measure = "weight")
pdf("compareInteractionss.pdf", width = 10, height = 6)
gg1 + gg2

pdf("netVisual_heatmap.pdf", width = 10, height = 6)
We can also show differential number of interactions or interaction strength in a greater details using a heatmap. The top colored bar plot represents the sum of column of values displayed in the heatmap (incoming signaling). The right colored bar plot represents the sum of row of values (outgoing signaling). In the colorbar, $\color{red}{\text{red}}$ (or $\color{blue}{\text{blue}}$) represents $\color{red}{\text{increased}}$ (or $\color{blue}{\text{decreased}}$) signaling in the second dataset compared to the first one.
```{r, fig.width=10,fig.height = 5, fig.wide = TRUE, fig.align = "center"}
gg1 <- netVisual_heatmap(cellchat.merge)
gg2 <- netVisual_heatmap(cellchat.merge, measure = "weight")
gg1 + gg2
dev.off()

```
tiff("WEIGHT.tiff",width = 9000, height = 7000,res = 600)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

# High-resolution TIFF
tiff("circle_plots_all_cancers.tiff",
     width = 9000, height = 7000,
     res = 600)
# 2 rows x 2 columns
par(mfrow = c(2,2), mar = c(5,5,4,2), xpd = TRUE)  # Adjust margins
# Loop through all cancers
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count,
                   weight.scale = TRUE,
                   label.edge = FALSE,
                   edge.weight.max = weight.max[2],
                   edge.width.max = 12,
                   title.name = names(object.list)[i])  # Title is cancer name
}
dev.off()
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat.merge, slot.name = "netP") 
