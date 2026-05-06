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
# 6. Visualizations 
############################################

# Compare overall signaling strength
compareInteractions(cellchat.merge, show.legend = TRUE)

# Compare number of interactions
compareInteractions(cellchat.merge, show.legend = TRUE, group = 2)
