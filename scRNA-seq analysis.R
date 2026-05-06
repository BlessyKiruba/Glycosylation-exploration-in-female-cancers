############################
LOAD LIBRARIES
############################
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)

############################
 FUNCTIONS
############################

# Create Seurat object
create_seurat <- function(path, project_name, type_label) {
  data <- Read10X(path)
  
  obj <- CreateSeuratObject(
    counts = data,
    project = project_name,
    min.cells = 3,
    min.features = 200
  )
  
  obj$Type <- type_label
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  return(obj)
}

# QC filtering
qc_filter <- function(obj) {
  gene.counts <- obj$nFeature_RNA
  umi.counts  <- obj$nCount_RNA
  
  obj <- subset(
    obj,
    subset =
      nFeature_RNA > (mean(gene.counts) - 2 * sd(gene.counts)) &
      nFeature_RNA < (mean(gene.counts) + 2 * sd(gene.counts)) &
      nCount_RNA  > (mean(umi.counts)  - 2 * sd(umi.counts)) &
      nCount_RNA  < (mean(umi.counts)  + 2 * sd(umi.counts)) &
      percent.mt < 10
  )
  
  return(obj)
}

# Merge helper
merge_list <- function(seurat_list) {
  merge(x = seurat_list[[1]], y = seurat_list[-1])
}

# Add sample label
add_sample <- function(obj, label) {
  obj$sample <- label
  return(obj)
}

############################
 PATHS
############################

### Cervical Cancer
cc_paths <- list(
  CN1 = c("C:/.../GSE168652_RAW/Normal", "Cervical Control"),
  SCC1 = c("C:/.../GSE168652_RAW/Tumor", "Cervical Cancer"),
  
  ADC1 = c("C:/.../GSE197461_RAW/adc1", "Cervical Cancer"),
  ADC2 = c("C:/.../adc2", "Cervical Cancer"),
  ADC3 = c("C:/.../adc3", "Cervical Cancer"),
  ADC4 = c("C:/.../adc4", "Cervical Cancer"),
  ADC5 = c("C:/.../adc5", "Cervical Cancer"),
  
  SCC2 = c("C:/.../scc1", "Cervical Cancer"),
  SCC3 = c("C:/.../scc2", "Cervical Cancer"),
  SCC4 = c("C:/.../scc3", "Cervical Cancer"),
  
  SCC5 = c("C:/.../GSE208653_RAW/scc1", "Cervical Cancer"),
  SCC6 = c("C:/.../GSE208653_RAW/scc2", "Cervical Cancer"),
  
  ADC6 = c("C:/.../GSE208653_RAW/adc1", "Cervical Cancer"),
  
  CN2 = c("C:/.../GSE208653_RAW/n1", "Cervical Control"),
  CN3 = c("C:/.../GSE208653_RAW/n2", "Cervical Control")
)

### Endometrial Cancer
ec_paths <- list(
  EC1 = c("C:/.../c1", "Endometrial Cancer"),
  EC2 = c("C:/.../c2", "Endometrial Cancer"),
  EC3 = c("C:/.../c3", "Endometrial Cancer"),
  EC4 = c("C:/.../c4", "Endometrial Cancer"),
  EC5 = c("C:/.../c5", "Endometrial Cancer"),
  
  EN1 = c("C:/.../EN1", "Endometrial Control"),
  EN2 = c("C:/.../EN2", "Endometrial Control"),
  EN3 = c("C:/.../EN3", "Endometrial Control"),
  EN4 = c("C:/.../EN4", "Endometrial Control"),
  EN5 = c("C:/.../EN5", "Endometrial Control"),
  EN6 = c("C:/.../EN6", "Endometrial Control"),
  EN7 = c("C:/.../EN7", "Endometrial Control")
)

### Ovarian & Breast already folder-based
ov_base <- "C:/.../OV/GSE184880_RAW/"
bc_base <- "C:/.../BC/GSE243526_RAW/"

############################
 SEURAT LISTS
############################

# Cervical
seurat_list_CC <- lapply(names(cc_paths), function(name) {
  create_seurat(cc_paths[[name]][1], name, cc_paths[[name]][2])
})
names(seurat_list_CC) <- names(cc_paths)

# Endometrial
seurat_list_EC <- lapply(names(ec_paths), function(name) {
  create_seurat(ec_paths[[name]][1], name, ec_paths[[name]][2])
})
names(seurat_list_EC) <- names(ec_paths)

# Ovarian
ov_samples <- list.dirs(ov_base, full.names = FALSE, recursive = FALSE)
seurat_list_OV <- lapply(ov_samples, function(sample) {
  obj <- create_seurat(file.path(ov_base, sample), sample, NA)
  
  if (grepl("^OVC", sample)) obj$Type <- "Ovarian Cancer"
  if (grepl("^OVN", sample)) obj$Type <- "Ovarian Control"
  
  return(obj)
})
names(seurat_list_OV) <- ov_samples

# Breast
bc_samples <- list.dirs(bc_base, full.names = FALSE, recursive = FALSE)
seurat_list_BC <- lapply(bc_samples, function(sample) {
  obj <- create_seurat(file.path(bc_base, sample), sample, NA)
  
  if (grepl("^BC", sample)) obj$Type <- "Breast Cancer"
  if (grepl("^BN", sample)) obj$Type <- "Breast Control"
  
  return(obj)
})
names(seurat_list_BC) <- bc_samples

############################
 QC 
############################

seurat_list_CC <- lapply(seurat_list_CC, qc_filter)
seurat_list_EC <- lapply(seurat_list_EC, qc_filter)
seurat_list_OV <- lapply(seurat_list_OV, qc_filter)
seurat_list_BC <- lapply(seurat_list_BC, qc_filter)

############################
MERGING
############################

CC_merged <- merge_list(seurat_list_CC)
EC_merged <- merge_list(seurat_list_EC)
OV_merged <- merge_list(seurat_list_OV)
BC_merged <- merge_list(seurat_list_BC)

# Add labels
CC_merged <- add_sample(CC_merged, "CC")
EC_merged <- add_sample(EC_merged, "EC")
OV_merged <- add_sample(OV_merged, "OV")
BC_merged <- add_sample(BC_merged, "BC")

# Final merge
merged.data <- merge(
  x = BC_merged,
  y = list(OV_merged, CC_merged, EC_merged),
  add.cell.ids = c("BC","OV","CC","EC")
)

############################
NORMALIZATION + INTEGRATION
############################

merged.data <- NormalizeData(merged.data)
merged.data <- FindVariableFeatures(merged.data, nfeatures = 2000)
merged.data <- ScaleData(merged.data)
merged.data <- RunPCA(merged.data, npcs = 50)

# Before integration
merged.data <- RunUMAP(merged.data, dims = 1:30)
DimPlot(merged.data) + plot_annotation(title = "Before integration")

# Harmony
merged.data <- RunHarmony(merged.data, "orig.ident")

merged.data <- RunUMAP(merged.data, reduction = "harmony", dims = 1:30)
merged.data <- FindNeighbors(merged.data, reduction = "harmony", dims = 1:30)
merged.data <- FindClusters(merged.data, resolution = 1)

# After integration
DimPlot(merged.data, reduction = "umap", label = TRUE) + 
  plot_annotation(title = "After integration (Harmony)")

############################
FINAL STEP
############################

merged.data <- JoinLayers(merged.data)
