library(dplyr)
library(Seurat)
library(patchwork)
library(presto)
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(harmony)
library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)

CN1 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\CC\\GSE168652_RAW\\Normal")
CN1.data <- CreateSeuratObject(CN1,project = "CN1", min.cells = 3, min.features = 200)
SCC1 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\CC\\GSE168652_RAW\\Tumor")
SCC1.data <- CreateSeuratObject(SCC1,project = "CCSCC1", min.cells = 3, min.features = 200)

ADC1 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\CC\\GSE197461_RAW\\adc1")
ADC1.data <- CreateSeuratObject(ADC1,project = "CCADC1", min.cells = 3, min.features = 200)
ADC2 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\CC\\GSE197461_RAW\\adc2")
ADC2.data <- CreateSeuratObject(ADC2,project = "CCADC2", min.cells = 3, min.features = 200)

ADC3 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\CC\\GSE197461_RAW\\adc3")
ADC3.data <- CreateSeuratObject(ADC3,project = "CCADC3", min.cells = 3, min.features = 200)

ADC4 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\CC\\GSE197461_RAW\\adc4")
ADC4.data <- CreateSeuratObject(ADC4,project = "CCADC4", min.cells = 3, min.features = 200)
ADC5 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\CC\\GSE197461_RAW\\adc5")
ADC5.data <- CreateSeuratObject(ADC5,project = "CCADC5", min.cells = 3, min.features = 200)


SCC2 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\CC\\GSE197461_RAW\\scc1")
SCC2.data <- CreateSeuratObject(SCC2,project = "CCSCC2", min.cells = 3, min.features = 200)

SCC3 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\CC\\GSE197461_RAW\\scc2")
SCC3.data <- CreateSeuratObject(SCC3,project = "CCSCC3", min.cells = 3, min.features = 200)

SCC4 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\CC\\GSE197461_RAW\\scc3")
SCC4.data <- CreateSeuratObject(SCC4,project = "CCSCC4", min.cells = 3, min.features = 200)

SCC5 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\CC\\GSE208653_RAW\\scc1")
SCC5.data <- CreateSeuratObject(SCC5,project = "CCSCC5", min.cells = 3, min.features = 200)

SCC6 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\CC\\GSE208653_RAW\\scc2")
SCC6.data <- CreateSeuratObject(SCC6,project = "CCSCC6", min.cells = 3, min.features = 200)

ADC6 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\CC\\GSE208653_RAW\\adc1")
ADC6.data <- CreateSeuratObject(ADC6,project = "CCADC6", min.cells = 3, min.features = 200)

CN2 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\CC\\GSE208653_RAW\\n1")
CN2.data <- CreateSeuratObject(CN2,project = "CN2", min.cells = 3, min.features = 200)

CN3 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\CC\\GSE208653_RAW\\n2")
CN3.data <- CreateSeuratObject(CN3,project = "CN3", min.cells = 3, min.features = 200)



CN1.data$Type <- "Cervical Control"
CN2.data$Type <- "Cervical Control"
CN3.data$Type <- "Cervical Control"
SCC1.data$Type <- "Cervical Cancer"
SCC2.data$Type <- "Cervical Cancer"
SCC3.data$Type <- "Cervical Cancer"
SCC4.data$Type <- "Cervical Cancer"
SCC5.data$Type <- "Cervical Cancer"
SCC6.data$Type <- "Cervical Cancer"
ADC1.data$Type <- "Cervical Cancer"
ADC2.data$Type <- "Cervical Cancer"
ADC3.data$Type <- "Cervical Cancer"
ADC4.data$Type <- "Cervical Cancer"
ADC5.data$Type <- "Cervical Cancer"
ADC6.data$Type <- "Cervical Cancer"

CN1.data[["percent.mt"]]  <- PercentageFeatureSet(CN1.data, pattern = "^MT-")
CN2.data[["percent.mt"]]  <- PercentageFeatureSet(CN2.data, pattern = "^MT-")
CN3.data[["percent.mt"]]  <- PercentageFeatureSet(CN3.data, pattern = "^MT-")
SCC1.data[["percent.mt"]]  <- PercentageFeatureSet(SCC1.data, pattern = "^MT-")
SCC2.data[["percent.mt"]]  <- PercentageFeatureSet(SCC2.data, pattern = "^MT-")
SCC3.data[["percent.mt"]]  <- PercentageFeatureSet(SCC3.data, pattern = "^MT-")
SCC4.data[["percent.mt"]]  <- PercentageFeatureSet(SCC4.data, pattern = "^MT-")
SCC5.data[["percent.mt"]]  <- PercentageFeatureSet(SCC5.data, pattern = "^MT-")
SCC6.data[["percent.mt"]]  <- PercentageFeatureSet(SCC6.data, pattern = "^MT-")
ADC1.data[["percent.mt"]]  <- PercentageFeatureSet(ADC1.data, pattern = "^MT-")
ADC2.data[["percent.mt"]]  <- PercentageFeatureSet(ADC2.data, pattern = "^MT-")
ADC3.data[["percent.mt"]]  <- PercentageFeatureSet(ADC3.data, pattern = "^MT-")
ADC4.data[["percent.mt"]]  <- PercentageFeatureSet(ADC4.data, pattern = "^MT-")
ADC5.data[["percent.mt"]]  <- PercentageFeatureSet(ADC5.data, pattern = "^MT-")
ADC6.data[["percent.mt"]]  <- PercentageFeatureSet(ADC6.data, pattern = "^MT-")


VlnPlot(CN1.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
VlnPlot(CN2.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
VlnPlot(CN3.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
VlnPlot(SCC1.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
VlnPlot(SCC2.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
VlnPlot(SCC3.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
VlnPlot(SCC4.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
VlnPlot(SCC5.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
VlnPlot(SCC6.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
VlnPlot(ADC1.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
VlnPlot(ADC2.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
VlnPlot(ADC3.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
VlnPlot(ADC4.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
VlnPlot(ADC5.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
VlnPlot(ADC6.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)


gene.counts <- CN2.data$nFeature_RNA
umi.counts <- CN2.data$nCount_RNA

upper.limit.genes <- mean(gene.counts) + 2 * sd(gene.counts)
lower.limit.genes <- mean(gene.counts) - 2 * sd(gene.counts)
upper.limit.umis  <- mean(umi.counts)  + 2 * sd(umi.counts)
lower.limit.umis  <- mean(umi.counts)  - 2 * sd(umi.counts)

CN2.data <- subset(
  CN2.data,
  subset = nFeature_RNA > lower.limit.genes & 
    nFeature_RNA < upper.limit.genes &
    nCount_RNA > lower.limit.umis  & 
    nCount_RNA < upper.limit.umis &
    percent.mt < 10
)


seurat_list_CC <- list(
CN1.data,CN2.data, CN3.data,
  SCC1.data, SCC2.data, SCC3.data, SCC4.data,SCC5.data, SCC6.data,
  ADC1.data, ADC2.data, ADC3.data, ADC4.data, ADC5.data, ADC6.data
)


##########################
#####ENDOMETRIAL CANCER##
#########################

EC1 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\EC\\GSE173682_RAW\\c1")
EC1.data <- CreateSeuratObject(EC1,project = "EC1", min.cells = 3, min.features = 200)

EC2 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\EC\\GSE173682_RAW\\c2")
EC2.data <- CreateSeuratObject(EC2,project = "EC2", min.cells = 3, min.features = 200)

EC3 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\EC\\GSE173682_RAW\\c3")
EC3.data <- CreateSeuratObject(EC3,project = "EC3", min.cells = 3, min.features = 200)

EC4 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\EC\\GSE173682_RAW\\c4")
EC4.data <- CreateSeuratObject(EC4,project = "EC4", min.cells = 3, min.features = 200)

EC5 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\EC\\GSE173682_RAW\\c5")
EC5.data <- CreateSeuratObject(EC5,project = "EC5", min.cells = 3, min.features = 200)



EN1 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\EC\\GSE214411\\EN1")
EN1.data <- CreateSeuratObject(EN1,project = "EN1", min.cells = 3, min.features = 200)

EN2 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\EC\\GSE214411\\EN2")
EN2.data <- CreateSeuratObject(EN2,project = "EN2", min.cells = 3, min.features = 200)

EN3 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\EC\\GSE214411\\EN3")
EN3.data <- CreateSeuratObject(EN3,project = "EN3", min.cells = 3, min.features = 200)

EN4 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\EC\\GSE214411\\EN4")
EN4.data <- CreateSeuratObject(EN4,project = "EN4", min.cells = 3, min.features = 200)

EN5 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\EC\\GSE214411\\EN5")
EN5.data <- CreateSeuratObject(EN5,project = "EN5", min.cells = 3, min.features = 200)

EN6 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\EC\\GSE214411\\EN6")
EN6.data <- CreateSeuratObject(EN6,project = "EN6", min.cells = 3, min.features = 200)

EN7 <- Read10X("C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\EC\\GSE214411\\EN7")
EN7.data <- CreateSeuratObject(EN7,project = "EN7", min.cells = 3, min.features = 200)



EC1.data$Type <- "Endometrial Cancer"
EC2.data$Type <- "Endometrial Cancer"
EC3.data$Type <- "Endometrial Cancer"
EC4.data$Type <- "Endometrial Cancer"
EC5.data$Type <- "Endometrial Cancer"

EN1.data$Type <- "Endometrial Control"
EN2.data$Type <- "Endometrial Control"
EN3.data$Type <- "Endometrial Control"
EN4.data$Type <- "Endometrial Control"
EN5.data$Type <- "Endometrial Control"
EN6.data$Type <- "Endometrial Control"
EN7.data$Type <- "Endometrial Control"

EC1.data[["percent.mt"]]  <- PercentageFeatureSet(EC1.data, pattern = "^MT-")
EC2.data[["percent.mt"]]  <- PercentageFeatureSet(EC2.data, pattern = "^MT-")
EC3.data[["percent.mt"]]  <- PercentageFeatureSet(EC3.data, pattern = "^MT-")
EC4.data[["percent.mt"]]  <- PercentageFeatureSet(EC4.data, pattern = "^MT-")
EC5.data[["percent.mt"]]  <- PercentageFeatureSet(EC5.data, pattern = "^MT-")

EN1.data[["percent.mt"]]  <- PercentageFeatureSet(EN1.data, pattern = "^MT-")
EN2.data[["percent.mt"]]  <- PercentageFeatureSet(EN2.data, pattern = "^MT-")
EN3.data[["percent.mt"]]  <- PercentageFeatureSet(EN3.data, pattern = "^MT-")
EN4.data[["percent.mt"]]  <- PercentageFeatureSet(EN4.data, pattern = "^MT-")
EN5.data[["percent.mt"]]  <- PercentageFeatureSet(EN5.data, pattern = "^MT-")
EN6.data[["percent.mt"]]  <- PercentageFeatureSet(EN6.data, pattern = "^MT-")
EN7.data[["percent.mt"]]  <- PercentageFeatureSet(EN7.data, pattern = "^MT-")

VlnPlot(EC1.data, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(EC2.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 3)
VlnPlot(EC3.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 3)
VlnPlot(EC4.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 3)
VlnPlot(EC5.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 3)

VlnPlot(EN1.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 3)
VlnPlot(EN2.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 3)
VlnPlot(EN3.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 3)
VlnPlot(EN4.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 3)
VlnPlot(EN5.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 3)
VlnPlot(EN6.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 3)
VlnPlot(EN7.data, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 3)


seurat_list_EC <- list(EC1.data, EC2.data, EC3.data, EC4.data, EC5.data,
   EN1.data, EN2.data, EN3.data,EN4.data,EN5.data, EN6.data, EN7.data
)

##########################
#####OVARIAN CANCER######
#########################

library(Seurat)
library(dplyr)

# Define the base directory containing all sample folders
base_dir <- "C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\OV\\GSE184880_RAW\\"

# List all folder names (assuming each folder is a sample like c1, c2, ..., n6)
sample_folders <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)

# Create a named list to store Seurat objects
seurat_list <- list()

# Loop through each folder and create Seurat object
for (sample in sample_folders) {
  sample_path <- file.path(base_dir, sample)
  data <- Read10X(sample_path)
  seurat_obj <- CreateSeuratObject(data, project = sample, min.cells = 3, min.features = 200)
  seurat_list[[sample]] <- seurat_obj
}
seurat_list[["OVC1"]]$Type <- "Ovarian Cancer"
seurat_list[["OVC2"]]$Type <- "Ovarian Cancer"
seurat_list[["OVC3"]]$Type <- "Ovarian Cancer"
seurat_list[["OVC4"]]$Type <- "Ovarian Cancer"
seurat_list[["OVC5"]]$Type <- "Ovarian Cancer"
seurat_list[["OVC6"]]$Type <- "Ovarian Cancer"
seurat_list[["OVC7"]]$Type <- "Ovarian Cancer"
seurat_list[["OVC8"]]$Type <- "Ovarian Cancer"
seurat_list[["OVC9"]]$Type <- "Ovarian Cancer"
seurat_list[["OVC10"]]$Type <- "Ovarian Cancer"
seurat_list[["OVC11"]]$Type <- "Ovarian Cancer"
seurat_list[["OVC12"]]$Type <- "Ovarian Cancer"

seurat_list[["OVN1"]]$Type <- "Ovarian Control"
seurat_list[["OVN2"]]$Type <- "Ovarian Control"
seurat_list[["OVN3"]]$Type <- "Ovarian Control"
seurat_list[["OVN4"]]$Type <- "Ovarian Control"
seurat_list[["OVN5"]]$Type <- "Ovarian Control"

seurat_list[["OVC1"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list[["OVC1"]], pattern = "^MT-")
seurat_list[["OVC2"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list[["OVC2"]], pattern = "^MT-")
seurat_list[["OVC3"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list[["OVC3"]], pattern = "^MT-")
seurat_list[["OVC4"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list[["OVC4"]], pattern = "^MT-")
seurat_list[["OVC5"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list[["OVC5"]], pattern = "^MT-")
seurat_list[["OVC6"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list[["OVC6"]], pattern = "^MT-")
seurat_list[["OVC7"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list[["OVC7"]], pattern = "^MT-")
seurat_list[["OVC8"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list[["OVC8"]], pattern = "^MT-")
seurat_list[["OVC9"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list[["OVC9"]], pattern = "^MT-")
seurat_list[["OVC10"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list[["OVC10"]], pattern = "^MT-")
seurat_list[["OVC11"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list[["OVC11"]], pattern = "^MT-")
seurat_list[["OVC12"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list[["OVC12"]], pattern = "^MT-")

seurat_list[["OVN1"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list[["OVN1"]], pattern = "^MT-")
seurat_list[["OVN2"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list[["OVN2"]], pattern = "^MT-")
seurat_list[["OVN3"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list[["OVN3"]], pattern = "^MT-")
seurat_list[["OVN4"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list[["OVN4"]], pattern = "^MT-")
seurat_list[["OVN5"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list[["OVN5"]], pattern = "^MT-")

VlnPlot(seurat_list[["OVC1"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list[["OVC2"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list[["OVC3"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list[["OVC4"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list[["OVC5"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list[["OVC6"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list[["OVC7"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list[["OVC8"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list[["OVC9"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list[["OVC10"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list[["OVC11"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list[["OVC12"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)

VlnPlot(seurat_list[["OVN1"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list[["OVN2"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list[["OVN3"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list[["OVN4"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list[["OVN5"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)

seurat_list_OV <- seurat_list
##########################
##### BREAST CANCER#####
#########################

library(Seurat)
library(dplyr)

# Define the base directory containing all sample folders
base_dir <- "C:\\Users\\Admin\\Documents\\Gyn cancer_scRNA\\female_cancer_scRNA_seq\\BC\\GSE243526_RAW\\"

# List all folder names (assuming each folder is a sample like c1, c2, ..., n6)
sample_folders <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)

# Create a named list to store Seurat objects
seurat_list_BC <- list()

# Loop through each folder and create Seurat object
for (sample in sample_folders) {
  sample_path <- file.path(base_dir, sample)
  data <- Read10X(sample_path)
  seurat_obj_BC <- CreateSeuratObject(data, project = sample, min.cells = 3, min.features = 200)
  seurat_list_BC[[sample]] <- seurat_obj_BC
}

seurat_list_BC[["BC1"]]$Type <- "Breast Cancer"
seurat_list_BC[["BC2"]]$Type <- "Breast Cancer"
seurat_list_BC[["BC3"]]$Type <- "Breast Cancer"
seurat_list_BC[["BC4"]]$Type <- "Breast Cancer"
seurat_list_BC[["BC5"]]$Type <- "Breast Cancer"
seurat_list_BC[["BC6"]]$Type <- "Breast Cancer"
seurat_list_BC[["BC7"]]$Type <- "Breast Cancer"
seurat_list_BC[["BC8"]]$Type <- "Breast Cancer"
seurat_list_BC[["BC9"]]$Type <- "Breast Cancer"
seurat_list_BC[["BC10"]]$Type <- "Breast Cancer"
seurat_list_BC[["BC11"]]$Type <- "Breast Cancer"
seurat_list_BC[["BC12"]]$Type <- "Breast Cancer"

seurat_list_BC[["BN1"]]$Type <- "Breast Control"
seurat_list_BC[["BN2"]]$Type <- "Breast Control"
seurat_list_BC[["BN3"]]$Type <- "Breast Control"
seurat_list_BC[["BN4"]]$Type <- "Breast Control"


seurat_list_BC[["BC1"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list_BC[["BC1"]], pattern = "^MT-")
seurat_list_BC[["BC2"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list_BC[["BC2"]], pattern = "^MT-")
seurat_list_BC[["BC3"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list_BC[["BC3"]], pattern = "^MT-")
seurat_list_BC[["BC4"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list_BC[["BC4"]], pattern = "^MT-")
seurat_list_BC[["BC5"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list_BC[["BC5"]], pattern = "^MT-")
seurat_list_BC[["BC6"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list_BC[["BC6"]], pattern = "^MT-")
seurat_list_BC[["BC7"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list_BC[["BC7"]], pattern = "^MT-")
seurat_list_BC[["BC8"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list_BC[["BC8"]], pattern = "^MT-")
seurat_list_BC[["BC9"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list_BC[["BC9"]], pattern = "^MT-")
seurat_list_BC[["BC10"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list_BC[["BC10"]], pattern = "^MT-")
seurat_list_BC[["BC11"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list_BC[["BC11"]], pattern = "^MT-")
seurat_list_BC[["BC12"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list_BC[["BC12"]], pattern = "^MT-")

seurat_list_BC[["BN1"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list_BC[["BN1"]], pattern = "^MT-")
seurat_list_BC[["BN2"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list_BC[["BN2"]], pattern = "^MT-")
seurat_list_BC[["BN3"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list_BC[["BN3"]], pattern = "^MT-")
seurat_list_BC[["BN4"]][["percent.mt"]]  <- PercentageFeatureSet(seurat_list_BC[["BN4"]], pattern = "^MT-")

VlnPlot(seurat_list_BC[["BC1"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list_BC[["BC2"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list_BC[["BC3"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list_BC[["BC4"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list_BC[["BC5"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list_BC[["BC6"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list_BC[["BC7"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list_BC[["BC8"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list_BC[["BC9"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list_BC[["BC10"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list_BC[["BC11"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list_BC[["BC12"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)

VlnPlot(seurat_list_BC[["BN1"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list_BC[["BN2"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list_BC[["BN3"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
VlnPlot(seurat_list_BC[["BN4"]], features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)


#################
###MERGE ALL####
################

BC_merged <- merge(x = seurat_list_BC[[1]], y = seurat_list_BC[-1])
CC_merged <- merge(x = seurat_list_CC[[1]], y = seurat_list_CC[-1])
OV_merged <- merge(x = seurat_list_OV[[1]], y = seurat_list_OV[-1])
EC_merged <- merge(x = seurat_list_EC[[1]], y = seurat_list_EC[-1])

merged.data <- merge(
  x = BC_merged,
  y = list(OV_merged, CC_merged, EC_merged),
  add.cell.ids = c("BC","OV","CC","EC")
)

# 1. Add sample identifiers
BC_merged$sample <- "BC"
CC_merged$sample <- "CC"
EC_merged$sample <- "EC"
OV_merged$sample <- "OV"

merged.data <- NormalizeData(merged.data, verbose = T)
merged.data <- FindVariableFeatures(merged.data, selection.method = "vst", nfeatures = 2000, verbose = T)
merged.data <- ScaleData(merged.data, verbose = T)
merged.data <- RunPCA(merged.data, npcs = 50, verbose = T)
ElbowPlot(merged.data)
merged.data <- RunUMAP(merged.data, reduction = "pca", dims = 1:30, verbose = T)

DimPlot(merged.data,reduction = "umap") + plot_annotation(title = "Before integration")

merged.data <- merged.data %>% RunHarmony("orig.ident", plot_convergence = T)
harmony_embeddings <- Embeddings(merged.data, 'harmony')
harmony_embeddings[1:5, 1:5]
p1 <- DimPlot(object = merged.data, reduction = "harmony", pt.size = .1, group.by = "orig.ident") + NoLegend()
p2 <- VlnPlot(object = merged.data, features = "harmony_1", group.by = "orig.ident", pt.size = .1) + NoLegend()
library(cowplot)
plot_grid(p1, p2)
ElbowPlot(merged.data)

merged.data <- merged.data %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:30) %>% 
  FindClusters() %>% 
  identity()

merged.data <- SetIdent(merged.data,value = "orig.ident")
DimPlot(merged.data,reduction = "umap") + plot_annotation(title = "After integration (Harmony)")

DimPlot(merged.data, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident') + NoLegend()

merged.data <- SetIdent(merged.data,value = "seurat_clusters")

merged.data <- FindClusters(merged.data, verbose = F, algorithm = 4, resolution = 1)

DimPlot(merged.data,label = T) 
source("./data/cellgeni/custom_seurat_functions.R")
plot_integrated_clusters(merged.data)

###########################
#######join################
###########################

merged.data <- JoinLayers(merged.data)
