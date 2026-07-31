#2.comparing_cluster_annotations_hpc.R


.libPaths(c("/home/ssaab/links/projects/def-itobias/R_libs/4.4.0", .libPaths()))
.libPaths()

#Load libraries ----------------------------------------------------------------
library(SingleR)
library(BiocParallel)
library(Seurat)
library(SingleCellExperiment)
library(org.Mm.eg.db)
library(SummarizedExperiment)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(ggridges)
library(plotly)

#Load objects-------------------------------------------------------------------

#SingleR annotated object -> Full dataset
singleR_rds <- readRDS("~/links/scratch/data/singleR_fulldata_results.rds")
cat("the manual clusters are ", unique(singleR_rds$seurat_clusters)) #NULL --> No manual clusters

#Data Seurat object
srt_rds <- readRDS("~/links/scratch/data/merged_tumour_only.rds")
cat("the manual clusters are ", length(unique(srt_rds$seurat_clusters)))

#Adding Manual cluster labels---------------------------------------------------
# Make sure you cover all levels in your seurat object(for me, I had the object "merged_all")

new_cluster_ids <- c(
  "0" = "Upper Layer Cortical Glutamatergic Neuron",
  "1" = "Cycling Early Neuroepithelial Stem Cell",
  "2" = "Mesenchymal Stromal Cell",
  "3" = "Radial Glia Cell",
  "4" = "Immature GABAergic Interneuron",
  "5" = "Primordial Germ Cell",
  "6" = "Immature Neuron",
  "7" = "Deep Layer Cortical Glutamatergic Neuron",
  "8" = "Mature GABAergic Interneuron",
  "9" = "Peptidergic GABA Neuron",
  "10" = "Floor Plate Progenitor Cell 1",
  "11" = "Floor Plate Progenitor Cell 2",
  "12" = "M1-Like Macrophage",
  "13" = "Immature Dopaminergic Neuron",
  "14" = "Immature Autonomic Neuron",
  "15" = "Ventral Diencephalic Neural Progenitor Cell",
  "16" = "Cycling Primordial Germ Cell",
  "17" = "Forebrain Radial Glia Cell",
  "18" = "Anterior Epiblast-Like Cell",
  "19" = "Visceral Endoderm Epithelial Cell",
  "20" = "Mature Astrocyte",
  "21" = "Fibroblasts",
  "22" = "Myelinating Oligodendrocyte",
  "23" = "Luminal/Glandular Epithelial Cell",
  "24" = "Skeletal Myocyte Progenitor Cell",
  "25" = "Vascular Endothelial Cell",
  "26" = "Pericyte",
  "27" = "Type II Skeletal Muscle 1",
  "28" = "Immature Skeletal Myocyte",
  "29" = "Oligodendrocyte Progenitor Cell",
  "30" = "Adventitial Fibroblast Progenitor Cell",
  "31" = "Type II Skeletal Muscle 2",
  "32" = "Motile Multiciliated Epithelial Cell",
  "33" = "M2-Like Macrophage",
  "34" = "Vascular Smooth Muscle Cell",
  "35" = "White/Beige Adipocyte",
  "36" = "Choroid Plexus Epithelial Cell",
  "37" = "Immature Cardiomyocyte",
  "38" = "Osteoblast"
)

#Ensuring seurat sees the identities as the cluster numbers
Idents(srt_rds) <- "seurat_clusters"

#Setting new labels to the seurat object
srt_rds <- RenameIdents(object = srt_rds, new_cluster_ids)

# Store cell type labels in metadata
srt_rds$cell_type <- Idents(srt_rds)

#Adding Metadata----------------------------------------------------------------

#First check if cell barcodes are in the same order
#E.g. cell barcode: WT_REP1_WT_REP1_1A_CTATTGAGTTTGTCCCAGTAGGCT

#Check they're the same order in test and singleR annotations
cell_barcodes_og <- Cells(srt_rds)
cat("length of unique UMI is ", length(unique(cell_barcodes_og)))

cell_barcodes_singleR <- rownames(singleR_rds)
cat("length of the cell barcodes in the SingleR object is ", length(cell_barcodes_singleR))

#Check if they're in the same order
cat("Are the cell barcodes in the same order?", identical(cell_barcodes_og, cell_barcodes_singleR))


#Add the cluster labels to the metadata
singleR_labels <- singleR_rds$labels
srt_rds[["singleR_labels"]] <- singleR_labels

#Save the RDS with singleR labels
saveRDS(srt_rds, "~/links/scratch/data/merged_labels_seurat_object_full.rds")
