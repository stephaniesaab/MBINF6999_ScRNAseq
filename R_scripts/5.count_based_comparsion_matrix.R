#5.count_based_comparsion_matrix.R

#Set output paths---------------------------------------------------------------

.libPaths(c("/home/ssaab/links/projects/def-itobias/R_libs/4.4.0", .libPaths()))
.libPaths()

out_dir <- "~/links/scratch/data/"
#Load libraries ----------------------------------------------------------------
library(SingleR)
library(BiocParallel)
library(Seurat)
library(SingleCellExperiment)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
#Load objects-------------------------------------------------------------------

#Data Seurat object
srt_rds <- readRDS("~/links/scratch/data/merged_labels_seurat_object_lean.rds")

#SingleR object
singleR_rds <- readRDS("~/links/scratch/data/singleR_fulldata_results_V2.rds")

#Making matrix -----------------------------------------------------------------

#Matrix has columns of SingleR labels
#Matrix has rows of manual labels
labels_count_mat <- as.matrix(
  
)