#4.making_loupe_objects.R


#Set output paths---------------------------------------------------------------

.libPaths(c("/home/ssaab/links/projects/def-itobias/R_libs/4.4.0", .libPaths()))
.libPaths()

out_dir <- "~/links/scratch/data/"
#Load libraries ----------------------------------------------------------------
library(SingleR)
library(BiocParallel)
library(Seurat)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(GenomicRanges)
library(dplyr)
library(loupeR)

#Load objects--------------------------------------------

#OG seurat object (no singleR labels)
#srt_rds <- readRDS("~/links/scratch/data/merged_tumour_only.rds")

#seurat object with singleR labels
#srt_rds_labeled <- readRDS("~/links/scratch/data/merged_labels_seurat_object_full.rds")

#seurat object with singleR labels after editing
srt_rds_lean <- readRDS("~/links/scratch/data/merged_labels_seurat_object_lean.rds")

#Making Loupe object--------------------------------------------------
#create_loupe_from_seurat(srt_rds, 
 #                        output_dir = out_dir,
  #                       output_name = "loupe_from_seurat_no_labels",
#			 force = TRUE
 #                        )
#create_loupe_from_seurat(srt_rds_labeled, 
 #                        output_dir = out_dir,
  #                       output_name = "loupe_from_seurat_labeled"
#			)
create_loupe_from_seurat(srt_rds_lean,
output_dir = out_dir,
output_name = "loupe_from_seurat_lean"
)
