## Refining cluster annotations based on known markers
# The script was getting a little crowded

library(Seurat)
library(SeuratDisk)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(stringr)
library(here)
library(viridisLite)
library(readxl)
library(reshape2)
library(pheatmap)
library(tibble)
library(data.table)
library(future)
plan("multisession", workers = 4)
options(future.globals.maxSize= 20000*1024^2)

setwd("/mnt/isilon/cscb/Projects/Merscope/marco")

seu <- readRDS("processed/seu_merge_annotated_filtered.rds")
seu@meta.data <- seu@meta.data %>% select(!starts_with("Cluster"))

'''

A bunch of other code....

'''


# EXPORT
seu.out <- readRDS("processed/seu_merge_annotated_filtered.rds")
seu.out@meta.data <- seu.out@meta.data %>% 
  select(orig.ident, nCount_RNA, nFeature_RNA, center_x, center_y, perimeter_area_ratio, 
        anisotropy, solidity, Rab11b.Rat_raw, Rab11b.Rat_high_pass,
        Sun1GFP_raw, Sun1GFP_high_pass, EEA1.Chicken_raw, EEA1.Chicken_high_pass,
        Lamp1.Goat_raw, Lamp1.Goat_high_pass, CRE_raw, CRE_high_pass, Cas9_raw,
        Cas9_high_pass, seurat_clusters, predicted.id_full, predicted.id.e12_full,
        CRE_pos, ref_ct_l1.full, ref_ct_l2.full, cluster_full, anatomy) %>%
  mutate_if(is.factor, as.character)

space_coords <- seu.out@meta.data %>% 
  dplyr::select(center_x, center_y) %>%
  mutate_all(as.integer) %>%
  as.matrix
colnames(space_coords) <- c("space1","space2")
space_coords <- CreateDimReducObject(space_coords, key = "space_")
seu.out@reductions[["space"]] <- space_coords

Idents(seu.out) <- "seurat_clusters"
DefaultAssay(seu.out) <- "RNA"
seu.out.slim <- DietSeurat(seu.out, assays=c("RNA"), dimreducs = c("umap"))

seu.out.slim[["RNA"]] <- as(object = seu.out.slim[["RNA"]], Class = "Assay")
SaveH5Seurat(seu.out.slim, filename = paste0("outs/cellxgene/merged_data.h5Seurat"), overwrite = T)
Convert(paste0("outs/cellxgene/merged_data.h5Seurat"), dest = "h5ad", overwrite = T)

seu_objs <- SplitObject(seu.out, split.by = "orig.ident")
# individual slides, physical space
for(x in seu_objs){
  test <- x
  sname <- test$orig.ident[1]
  SaveH5Seurat(test, filename = paste0("outs/cellxgene/",sname,".h5Seurat"), overwrite = T)
  Convert(paste0("outs/cellxgene/",sname,".h5Seurat"), dest = "h5ad", overwrite = T)
}


##############
# For Loupe
library(loupeR)
make_cloupe_DGS <- function (obj, output_dir = NULL, output_name = NULL, dedup_clusters = FALSE, 
                             feature_ids = NULL, executable_path = NULL, force = FALSE) 
{
  # v <- needs_setup(executable_path)
  # if (!v$success) {
  #   stop(v$msg)
  # }
  if (!is(obj, "Seurat")) {
    stop(validation_err("input object was not a Seurat object", 
                        "Seurat Object"))
  }
  print("extracting matrix, clusters, and projections")
  namedAssay <- select_assay_DGS(obj)
  if (is.null(namedAssay)) {
    stop(validation_err("could not find a usable count matrix", 
                        "Seurat Object"))
  }
  assay_name <- names(namedAssay)
  assay <- namedAssay[[1]]
  counts <- assay$counts
  clusters <- select_clusters(obj, dedup = dedup_clusters)
  projections <- select_projections(obj)
  print(paste0("selected assay:", assay_name))
  print(paste0("selected clusters:", names(clusters)))
  print(paste0("selected projections:", names(projections)))
  seurat_obj_version <- NULL
  if (!is.null(obj@version)) {
    seurat_obj_version <- as.character(obj@version)
  }
  success <- create_loupe(counts, clusters = clusters, projections = projections, 
                          output_dir = output_dir, output_name = output_name, feature_ids = feature_ids, 
                          executable_path = executable_path, force = force, seurat_obj_version = seurat_obj_version)
  invisible(success)
}

select_assay_DGS <- function (obj) 
{
  assay_priority <- c()
  for (name in Seurat::Assays(obj)) {
    if (identical(name, obj@active.assay)) {
      priority <- 1
    }
    else if (grepl("rna", name, ignore.case = TRUE)) {
      priority <- 2
    }
    else {
      priority <- 3
    }
    assay_priority[name] <- priority
  }
  assay_priority <- sort(assay_priority)
  assay <- NULL
  for (i in seq_along(assay_priority)) {
    name <- names(assay_priority[i])
    assay <- Seurat::GetAssay(obj, assay = name)
    counts <- assay$counts
    if (length(counts) > 0) {
      result = list()
      result[[name]] = assay
      return(result)
    }
  }
  NULL
}
#################
# The actual export

test <- seu_objs$run1_reg1
test$CRE_pos <- test$CRE_raw > 5e8

space_coords <- test@meta.data %>% 
  dplyr::select(center_x, center_y) %>%
  mutate_all(as.integer) %>%
  as.matrix
colnames(space_coords) <- c("space1","space2")
space_coords <- CreateDimReducObject(space_coords, key = "space_")

test@reductions[["space"]] <- space_coords
out_lite <- DietSeurat(test, assays = c("RNA"), dimreducs = c("space"))

make_cloupe_DGS(
  out_lite,
  output_dir = "outs/Loupe",
  output_name = "run1_reg1",
  force=T
)










