get_percent <- function(seuratObj, species = "human", exp_type = "GEX") {
  assay_name <- "RNA"
  if(exp_type == "spatial"){
    assay_name <- "Spatial"
  }
  
  if(species == "human"){
    seuratObj[["percent.mt"]] <-
      Seurat::PercentageFeatureSet(object = seuratObj,
                                   pattern = "^MT-",
                                   assay = assay_name)
    
    seuratObj[["percent.rps"]] <-
      Seurat::PercentageFeatureSet(object = seuratObj,
                                   pattern = "^RP[LS]|^MRPL",
                                   assay = assay_name)
  } else if(species == "mouse"){
    seuratObj[["percent.mt"]] <-
      Seurat::PercentageFeatureSet(object = seuratObj,
                                   pattern = "^Mrp",
                                   assay = assay_name)
    
    seuratObj[["percent.rps"]] <-
      Seurat::PercentageFeatureSet(object = seuratObj,
                                   pattern = "^Rp[ls]",
                                   assay = assay_name)
  } else{
    seuratObj[["percent.mt"]] <- 0
    seuratObj[["percent.rps"]] <- 0
  }
  
  return(seuratObj)
}

LoadPackages_dgs <- function (packages, reinstall = FALSE, list.source = "https://raw.githubusercontent.com/zhezhangsh/RoCA/master/document/packages.txt") 
{
  require(devtools)
  require(RCurl)
  if (!reinstall) {
    installed <- installed.packages()
    pkgs <- packages[!(packages %in% rownames(installed))]
  }
  else pkgs <- packages
  if (length(pkgs) > 0) {
    known <- read.table(list.source, stringsAsFactors = FALSE)
    pkgs <- pkgs[!sapply(pkgs, function(p) {
      ind <- which(known[, 2] == p)[1]
      if (is.na(ind)) 
        FALSE
      else {
        cat("Installing package: ", p, "\n")
        if (known[ind, 1] == "github") 
          install_github(known[ind, 3], quiet = TRUE)
        else install_url(known[ind, 3], quite = TRUE)
        TRUE
      }
    })]
  }
  if (length(pkgs) > 0) {
    try(install.packages(pkgs))
    pkgs <- pkgs[!(pkgs %in% rownames(installed.packages()))]
    if (length(pkgs) > 0) {
      # source("http://bioconductor.org/biocLite.R")
      # biocLite(pkgs, type = "source", suppressUpdates = TRUE)
      BiocManager::install(pkgs)
      pkgs <- pkgs[!(pkgs %in% rownames(installed.packages()))]
      if (length(pkgs) > 0) {
        i <- try(sapply(pkgs, function(p) install_github(paste("https://github.com", 
                                                               p, sep = "/"))))
      }
    }
  }
  loaded <- sapply(packages, function(p) require(p, character.only = TRUE))
  loaded
}

seurat_process_v2 <-
  function(counts,
           counts_name,
           feature.cutoff,
           count.cutoff,
           # can we store these options in the seurat object?
           mito.cutoff = NULL,
           rps.cutoff = 100,
           cc_adjust = FALSE,
           sct.med = NULL) {
    if (class(counts)[1] != "Seurat") {
      seuratObj <- CreateSeuratObject(
        counts = counts,
        project = counts_name,
        min.cells = 10,
        min.features = 200
      )
    } else {
      # if `counts` is already a Seurat object:
      seuratObj <- counts
      DefaultAssay(seuratObj) <- "Spatial"
    }
    
    # Absolutely imperative that these thresholds are set PER SAMPLE
    if(!is.null(mito.cutoff)){
      seuratObj <-
        subset(
          seuratObj,
          subset = nFeature_Spatial > feature.cutoff[1] & nFeature_Spatial < feature.cutoff[2] &
            nCount_Spatial < count.cutoff[2] & nCount_Spatial > count.cutoff[1] &
            percent.mt < mito.cutoff &
            percent.rps < rps.cutoff
        )
    }
    # seuratObj@assays$RNA@misc <-
    #   data.frame (
    #     parameter  = c("Count Cutoff", "Mito Cutoff", "Ribo Cutoff", "CC Adjust?"),
    #     value = c(count.cutoff, mito.cutoff, rps.cutoff, cc_adjust)
    #   )
    
    # normalize
    
    # TODO expanding limit to prevent some arbitrary error msg
    # options(future.globals.maxSize = 20971520000)
    if (cc_adjust) {
      seuratObj <-
        CellCycleScoring(
          seuratObj,
          s.features = cc.genes$s.genes,
          g2m.features = cc.genes$g2m.genes,
          set.ident = TRUE
        )
      seuratObj <-
        FindVariableFeatures(seuratObj,
                             selection.method = "vst",
                             nfeatures = 2000)
      if (dim(seuratObj)[2] < 10000) {
        seuratObj <-
          SCTransform(
            seuratObj,
            ncells = 2000,
            vars.to.regress = c("S.Score", "G2M.Score")
          )
      }
      else{
        seuratObj <-
          SCTransform(
            seuratObj,
            ncells = 0.25 * ncol(seuratObj),
            vars.to.regress = c("S.Score", "G2M.Score")
          )
      }
    }
    else {
      if(!is.null(sct.med)){
        seuratObj <- SCTransform(seuratObj, assay = "Spatial", vst.flavor = "v2", scale_factor = sct.med)
      }
      else{
        seuratObj <- SCTransform(seuratObj, assay = "Spatial", vst.flavor = "v2")
      }
    }
    seuratObj <- RunPCA(seuratObj)
    
    eb <- ElbowPlot(seuratObj, reduction = "pca")
    
    seuratObj <-
      FindNeighbors(seuratObj, dims = 1:10) # we pick this based on elbow plot
    seuratObj <- FindClusters(seuratObj, resolution = 0.5)
    seuratObj <- RunUMAP(seuratObj, dims = 1:10)
    return(seuratObj)
  }

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