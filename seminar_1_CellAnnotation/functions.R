get_percent <- function(seuratObj) {
  seuratObj[["percent.mt"]] <-
    Seurat::PercentageFeatureSet(object = seuratObj,
                                 pattern = "^MT-",
                                 assay = "RNA")
  
  seuratObj[["percent.rps"]] <-
    Seurat::PercentageFeatureSet(object = seuratObj,
                                 pattern = "^RP[LS]|^MRPL",
                                 assay = "RNA")
  
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
           rps.cutoff = NULL,
           cc_adjust = FALSE) {
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
      DefaultAssay(seuratObj) <- "RNA"
    }
    
    # Absolutely imperative that these thresholds are set PER SAMPLE
    if(!is.null(mito.cutoff)){
      seuratObj <-
        subset(
          seuratObj,
          subset = nFeature_RNA > feature.cutoff[1] & nFeature_RNA < feature.cutoff[2] &
            nCount_RNA < count.cutoff[2] & nCount_RNA > count.cutoff[1] &
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
            vst.flavor = "v2",
            vars.to.regress = c("S.Score", "G2M.Score")
          )
      }
      else{
        seuratObj <-
          SCTransform(
            seuratObj,
            ncells = 0.25 * ncol(seuratObj),
            vst.flavor = "v2",
            vars.to.regress = c("S.Score", "G2M.Score")
          )
      }
    }
    else {
      seuratObj <- SCTransform(seuratObj, ncells = 0.25 * ncol(seuratObj), vst.flavor = "v2")
    }
    seuratObj <- RunPCA(seuratObj)
    
    eb <- ElbowPlot(seuratObj, reduction = "pca")
    
    seuratObj <-
      FindNeighbors(seuratObj, dims = 1:10) # we pick this based on elbow plot
    seuratObj <- FindClusters(seuratObj, resolution = 0.5)
    seuratObj <- RunUMAP(seuratObj, dims = 1:10)
    return(seuratObj)
  }
