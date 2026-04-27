sum_sweep_v2 <- function (sample, cores = 1){
  message("Parameter sweep...")
  system.time({
    sweep.res.sample <- paramSweep(sample)
  })
  sweep.stats.sample <- summarizeSweep(sweep.res.sample, GT = FALSE)
  return(sweep.stats.sample)
}

doubFinder_v2 <-
  function(sample, sweep.stats.sample, rm.doubs = FALSE) {
    bcmvn_sample <- find.pK(sweep.stats.sample) # save this plot
    # 
    message("Selecting pK Value...")
    # alternative method
    maxima <- which(diff(sign(diff(
      bcmvn_sample$BCmetric
    ))) == 2)
    if (length(maxima) == 0) {
      maxima <- which.max(bcmvn_sample$BCmetric)
    } else {
      maxima <-
        maxima[which(bcmvn_sample$BCmetric[maxima] > min(bcmvn_sample$BCmetric))]
    }
    pK_val <- as.numeric(levels(bcmvn_sample$pK)[maxima[1]])
    # 
    # 
    # # peaks <-
    # #   bcmvn_sample[bcmvn_sample$BCmetric > (max(bcmvn_sample$BCmetric) / 2),]
    # 
    # # # consider adding trycatch
    # # # https://stackoverflow.com/questions/12193779/how-to-write-trycatch-in-r
    # 
    # # # pay special attention to the different types of boolean operators in R
    # # # https://medium.com/biosyntax/single-or-double-and-operator-and-or-operator-in-r-442f00332d5b
    # # for (v in 1:nrow(peaks)) {
    # #   if (nrow(peaks) == 1) {
    # #     pK_val <- as.numeric(as.vector(peaks[v, ]$pK))
    # #     break
    # #   }
    # #   else if (peaks[v, ]$BCmetric > peaks[v + 1, ]$BCmetric |
    # #            is.na(peaks[v + 1, ]$BCmetric)) {
    # #     pK_val <- as.numeric(as.vector(peaks[v, ]$pK))
    # #     break
    # #   }
    # #   else {
    # #     next
    # #   }
    # # }
    # 
    message("Homotypic Doublet Proportion Estimation")
    homotypic.prop <-
      modelHomotypic(sample@meta.data$SCT_snn_res.0.5) ## ex: annotations <- sample@meta.data$ClusteringResults
    assumed_db <- (0.0008*ncol(sample) - 0.067)/100
    nExp_poi <-
      round(assumed_db * nrow(sample@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    # 
    # ## Run DoubletFinder with varying classification stringencies
    # ## Remember to update with appropriate params
    # # options(future.globals.maxSize = 2097152000)
    # 
    message("Running doubletFinder...")
    
    if(packageVersion("DoubletFinder") == '2.0.4'){
      # stop("doubFinder_v2 has detected DoubletFinder package version 2.0.4. Install DoubletFinder 2.0.6 in order to use doubFinder_v2 helper function.")
      sample <-
        doubletFinder(
          sample,
          PCs = 1:10,
          pN = 0.25,
          pK = pK_val,
          nExp = nExp_poi,
          reuse.pANN = F,
          sct = FALSE
        )
      
      
      # grep for the column name that has pANN at the start
      pANN_col <-
        colnames(sample@meta.data)[grep("pANN*", colnames(sample@meta.data))]
      
      sample <-
        doubletFinder(
          sample,
          PCs = 1:10,
          pN = 0.25,
          pK = pK_val,
          nExp = nExp_poi.adj,
          reuse.pANN = pANN_col, 
          sct = FALSE
        )
      
    } else {
      print("DoubletFinder package version: ")
      print(packageVersion("DoubletFinder"))
      sample <-
        doubletFinder(
          sample,
          PCs = 1:10,
          pN = 0.25,
          pK = pK_val,
          nExp = nExp_poi,
          reuse.pANN = NULL,
          sct = FALSE
        )
      
      
      # grep for the column name that has pANN at the start
      pANN_col <-
        colnames(sample@meta.data)[grep("pANN*", colnames(sample@meta.data))]
      
      sample <-
        doubletFinder(
          sample,
          PCs = 1:10,
          pN = 0.25,
          pK = pK_val,
          nExp = nExp_poi.adj,
          reuse.pANN = pANN_col, 
          sct = FALSE
        )
    }
    
    
    class_col <-
      colnames(sample@meta.data)[grep("*classification*", colnames(sample@meta.data))]
    
    class_col <- class_col[length(class_col)]
    
    colnames(sample@meta.data)[ncol(sample@meta.data)] <- "doublet"
    
    if (rm.doubs) {
      sample <- subset(sample, subset = doublet == "Singlet")
    }
    
    return(sample)
    
  }
