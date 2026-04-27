# helper functions



load10X_dev = function(dataDir,cellIDs=NULL,channelName=NULL,readArgs=list(),includeFeatures=c('Gene Expression'),verbose=TRUE,...){
  #Work out which version we're dealing with
  isV3 = dir.exists(file.path(dataDir,'raw_feature_bc_matrix'))
  isV7 = dir.exists(file.path(dataDir,'analysis','clustering','gene_expression_graphclust'))
  isMulti = dir.exists(file.path(dataDir,'per_sample_outs'))
  tgt = file.path(dataDir,
                  ifelse(isV3,'raw_feature_bc_matrix','raw_gene_bc_matrices'))
  #Add the reference genome for the non-V3 ones
  if(!isV3){
    sname <- unlist(str_split(dataDir,"/"))
    sname <- sname[length(sname)-1]
    tgt = file.path(dataDir,"multi","count","raw_feature_bc_matrix")
  }
  if(verbose)
    message(sprintf("Loading raw count data"))
  dat = do.call(Read10X,c(list(data.dir=tgt),readArgs))
  if(verbose)
    message(sprintf("Loading cell-only count data"))
  if(!is.null(cellIDs)){
    #This is now redundant as we require the version of Seurat that does not strip the suffix
    ####Do the same sripping that Seurat does on IDs
    ###if(all(grepl('\\-1$',cellIDs)))
    ###  cellIDs = gsub('\\-1$','',cellIDs)
    #Check we have the IDs
    if(!all(cellIDs %in% colnames(dat)))
      stop("Not all supplied cellIDs found in raw data.")
    datCells = dat[,match(cellIDs,colnames(dat))]
  }else{
    #Work out which ones contain cells
    tgt = file.path(dataDir,
                    ifelse(isV3,'filtered_feature_bc_matrix','filtered_gene_bc_matrices'))
    if(!isV3)
      tgt = file.path(dataDir,"per_sample_outs", sname, "count", "sample_filtered_feature_bc_matrix")
    datCells = do.call(Read10X,c(list(data.dir=tgt),readArgs))
    #If it's a list of multiple types, have to decide what to include and collapse to one matrix.
    if(is.list(dat)){
      dat = do.call(rbind,dat[includeFeatures])
      datCells = do.call(rbind,datCells[includeFeatures])
    }
  }
  if(verbose)
    message(sprintf("Loading extra analysis data where available"))
  #Get the cluster annotation if available
  mDat = NULL
  #What needs to be added to make V7 directory structure work
  v7Prefix=ifelse(isV7,'gene_expression_','')
  tgt = ifelse(isMulti,
               file.path(dataDir,"per_sample_outs", sname, "count",'analysis','clustering', 'gene_expression_graphclust','clusters.csv'),
               file.path(dataDir,'analysis','clustering','graphclust','clusters.csv')
  )
  if(file.exists(tgt)){
    clusters = read.csv(tgt)
    mDat = data.frame(clusters=clusters$Cluster,row.names=clusters$Barcode)
  }
  #Add fine grained clusters too if present
  tgt = ifelse(isMulti,
               file.path(dataDir,"per_sample_outs", sname, "count",'analysis','clustering', 'gene_expression_kmeans_10_clusters','clusters.csv'),
               file.path(dataDir,'analysis','clustering','kmeans_10_clusters','clusters.csv')
  )
  if(file.exists(tgt)){
    clusters = read.csv(tgt)
    mDat$clustersFine = clusters$Cluster
  }
  #Get tSNE if available and point to it
  tgt = ifelse(isMulti,
               file.path(dataDir,"per_sample_outs", sname, "count",'analysis','tsne','gene_expression_2_components','projection.csv'),
               file.path(dataDir,'analysis','tsne','2_components','projection.csv')
  )
  if(file.exists(tgt)){
    tsne = read.csv(tgt)
    if(is.null(mDat)){
      mDat = data.frame(tSNE1=tsne$TSNE.1,tSNE2=tsne$TSNE.2,row.names=tsne$Barcode)
    }else{
      mDat$tSNE1 = tsne$TSNE.1[match(rownames(mDat),tsne$Barcode)]
      mDat$tSNE2 = tsne$TSNE.2[match(rownames(mDat),tsne$Barcode)]
    }
    DR = c('tSNE1','tSNE2')
  }else{
    DR=NULL
  }
  #Ensure rownames of metadata match column names of counts
  if(!is.null(mDat) && any(rownames(mDat)!=colnames(datCells))){
    rownames(mDat) = gsub('-1$','',rownames(mDat))
    if(any(rownames(mDat)!=colnames(datCells)))
      stop("Error matching meta-data to cell names.")
  }
  #Get a name for the channel
  if(is.null(channelName))
    channelName = ifelse(is.null(names(dataDir)),dataDir,names(dataDir))
  channel = SoupChannel(tod = dat,
                        toc = datCells,
                        metaData = mDat,
                        channelName = channelName,
                        dataDir = dataDir,
                        dataType='10X',
                        isV3=isV3,
                        DR=DR,
                        ...
  )
  return(channel)
}
