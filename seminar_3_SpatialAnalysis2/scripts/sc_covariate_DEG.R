# tradeSeq analysis using latent time discovered from scVelo (all samples)
#
# OBJECTIVES:
#   - Test differnetial expression along latent time, split analysis among terminal states

library(Seurat)
library(ggplot2)
library(reshape2)
# library(tradeSeq)
# library(RColorBrewer)
# library(SingleCellExperiment)
# library(monocle)
# library(slingshot)
library(future)
library(parallel)
library(dplyr)
library(tidyr)
library(stats)
library(MAST)
library(emmeans)
plan(multisession, workers = 8)
options(future.globals.maxSize= 50000*1024^2)

setwd("/mnt/isilon/cscb/Projects/10x/chous/SCTC-KT-67")
seu <- readRDS("/mnt/isilon/cscb/Projects/10x/chous/SCTC-KT-67/processed/seurat/merged_dat_noCCregression_scaled.RDS")

# metadata from scvelo
meta <- read.csv("code/04-scVelo/scVelo_metadata.csv")
rownames(meta) <- meta$orig.bc
meta <- meta[Cells(seu),]

seu$latent_time <- meta$latent_time
seu$genome <- meta$genome
seu$gata1 <- meta$GATA1
seu$combo <- paste(seu$genome, seu$gata1, sep = "-")

#########
# Trying seurat's FindMarkers, with latent time
Idents(seu) <- "cell_type"

test.genes <- row.names(seu@assays$SCT@scale.data)
test.genes <- test.genes[!grepl("^RP[LS]|^MRPL", test.genes)]
test.genes <- test.genes[!grepl("^MT-", test.genes)]

DefaultAssay(seu) <- "RNA"



# # have to move some data around for 'SCT' to work...
# DefaultAssay(seu) <- "SCT"
# seu[["newSCT"]] <- CreateAssayObject(counts = GetAssayData(seu, assay="SCT", slot="data"))
# seu <- SetAssayData(seu, assay = "newSCT", slot = "scale.data", new.data = GetAssayData(seu, assay="SCT", slot="scale.data")); gc()
# DefaultAssay(seu) <- "newSCT"
# seu <- DietSeurat(seu, assays = c("newSCT","RNA"), counts=T, data=T, scale.data=T)
# seu <- RenameAssays(seu, newSCT = 'SCT')
## Still didn't work, see below to run MAST separately

# seu <- PrepSCTFindMarkers(seu)
DefaultAssay(seu) <- "RNA"
t21.deg <- FindMarkers(seu,
                       test.use = "poisson",
                       group.by = "combo",
                       ident.1 = "T21-WT",
                       ident.2 = "T21-S",
                       subset.ident = "early_Ery",
                       latent.vars = "latent_time",
                       features = test.genes)

# This isn't working...going to manually run MAST
cellClust <- unique(seu$cell_type)
cellClust <- cellClust[cellClust!="early_Ery"]

DefaultAssay(seu) <- "SCT"
Idents(seu) <- "cell_type"

res <- list()

for(ct in cellClust){
  print(ct)
  test.dat <- subset(seu, idents = ct)
  cdat <- data.frame(list(genome=test.dat$genome, gata = test.dat$gata1, latentvar = test.dat$latent_time))
  cdat$wellKey <- Cells(test.dat)
  fdat <- data.frame(list(primerid=test.genes), row.names = test.genes)
  
  sca <- FromMatrix(
    exprsArray = as.matrix(GetAssayData(test.dat)[test.genes,]),
    check_sanity = FALSE,
    cData = cdat,
    fData = fdat
  )
  cond <- factor(x = colData(sca)$genome)
  cond <- relevel(x = cond, ref = "Euploid")
  colData(sca)$varGenome <- cond
  cond <- factor(x = colData(sca)$gata)
  cond <- relevel(x = cond, ref = "WT")
  colData(sca)$varGata <- cond
  
  fmla <- as.formula("~ varGenome*varGata + latentvar")
  zlmCond <- zlm(formula = fmla, sca = sca, parallel=TRUE)
  summaryCond <- summary(object = zlmCond, doLRT = 'varGenomeT21:varGataS', parallel=TRUE)
  summaryH <- summaryCond$datatable %>% filter(component == "H") %>% select(primerid, `Pr(>Chisq)`)
  summaryFC <- summaryCond$datatable %>% filter(component == "logFC") %>% 
    select(primerid, contrast, coef) %>% pivot_wider(names_from = contrast, values_from = coef)
  
  summaryDt <- left_join(summaryH, summaryFC, by=join_by(primerid)) %>% arrange(`Pr(>Chisq)`)
  write.csv(summaryDt, paste0("outs/Seurat/documents/MAST_DEG_itxn_",ct,".csv"))
  res[[ct]] <- summaryDt
}

summaryCond <- summary(object = zlmCond)


res$early_Ery <- summary_earlyery
saveRDS(res, "code/04-tradeSeq/MAST_itxn_res.RDS")

#########
# Trying seurat's FindAllMarkers, no latent time

DefaultAssay(seu) <- "RNA"
Idents(seu) <- "cell_type"
res_seurat <- list()
cellClust <- unique(seu$cell_type)
for(ct in cellClust){
  test.dat <- subset(seu, idents = ct)
  Idents(test.dat) <- "combo"
  res_seurat[[ct]] <- FindAllMarkers(test.dat, assay = 'SCT', features = test.genes)
}

saveRDS(res_seurat, "code/04-tradeSeq/seurat_ct_FindAllMarkers_res.RDS")

res <- readRDS("code/04-tradeSeq/seurat_ct_FindAllMarkers_res.RDS")

for(r in seq_along(res)){
  write.csv(res[[r]], paste0("outs/Seurat/documents/Wilcox_DEG_",names(res)[r],".csv"))
}

#########
# Some plotting
res_ee <- res$early_Ery
test_genes <- res_ee[res_ee$cluster=="T21-WT" & res_ee$avg_log2FC>0 & res_ee$p_val_adj < 0.1,]
test_genes <- test_genes[order(test_genes$avg_log2FC, decreasing = T),]
test_genes <- unlist(sapply(head(rownames(test_genes),5), function(x) strsplit(x,"\\.")[[1]][1]))
VlnPlot(seu, idents = "early_Ery", features = test_genes, 
        group.by = "combo", assay='RNA', pt.size = 0)

res.mast <- readRDS("code/04-tradeSeq/MAST_itxn_res.RDS")
res_ee.mast <- res.mast$early_Ery
res_ee.mast$p.adj <- p.adjust(res_ee.mast$`Pr(>Chisq)`, method="bonferroni")
res_ee.mast <- res_ee.mast[res_ee.mast$p.adj<0.1,]

test_genes <- res_ee.mast[res_ee.mast$`varGenomeT21:varGataS`>0,]
test_genes <- test_genes[order(test_genes$p.adj, decreasing = F),]

test_genes <- res_ee.mast[res_ee.mast$varGataS>0,]
# test_genes <- res_ee.mast[res_ee.mast$p.adj < 1e-10,]
test_genes <- test_genes[order(test_genes$varGataS, decreasing = T),]

DefaultAssay(seu) <- "SCT"
test_dat <- subset(seu, idents="early_Ery")
pdata <- GetAssayData(test_dat)

pp <- data.frame(list(latent_time=test_dat$latent_time,
                      combo=test_dat$combo))
g <- c("AGAP1", "APEX1", "APRT", "CD24", "CYTOR", "AUTS2")
pp <- cbind(pp, t(pdata[g,]))
# pp <- cbind(pp, t(pdata[test_genes$primerid[7:12],]))

pp <- melt(pp, id.vars = c("latent_time","combo"))
ggplot(data=pp, aes(x=latent_time, y=value, colour=combo))+
  geom_smooth(se=F) + facet_wrap('variable', ncol = 2, scales = 'free')

## Single gene
pp <- data.frame(list(latent_time=test_dat$latent_time,
                      combo=test_dat$combo))
pp$gene <- pdata["COL24A1",]
ggplot(data=pp, aes(x=latent_time, y=gene, colour=combo))+
  geom_smooth(se=F)

VlnPlot(seu, idents = "early_Ery", features = test_genes$primerid[1], 
        group.by = "combo", assay='SCT', pt.size = 0)

VlnPlot(test.dat, features = head(summaryDt$primerid), group.by = "combo", assay = "SCT")
VlnPlot(seu, features = "HSPA8", group.by = "time", idents = "Ery", 
        assay = "RNA", split.by = 'combo', slot = 'counts')
VlnPlot(seu, features = "HSPA8", group.by = "gata1", idents = "Ery", 
        assay = "SCT", slot = 'data')
VlnPlot(seu, features = head(row.names(res_seurat$early_Ery)), idents = "early_Ery", group.by = "combo", assay = "SCT")
FeaturePlot(seu, features = head(rownames(t21.deg)), raster=F)
FeaturePlot(seu, features = "MALAT1", raster=F)

Idents(plot.data) <- "combo"
plot.data <- subset(plot.data, idents = c("T21-WT","T21-S"))



#cds <- as.CellDataSet(seu)
sce <- as.SingleCellExperiment(seu)

#########
# Predict lineages using Monocle
#   This is annoying and conversion from seurat isn't working well
# set.seed(200)
# 
# sct.counts <- GetAssayData(seu)
# meta <- new("AnnotatedDataFrame",seu@meta.data)
# varmeta <- new("AnnotatedDataFrame",data.frame(list(gene_short_name=rownames(seu)), row.names = rownames(seu)))
# cds <- newCellDataSet(as.sparse(sct.counts), phenoData = meta, featureData = varmeta)
# 
# # cds <- newCellDataSet(counts, phenoData = new("AnnotatedDataFrame", data = pd),
# #                       featureData = new("AnnotatedDataFrame", data = fd))
# # cds <- estimateSizeFactors(cds)
# # cds <- reduceDimension(cds, max_components = 2)

#########
# Predict lineages using SlingShot
sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP')
lin1 <- getLineages(sce, sce$seurat_clusters, start.clus = 1, reducedDim = 'UMAP')
sling <- SlingshotDataSet(sce)

plotGeneCount(sling, sce, clusters=sce$seurat_clusters)

seu$lin_time <- sce$slingPseudotime_1
FeaturePlot(seu, features = "lin_time", raster = F)

####
# TradeSeq steps
set.seed(3)
icMat <- evaluateK(counts = GetAssayData(seu),
                   sds = sling,
                   nGenes = 300,
                   k = 3:7,
                   BPPARAM = BiocParallel::bpparam(4))

set.seed(7)
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
              nknots = 6, verbose = FALSE)

save.image("checkpoint.RData")
