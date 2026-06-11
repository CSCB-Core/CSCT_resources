setwd("/mnt/isilon/cscb/Projects/10x_Visium/bennetf/SCTC-GL-156")
library(Seurat)
library(ggalluvial)
library(dplyr)
library(stringr)
library(tidyr)
library(writexl)
library(SPARK)
library("ReactomePA")
library(future)
library(pheatmap)
library(cowplot)
library(org.Mm.eg.db)
library(ggplot2)
library(ggridges)
# library(ComplexHeatmap)
library(viridisLite)
plan("multisession", workers = 8)
options(future.globals.maxSize= 10000*1024^2)

seu <- readRDS("processed/seu_merge_8micron_edit.RDS")

#####
# sankey plot

pdata <- seu@meta.data %>%
  group_by(orig.ident, predicted_ct_3, predicted_ct_2, predicted_ct) %>%
  summarise(n=n())
pdata$exp <- "WT"
pdata$exp[pdata$orig.ident == "2TWI-L-Brain-26_PAS"] <- "Twitch"
ggplot(pdata, aes(y=n, axis1=predicted_ct_3, axis2=predicted_ct_2, axis3=predicted_ct))+
  geom_alluvium(aes(fill=exp))+
  geom_stratum(alpha = .25, width = 1/8, reverse = FALSE)+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE, size=2)+
  theme_minimal()

##########################################
'''
4.	Genes enriched or de-enriched at a range of distances  from globoid cells. Then pathway analysis
  1.	For this, we should concentrate on brainstem
  2. the "range" of distances are the shaped splines

subset to brainstem! (midbrain + pons)

Maybe fit 5 shapes:
1. [increasing] cosine?
2. [decreasing] sine?
3. [peak at various distances] various gaussians

just take top 15% for enrichment

'''

seu.twitcher <- subset(seu, subset = orig.ident == "2TWI-L-Brain-26_PAS")
Idents(seu.twitcher) <- "Anatomy"
seu.twitcher <- subset(seu.twitcher, idents = c("Midbrain","Pons"))
# redo the gene fitting
pas_cells <- seu.twitcher$PAS_pos & seu.twitcher$predicted_ct %in% c("Perivascular macrophages", "Microglia")
pas_coords <- seu.twitcher@images$slice_2TWI.L.Brain.26_PAS$centroids@coords[pas_cells,]
all_coords <- seu.twitcher@images$slice_2TWI.L.Brain.26_PAS$centroids@coords

pas_norm <- apply(pas_coords, 1, function(x) sum(x**2))
all_norm <- apply(all_coords, 1, function(x) sum(x**2))

PAS_dist <- sqrt(outer(pas_norm, all_norm, FUN = "+") - 2*pas_coords%*%t(all_coords))

# filter out the very far cells
PAS_dmin <- apply(PAS_dist, 2, min)
plot(density(PAS_dmin))
keep_cells <- PAS_dmin < 400
# PAS_dist_filter <- PAS_dist[,keep_cells]

seu_test <- subset(seu.twitcher, cells = Cells(seu.twitcher)[keep_cells])
seu_test$PAS_dist <- PAS_dmin[keep_cells]

####
# getting the spatially relevant genes
test.ex <- GetAssayData(seu.twitcher, layer = "counts")
test.coords <- seu.twitcher@images$slice_2TWI.L.Brain.26_PAS$centroids@coords
rownames(test.coords) <- colnames(test.ex)
sparkX <- sparkx(test.ex,
                 data.frame(test.coords),
                 numCores=8,
                 option="mixture")
sum(sparkX$res_mtest$adjustedPval == 0)
spat.varg <- rownames(sparkX$res_mtest)[sparkX$res_mtest$adjustedPval == 0]
# test_genes <- intersect(VariableFeatures(seu.twitcher), spat.varg)

#  correcting for zero inflation
test.ex <- GetAssayData(seu_test, layer = "data")
g.iqr <- future.apply::future_sapply(spat.varg, function(g) quantile(test.ex[g,], c(0.95,0.5)))
g.iqr <- data.frame(t(g.iqr))
# g.iqr$thresh <- apply(g.iqr, 1, function(x) abs(diff(x)))

test_genes <- spat.varg
save(test_genes, g.iqr, seu_test, file="code/09-gary_objectives/outs/histo_DEG_checkpoint.RData")
load("code/09-gary_objectives/outs/histo_DEG_checkpoint.RData")
####
# next step
parm.resolution <- 16 # bin size
parm.distance <- 400 # total distance covered
parm.completeness <- 1 # fraction of available cells from total cells that would be assumed
parm.alpha <- parm.resolution / (parm.distance*parm.completeness)

g.loess_fit <- future.apply::future_sapply(test_genes, function(g) {
  # g.loess_fit <- list()
  # for(g in test_genes){
  #   print(which(test_genes==g)/length(test_genes))
  gdata <- data.frame(list("dist"=seu_test$PAS_dist,
                           "gene"=test.ex[g,])) %>%
    mutate("dist_bins"=cut(dist, breaks = seq(8,400,8))) # for rarefaction
  samp.max <- min(table(gdata$dist_bins))
  gdata <- gdata %>% 
    group_by(dist_bins) %>%
    slice_sample(n=samp.max) %>%
    filter(gene < g.iqr[g,"X95."]*1.5 & gene >= g.iqr[g,"X50."] - g.iqr[g,"X50."]*1.5)
  gdata$dist_bins <- as.numeric(gdata$dist_bins)*8
  
  sum_levels <- gdata %>%
    group_by(dist_bins) %>%
    summarise(totals=sum(gene))
  if(sum(sum_levels$totals==0) > 0.1*nrow(sum_levels)){
    return(rep(NA, length(seq(8,392,8))))
  }else{
    if(samp.max < 100){
      return(rep(NA, length(seq(8,392,8))))
    }
    # gdata$gene <- (gdata$gene/median(gdata$gene))*0.5 # standarize the shape?
    loess.pred <- tryCatch(
      {
        g.loess <- loess(gene~dist_bins, data = gdata, span = 0.2)
        predict(g.loess, newdata = seq(8,392,8))
      },
      error = function(cond){
        rep(NA, length(seq(8,392,8)))
      }
    )
    return(loess.pred)
  }
},
future.seed=T
)

saveRDS(g.loess_fit, "code/09-gary_objectives/outs/brainstem_gene_loess_fit.RDS")
g.loess_fit <- readRDS("code/09-gary_objectives/outs/brainstem_gene_loess_fit.RDS")
empty.genes <- apply(g.loess_fit, 2, function(x) all(is.na(x)))
g.loess_fit <- g.loess_fit[,!empty.genes]
g.loess.scale <- apply(g.loess_fit, 2, function(x) (x-min(x))/(max(x)-min(x)))

g_bckd <- mapIds(org.Mm.eg.db, rownames(seu_test), column = "ENTREZID", keytype = "SYMBOL")
g_bckd <- g_bckd[!is.na(g_bckd)]

####
# running the enrichments

pdf("code/09-gary_objectives/outs/histology_DEA/pattern_summary.pdf",
    width = 10,
    height = 12)

# [decreasing] sine
x_ = seq(0, pi/2, length.out=49)
y_ = sin(x_)
fit.error <- abs(sweep(g.loess.scale, 1, y_, FUN = "-")) %>%
  colSums()
fit.error <- fit.error[order(fit.error)]
p1 <- ggplot(data.frame("dist"=x_, "y_"=y_), aes(dist,y_)) + geom_point() + theme_bw()+ theme(axis.text.x = element_blank())
p2 <- ggplot(data.frame("error"=fit.error), aes(error)) + geom_density() + theme_bw() + geom_vline(xintercept = 15, col="red")
top_genes <- names(fit.error)[fit.error < 15]

hdata <- g.loess.scale[,top_genes]
# p2 <- Heatmap(t(hdata), cluster_columns = F)
p3 <- grid::grid.grabExpr(print(pheatmap(t(hdata), cluster_cols = F, fontsize_row = 6)))

tg_etrz <- mapIds(org.Mm.eg.db, top_genes, column = "ENTREZID", keytype = "SYMBOL")
epath <- enrichPathway(tg_etrz,
                       organism = "mouse",
                       universe = g_bckd,
                       readable = T)
write.csv(epath@result, "code/09-gary_objectives/outs/histology_DEA/sine_pathway_enrichment.csv")
pdata <- epath@result %>% slice_head(n=30)
pdata$gr <- pdata$Count/length(top_genes)
pdata$Description_short <- sapply(pdata$Description, function(s) substring(s,1,35))
pdata$Description_short <- factor(pdata$Description_short, levels = pdata$Description_short[order(pdata$gr)])
p4 <- ggplot(pdata, aes(x=gr, y=Description_short)) + 
  geom_point(aes(size=Count, colour = -log(p.adjust))) +
  theme_bw() + theme(axis.text.y = element_text(size=8))
# p4 <- dotplot(epath, showCategory=30, includeAll=T) + theme(axis.text.y = element_text(size=8))

p_left <- plot_grid(plotlist = list(arrangeGrob(p1,p2, ncol=2),p3), nrow = 2, rel_heights = c(0.2,1))
plot_grid(p_left, p4, ncol = 2)

# [increasing] cosine
x_ = seq(0, pi/2, length.out=49)
y_ = cos(x_)
fit.error <- abs(sweep(g.loess.scale, 1, y_, FUN = "-")) %>%
  colSums()
fit.error <- fit.error[order(fit.error)]
p1 <- ggplot(data.frame("dist"=x_, "y_"=y_), aes(dist,y_)) + geom_point() + theme_bw() + theme(axis.text.x = element_blank())
p2 <- ggplot(data.frame("error"=fit.error), aes(error)) + geom_density() + theme_bw() + geom_vline(xintercept = 15, col="red")
top_genes <- names(fit.error)[fit.error < 15]

hdata <- g.loess.scale[,top_genes]
# p2 <- Heatmap(t(hdata), cluster_columns = F)
p3 <- grid::grid.grabExpr(print(pheatmap(t(hdata), cluster_cols = F, fontsize_row = 6)))

tg_etrz <- mapIds(org.Mm.eg.db, top_genes, column = "ENTREZID", keytype = "SYMBOL")
epath <- enrichPathway(tg_etrz,
                       organism = "mouse",
                       universe = g_bckd,
                       readable = T)
write.csv(epath@result, "code/09-gary_objectives/outs/histology_DEA/cosine_pathway_enrichment.csv")
pdata <- epath@result %>% slice_head(n=30)
pdata$gr <- pdata$Count/length(top_genes)
pdata$Description_short <- sapply(pdata$Description, function(s) substring(s,1,35))
pdata$Description_short <- factor(pdata$Description_short, levels = pdata$Description_short[order(pdata$gr)])
p4 <- ggplot(pdata, aes(x=gr, y=Description_short)) + 
  geom_point(aes(size=Count, colour = -log(p.adjust))) +
  theme_bw() + theme(axis.text.y = element_text(size=8))
# p4 <- dotplot(epath, showCategory=30, includeAll=T) + theme(axis.text.y = element_text(size=8))

p_left <- plot_grid(plotlist = list(arrangeGrob(p1,p2, ncol=2),p3), nrow = 2, rel_heights = c(0.2,1))
plot_grid(p_left, p4, ncol = 2)

# [increasing] normal
x_ = seq(49)
y_ = dnorm(x_, sd= 5, log=F)
y_ <- y_/max(y_)

fit.error <- abs(sweep(g.loess.scale, 1, y_, FUN = "-")) %>%
  colSums()
fit.error <- fit.error[order(fit.error)]
p1 <- ggplot(data.frame("dist"=x_, "y_"=y_), aes(dist,y_)) + geom_point() + theme_bw() + theme(axis.text.x = element_blank())
p2 <- ggplot(data.frame("error"=fit.error), aes(error)) + geom_density() + theme_bw() + geom_vline(xintercept = 15, col="red")
top_genes <- names(fit.error)[fit.error < 15]

hdata <- g.loess.scale[,top_genes]
# p2 <- Heatmap(t(hdata), cluster_columns = F)
p3 <- grid::grid.grabExpr(print(pheatmap(t(hdata), cluster_cols = F, fontsize_row = 6)))

tg_etrz <- mapIds(org.Mm.eg.db, top_genes, column = "ENTREZID", keytype = "SYMBOL")
epath <- enrichPathway(tg_etrz,
                       organism = "mouse",
                       universe = g_bckd,
                       readable = T)
write.csv(epath@result, "code/09-gary_objectives/outs/histology_DEA/GloboidSpikeNormal_pathway_enrichment.csv")
pdata <- epath@result %>% slice_head(n=30)
pdata$gr <- pdata$Count/length(top_genes)
pdata$Description_short <- sapply(pdata$Description, function(s) substring(s,1,35))
pdata$Description_short <- factor(pdata$Description_short, levels = pdata$Description_short[order(pdata$gr)])
p4 <- ggplot(pdata, aes(x=gr, y=Description_short)) + 
  geom_point(aes(size=Count, colour = -log(p.adjust))) +
  theme_bw() + theme(axis.text.y = element_text(size=8))
# p4 <- dotplot(epath, showCategory=30, includeAll=T) + theme(axis.text.y = element_text(size=8))

p_left <- plot_grid(plotlist = list(arrangeGrob(p1,p2, ncol=2),p3), nrow = 2, rel_heights = c(0.2,1))
plot_grid(p_left, p4, ncol = 2)

# [middle] normal, 15
x_ = seq(49)
y_ = dnorm(x_, mean=15, sd= 10, log=F)
y_ <- y_/max(y_)

fit.error <- abs(sweep(g.loess.scale, 1, y_, FUN = "-")) %>%
  colSums()
fit.error <- fit.error[order(fit.error)]
p1 <- ggplot(data.frame("dist"=x_, "y_"=y_), aes(dist,y_)) + geom_point() + theme_bw() + theme(axis.text.x = element_blank())
p2 <- ggplot(data.frame("error"=fit.error), aes(error)) + geom_density() + theme_bw() + geom_vline(xintercept = 15, col="red")
top_genes <- names(fit.error)[fit.error < 15]

hdata <- g.loess.scale[,top_genes]
# p2 <- Heatmap(t(hdata), cluster_columns = F)
p3 <- grid::grid.grabExpr(print(pheatmap(t(hdata), cluster_cols = F, fontsize_row = 6)))

tg_etrz <- mapIds(org.Mm.eg.db, top_genes, column = "ENTREZID", keytype = "SYMBOL")
epath <- enrichPathway(tg_etrz,
                       organism = "mouse",
                       universe = g_bckd,
                       readable = T)
write.csv(epath@result, "code/09-gary_objectives/outs/histology_DEA/peakLateNormal_pathway_enrichment.csv")
pdata <- epath@result %>% slice_head(n=30)
pdata$gr <- pdata$Count/length(top_genes)
pdata$Description_short <- sapply(pdata$Description, function(s) substring(s,1,35))
pdata$Description_short <- factor(pdata$Description_short, levels = pdata$Description_short[order(pdata$gr)])
p4 <- ggplot(pdata, aes(x=gr, y=Description_short)) + 
  geom_point(aes(size=Count, colour = -log(p.adjust))) +
  theme_bw() + theme(axis.text.y = element_text(size=8))
# p4 <- dotplot(epath, showCategory=30, includeAll=T) + theme(axis.text.y = element_text(size=8))

p_left <- plot_grid(plotlist = list(arrangeGrob(p1,p2, ncol=2),p3), nrow = 2, rel_heights = c(0.2,1))
plot_grid(p_left, p4, ncol = 2)

# [middle] normal, 35
x_ = seq(49)
y_ = dnorm(x_, mean=35, sd= 10, log=F)
y_ <- y_/max(y_)

fit.error <- abs(sweep(g.loess.scale, 1, y_, FUN = "-")) %>%
  colSums()
fit.error <- fit.error[order(fit.error)]
p1 <- ggplot(data.frame("dist"=x_, "y_"=y_), aes(dist,y_)) + geom_point() + theme_bw() + theme(axis.text.x = element_blank())
p2 <- ggplot(data.frame("error"=fit.error), aes(error)) + geom_density() + theme_bw() + geom_vline(xintercept = 15, col="red")
top_genes <- names(fit.error)[fit.error < 15]

hdata <- g.loess.scale[,top_genes]
# p2 <- Heatmap(t(hdata), cluster_columns = F)
p3 <- grid::grid.grabExpr(print(pheatmap(t(hdata), cluster_cols = F, fontsize_row = 6)))

tg_etrz <- mapIds(org.Mm.eg.db, top_genes, column = "ENTREZID", keytype = "SYMBOL")
epath <- enrichPathway(tg_etrz,
                       organism = "mouse",
                       universe = g_bckd,
                       readable = T)
write.csv(epath@result, "code/09-gary_objectives/outs/histology_DEA/peakEarlyNormal_pathway_enrichment.csv")
pdata <- epath@result %>% slice_head(n=30)
pdata$gr <- pdata$Count/length(top_genes)
pdata$Description_short <- sapply(pdata$Description, function(s) substring(s,1,35))
pdata$Description_short <- factor(pdata$Description_short, levels = pdata$Description_short[order(pdata$gr)])
p4 <- ggplot(pdata, aes(x=gr, y=Description_short)) + 
  geom_point(aes(size=Count, colour = -log(p.adjust))) +
  theme_bw() + theme(axis.text.y = element_text(size=8))
# p4 <- dotplot(epath, showCategory=30, includeAll=T) + theme(axis.text.y = element_text(size=8))

p_left <- plot_grid(plotlist = list(arrangeGrob(p1,p2, ncol=2),p3), nrow = 2, rel_heights = c(0.2,1))
plot_grid(p_left, p4, ncol = 2)

dev.off()


##########################################
'''
5.	Could one generate a “spectrum” of macrophage reactivity states spreading out from globoid cells??
'''

seu.twitcher$PAS_dist <- PAS_dmin

Idents(seu.twitcher) <- "predicted_ct"
seu.twitcher <- FindSubCluster(seu.twitcher,
                               cluster="Perivascular macrophages",
                               subcluster.name = "macrophage_cl",
                               graph.name = "Spatial_snn",
                               resolution = 0.3)
pdata <- seu.twitcher@meta.data %>%
  filter(predicted_ct == "Perivascular macrophages")
cl_count <- table(pdata$macrophage_cl)
pdata <- pdata %>% filter(macrophage_cl %in% names(cl_count)[cl_count > 50])

p1 <- ggplot(pdata, aes(x=PAS_dist, y=macrophage_cl))+
  geom_density_ridges()+
  theme_minimal()

m_cells <- Cells(seu.twitcher)[seu.twitcher$predicted_ct == "Perivascular macrophages"]
p2 <- DimPlot(seu.twitcher, cells=m_cells, group.by = "macrophage_cl") + theme(legend.position = "none")
p3 <- FeaturePlot(seu.twitcher, cells=m_cells, features = "PAS_dist") + theme(legend.position = "none")

plot_grid(p1,p2,p3, nrow=3)


