setwd("/mnt/isilon/cscb/Projects/10x_Visium/zhangj/SCTC-YZ-149")

library(Seurat)
library(SeuratDisk)
library(loupeR)
library(dplyr)
library(ggplot2)

seu <- readRDS("processed/seurat/integrated_dat_cellpred.RDS")
new_meta <- read.csv("Code/05-extra_analysis/new_meta.csv", row.names = 1)

seu$center_dist <- new_meta$center_dist

###
# making some plots
#   merge_if_norm > 0.1 

pdf(
  "outs/Seurat/spatial_overlay_merge.pdf",
  height = 8,
  width = 11
)
snames <- unique(seu$orig.ident)
for(i in 1:length(snames)){
  if(i == 1){
    print(SpatialFeaturePlot(seu, features = c("merge_IF","merge_IF_norm"),
                             ncol = 2,alpha = 0.5,
                             images = paste("slice",snames[i], sep = "_")) +
            ggtitle(snames[i]))
  }
  else{
    print(SpatialFeaturePlot(seu, features = c("merge_IF","merge_IF_norm"),
                             ncol = 2,alpha = 0.5,
                             images = paste0("slice_",snames[i], ".", i)) +
            ggtitle(snames[i]))
  }
  
}

dev.off()

###
# Examining distributions
plot(density(seu$merge_IF_norm))
table(seu$var1[seu$merge_IF_norm>0.1])

pdata <- seu@meta.data
pdata <- pdata[pdata$merge_IF_norm > 0.1,]
pdata$var1 <- factor(pdata$var1, levels = c("Control", "aCD40 24hr", "aCD40 1wk", "ICB 1wk", "aCD40/ICB 1wk"))
ggplot(pdata, aes(x=var1, y=center_dist, color=var1))+
  geom_boxplot(outlier.size = 0) +
  geom_jitter(width = 0.1, size=pdata$merge_IF_norm)+
  theme_bw()

pdata$quant <- cut(pdata$merge_IF_norm, 
             breaks = quantile(pdata$merge_IF_norm, seq(0,1,0.25)),
             labels = paste0("Q", seq(4)),
             include.lowest = T)

ggplot(pdata, aes(x=var1, y=center_dist, color=var1))+
  geom_boxplot(outlier.size = 0) +
  geom_jitter(width = 0.1, size=pdata$merge_IF_norm)+
  theme_bw()+
  facet_wrap(~quant, scales = "free_y")

ggplot(pdata, aes(x=var1, y=center_dist, color=var1))+
  geom_violin() +
  geom_jitter(width = 0.1, size=pdata$merge_IF_norm)+
  theme_bw()+
  facet_wrap(~quant, scales = "free_y")

# what if we scale to maximum center distance
sum_stats <- seu@meta.data %>% group_by(sample.id) %>% summarise("max"=max(center_dist),
                                                                 "min"=min(center_dist),
                                                                 "n" = n(),
                                                                 "quant_max"=quantile(center_dist, 0.95))

# ...I'm wondering if some of the double tissue labels are wrong...
scaler <- sum_stats$quant_max; names(scaler) <- sum_stats$sample.id
pdata$scaled_dist <- pdata$center_dist/scaler[pdata$sample.id]

ggplot(pdata, aes(x=var1, y=scaled_dist, color=var1))+
  geom_boxplot(outlier.size = 0) +
  geom_jitter(width = 0.1, size=pdata$merge_IF_norm)+
  theme_bw()

ggplot(pdata, aes(x=var1, y=scaled_dist, color=var1))+
  geom_violin() +
  geom_jitter(width = 0.1, size=pdata$merge_IF_norm)+
  theme_bw()+
  facet_wrap(~quant)

###
# Exporting
seu$scaled_dist <- seu$center_dist/scaler[seu$sample.id]

saveRDS(seu, "processed/seurat/integrated_dat_cellpred.RDS")

# metadata
out <- seu@meta.data
write.csv(out, "processed/seurat/metadata.csv")

# for cellxgene
out <- seu
out[["Spatial"]] <- as(object = out[["Spatial"]], Class = "Assay")
out[["SCT"]] <- CreateSCTAssayObject(data = out@assays$SCT$data)
DefaultAssay(out) <- "SCT"
SaveH5Seurat(out, filename = "processed/seurat/integrated_dat_cellPred2.h5Seurat", overwrite = T)
Convert("processed/seurat/integrated_dat_cellPred2.h5Seurat", dest = "h5ad", overwrite = T)

# for Loupe
DefaultAssay(seu) <- "Spatial"
out <- DietSeurat(seu, assays = c("Spatial"), dimreducs = c("umap"))
# out[["SCT"]] <- CreateSCTAssayObject(data = out@assays$SCT$data)
# out$seurat_clusters <- as.character(out$seurat_clusters)
# DefaultAssay(out) <- "SCT"

create_loupe_from_seurat(
  out,
  output_dir = "outs/Loupe",
  output_name = "integrated_all_slides_v3_raw",
  force=T
)

# metadata for the spatial cloupes

