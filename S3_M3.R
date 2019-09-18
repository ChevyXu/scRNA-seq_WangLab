```r
library(monocle3)
library(tidyverse)
library(imputeTS)
library(Seurat)
library(RColorBrewer)
library(Matrix)
library(cowplot)
library(future)
require(scales)
temp <- list.files(path = "../Raw_data/",pattern="*_gene_exon_tagged.dge.txt.gz")
temp
temp <- temp[c(7,8)]

# ===========================================================================================data readin
filter_umi <- 450
name <- character()
for(i in 1:length(temp)){
  name[i] <- unlist(strsplit(temp[i],"_out_gene_exon_tagged.dge.txt.gz"))[1]
  message(paste(name[i], "is loading"))
  
  tmpvalue<-read.table(paste0("../Raw_data/", temp[i]), sep = "\t", quote = "", row.names = 1, header = T)
  message(paste(name[i], "is loaded, now is adding name"))
  
  colnames(tmpvalue) <- paste0(name[i], "-", colnames(tmpvalue))         
  message(paste0(name[i], "'s name added, now filtering ", filter_umi))
  
  tmpvalue <- tmpvalue[,colSums(tmpvalue) >= filter_umi]
  message(paste(name[i], "cells above", filter_umi, "filtered"))
  
  assign(name[i], tmpvalue)
  rm(tmpvalue)
}
message("data loading done, and strat merge counts file")

# ===========================================================================================UMI summary
summary <- NULL
for (s in 1:length(temp)) {
  summary <- cbind(summary, as.matrix(summary(colSums(get(name[s])))))
}
colnames(summary) <- paste0(name, "'s UMI")
summary

# ===========================================================================================data merge
dge <- get(name[1])
for(p in 2:length(temp)) {
  dge_temp <- get(name[p])
  dge_merge <- dplyr::full_join(dge %>% mutate(name = rownames(dge)), 
                                dge_temp %>% mutate(name = rownames(dge_temp)))
  rownames(dge_merge) <- dge_merge$name
  dge <-  dplyr::select(dge_merge, -(name)) %>% na.replace(0)
  rm(dge_temp)
  rm(dge_merge)
}
message("data merge done")
saveRDS(dge, file = "4th_dge_raw_data.rds")
dge <- readRDS("4th_dge_raw_data.rds")

# ===========================================================================================pbmc
filter_gene = 400
filter_cell = 5
ngene = 2000
permito = 0.1

pbmc <- CreateSeuratObject(counts = dge, 
                           project = "AML1-ETO-MigR1", 
                           min.features = filter_gene, 
                           min.cells = filter_cell, 
                           meta.data = metadata)

levels <- levels(factor(pbmc@meta.data$tech))
prefix <- character()
for (n in 1:length(levels)){prefix <- paste0(levels[n], "_", prefix)}
table(pbmc@meta.data$tech)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
pbmc <- subset(pbmc, subset = nFeature_RNA > 350 & nFeature_RNA < 2000 & percent.mt < 8)
hist(pbmc$percent.mt)
# # NormalizeData
# pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# # FindVariableFeatures
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1000)
# 
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(pbmc), 10)
# VariableFeaturePlot(pbmc)
# LabelPoints(plot = plot1, points = top10, repel = TRUE)
# 
# # SCTransform
# pbmc <- ScaleData(pbmc, vars.to.regress = c("percent.mt", "tech"))

pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
pbmc <- SCTransform(object = pbmc, variable.features.n = 1500, vars.to.regress = c("percent.mt", "tech"))

# ----------------------------------------------------------------------------------------------------RUN PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
ElbowPlot(pbmc)
dims = 11
pbmc <- FindNeighbors(pbmc, dims = 1:dims)
pbmc <- FindClusters(pbmc, resolution = 0.8)

# ----------------------------------------------------------------------------------------------------RUN umap
pbmc <- RunUMAP(pbmc, dims = 1:dims)
DimPlot(pbmc, reduction = "umap")
CombinePlots(plots = list(
  DimPlot(object = pbmc, pt.size = 1, reduction = "umap", group.by = "tech"),
  DimPlot(object = pbmc, pt.size = 1, reduction = "umap", group.by = "SCT_snn_res.0.5")))

# ----------------------------------------------------------------------------------------------------RUN FlITSNE
pbmc <- RunTSNE(pbmc,
                reduction = "pca", 
                cells = NULL,
                dims = 1:dims, 
                features = NULL, 
                seed.use = 1, 
                tsne.method = "Rtsne",
                reduction.name = "tsne", 
                reduction.key = "tsne_")
CombinePlots(plots = list(
  DimPlot(object = pbmc, pt.size = 1, reduction = "tsne", group.by = "tech"),
  DimPlot(object = pbmc, pt.size = 1, reduction = "tsne", group.by = "SCT_snn_res.0.8", label = T)))

DimPlot(object = pbmc, pt.size = 1, reduction = "tsne", 
        group.by = "SCT_snn_res.0.5", split.by = "tech", label = T, repel = T, label.size = 5)

source("C:/Users/DELL/Downloads/FIt-SNE-master/examples/test.R")
pbmc <- RunTSNE(object = pbmc, 
                dims = 1:dims, 
                reduction.name = "FItSNE", 
                reduction.key = "FItSNE_",
                tsne.method = "FIt-SNE", 
                fast_tsne_path = "C:/Users/DELL/Downloads/FIt-SNE-master/bin/FItSNE.exe")
CombinePlots(plots = list(
  DimPlot(object = pbmc, pt.size = 1, reduction = "FItSNE", group.by = "tech"),
  DimPlot(object = pbmc, pt.size = 1, reduction = "FItSNE", group.by = "SCT_snn_res.0.5")))

DimPlot(object = pbmc, pt.size = 1, reduction = "FItSNE", 
        group.by = "SCT_snn_res.0.5", split.by = "tech", label = T, repel = T, label.size = 8)
# ----------------------------------------------------------------------------------------------------stat
stat <- data.frame(cellident= pbmc@meta.data$tech, clusterident=pbmc@meta.data$SCT_snn_res.0.5) %>%
  group_by(clusterident, cellident) %>%
  summarise(count=n()) %>%
  group_by(clusterident)

stat %>%  
  summarise(n()) %>%
  dplyr::filter(`n()` == 1)

table <- data.frame(table(pbmc@meta.data$tech))
info <- data.frame(clusterident =as.factor(0:(length(stat$clusterident)/length(levels)-1)))

for (i in 1:length(table$Var1)) {
  percentage <- left_join(stat %>% dplyr::filter(cellident == table$Var1[i]) %>% summarise(count / table$Freq[i] * 100),
                          stat %>% dplyr::filter(cellident == table$Var1[i]) %>% summarise(count), 
                          by = "clusterident")
  colnames(percentage) <- c("clusterident", paste0(table$Var1[i], "_ratio"), paste0(table$Var1[i], "_counts"))
  info <- left_join(info, percentage, by = "clusterident")
}

result <- info %>% mutate(ratio = info[,2]/info[,4])
result

summary(pbmc@meta.data$nFeature_RNA)
summary(pbmc@meta.data$nCount_RNA)

cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, ident.2 = 1, only.pos = T)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers.0 <- FindMarkers(object = pbmc, ident.1 = c(0,1), ident.2 = c(2:11), logfc.threshold = 0.1, only.pos = T)

for (c in 1:nrow(result)) {
  CombinePlots(plots = list(
    FeaturePlot(object = pbmc, features = pbmc.markers[pbmc.markers$cluster == c,]$gene[1], 
                reduction = "FItSNE", cols = c("grey", "red"), pt.size = 1, sort.cell = T),
    FeaturePlot(object = pbmc, features = pbmc.markers[pbmc.markers$cluster == c,]$gene[2], 
                reduction = "FItSNE", cols = c("grey", "red"), pt.size = 1, sort.cell = T),
    FeaturePlot(object = pbmc, features = pbmc.markers[pbmc.markers$cluster == c,]$gene[3], 
                reduction = "FItSNE", cols = c("grey", "red"), pt.size = 1, sort.cell = T),
    FeaturePlot(object = pbmc, features = pbmc.markers[pbmc.markers$cluster == c,]$gene[4], 
                reduction = "FItSNE", cols = c("grey", "red"), pt.size = 1, sort.cell = T)))
}


DimPlot(object = pbmc, pt.size = 2, reduction = "FItSNE", label = T, label.size = 10, cols = "Paired", 
        order = c(0,1,2,3,4,8,6,7,5,9)) + NoLegend()
ggsave("overviewplot.png", width = 12, height = 10)

DimPlot(object = pbmc, pt.size = 2, reduction = "FItSNE", label = T, label.size = 10, split.by = "tech", cols = "Paired",
        order = c(0,1,2,3,4,8,6,7,5,9)) + NoLegend()
ggsave("overview_splitplot.png", width = 24, height = 10)

pbmc <- RenameIdents(pbmc, 
                     `0` = "CMP-1 (myeloid)", 
                     `1` = "CMP-2 (myeloid)", 
                     `2` = "GMP (myeloid)", 
                     `3` = "Erythroid (myeloid)", 
                     `4` = "Promyelocyte-1 (myeloid)", 
                     `5` = "Progenitor cell", 
                     `6` = "MEP (myeloid)", 
                     `7` = "Promyelocyte-2 (myeloid)", 
                     `8` = "HSC", 
                     `9` = "Monocyte (myeloid)", 
                     `10` = "pDC (lymphoid)")

DimPlot(pbmc, reduction = "FItSNE", label = TRUE, pt.size = 2, label.size = 10, cols = "Paired", 
        order = c("CMP-1 (myeloid)", "CMP-2 (myeloid)", "GMP (myeloid)", "Erythroid (myeloid)","Promyelocyte-1 (myeloid)",
                  "HSC", "MEP (myeloid)", "Promyelocyte-2 (myeloid)", "Progenitor cell","Monocyte (myeloid)", "pDC (lymphoid)")) + NoLegend()
ggsave("annotation.png", width = 12, height = 10)


DimPlot(pbmc, reduction = "FItSNE", label = TRUE, pt.size = 2, label.size = 10, split.by = "tech", cols = "Paired", 
        order = c("CMP-1 (myeloid)", "CMP-2 (myeloid)", "GMP (myeloid)", "Erythroid (myeloid)","Promyelocyte-1 (myeloid)",
                  "HSC", "MEP (myeloid)", "Promyelocyte-2 (myeloid)", "Progenitor cell","Monocyte (myeloid)", "pDC (lymphoid)")) + NoLegend()
ggsave("annotation_split.png", width = 24, height = 10)

markers.to.plot <- c("TPSAB1", "TPSB2", "LGALS1", "SH3BGRL3", "TUBA1B", "H2AFZ", "HMGB2","MALAT1", "SOX4", "PLIN2", 
                     "CTSG", "CMA1", "CPA3", "TM4SF1", "CRHBP", "GNG11", "SOD2", "CHI3L1" ,"CD14", "C1QB", "LGMN")
DotPlot(pbmc, features = markers.to.plot, cols = "Spectral", dot.scale = 10) + RotatedAxis()
ggsave("dotplot.png", width = 24, height = 10)

DoHeatmap(subset(pbmc, downsample = 100), features = markers.to.plot, size = 3)

DimPlot(object = pbmc, reduction = "FItSNE", group.by = "tech", 
        cells = cells.ae, cols = c("red"), do.label = F, no.legend = F)
DimPlot(object = pbmc, reduction = "FItSNE", group.by = "tech", 
        cells = cells.migr1, cols = c("blue"), do.label = F, no.legend = F)


load("../Raw_data/pbmc_4th.Rdata")
# ==========================================================================================cds construct
# 1500 feature data: pbmc@assays$SCT@scale.data[1:10, 1:5]
hist(colSums(pbmc@assays$RNA@counts), breaks = 100)
hist(colSums(pbmc@assays$SCT@counts), breaks = 100)
hist(colSums(pbmc@assays$SCT@data), breaks = 100)
hist(colSums(pbmc@assays$SCT@scale.data), breaks = 100)

dge_new <- as.matrix(pbmc@assays$SCT@data)

cell_metadata <- matrix(unlist(strsplit(colnames(dge_new), split = "-")), 
                        ncol = 2, byrow = T, dimnames = list(colnames(dge_new), c("Treat", "UMI")))

gene_metadata <- data.frame(gene_short_name = rownames(dge_new), row.names = rownames(dge_new))

cds <- new_cell_data_set(expression_data = dge_new,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

# ==========================================================================================get group UMI
metadata <- dplyr::left_join(data.frame(pbmc@meta.data) %>% dplyr::mutate(full_UMI = rownames(data.frame(pbmc@meta.data))),
                             data.frame(pbmc@active.ident) %>% dplyr::mutate(full_UMI = rownames(data.frame(pbmc@active.ident))),
                             by = "full_UMI")
head(metadata)

cds_metadata <- data.frame(colData(cds))
head(cds_metadata)

metadata_join <- dplyr::left_join(cds_metadata %>% dplyr::mutate(full_UMI = rownames(cds_metadata)),
                                  metadata, by = "full_UMI") %>%
  dplyr::select(-c("UMI.x", "SCT_snn_res.1.2", "nCount_SCT", "nFeature_SCT"))
head(metadata_join)

colData(cds)$seurat <- metadata_join$seurat_clusters
colData(cds)$celltype <- metadata_join$pbmc.active.ident
head(colData(cds))

# ==========================================================================================PCA and norm
cds@reducedDims$PCA <-              pbmc@reductions$pca@cell.embeddings[,1:15]
cds@preprocess_aux$gene_loadings <- pbmc@reductions$pca@feature.loadings[,1:15]
cds@preprocess_aux$prop_var_expl <- pbmc@reductions$pca@stdev[1:15]

# cds <- preprocess_cds(cds,
#                       method = "PCA",
#                       use_genes = c(pbmc@assays$SCT@var.features),
#                       num_dim = 12,
#                       norm_method = "none",
#                       scaling = T)
plot_pc_variance_explained(cds)

# ==========================================================================================tSNE
# cds <- reduce_dimension(cds, reduction_method = "tSNE")
# plot_cells(cds,
#            reduction_method = "tSNE",
#            color_cells_by = "Treat", label_cell_groups = FALSE,
#            cell_size = 1)
# plot_cells(cds,
#            reduction_method = "tSNE",
#            genes = c("TM4SF1", "CRHBP", "GNG11", "RUNX1T1"),
#            cell_size = 1)
# 
# # ==========================================================================================cluster
# cds = cluster_cells(cds,
#                     reduction_method = "tSNE",
#                     k = 150)
# plot_cells(cds, reduction_method = "tSNE", cell_size = 1)
# plot_cells(cds, reduction_method = "tSNE", group_cells_by = "cluster")

# ==========================================================================================UMAP
cds@reducedDims$UMAP <- pbmc@reductions$umap@cell.embeddings

cds <- reduce_dimension(cds, 
                        cores = 1,
                        reduction_method = "UMAP", 
                        preprocess_method = "PCA")
plot_cells(cds, 
           reduction_method = "UMAP", 
           color_cells_by = "Treat", 
           label_cell_groups = FALSE,
           cell_size = 1)
plot_cells(cds,
           reduction_method = "UMAP",
           genes = c("TM4SF1"),
           cell_size = 1)

# ==========================================================================================cluster
cds = cluster_cells(cds, 
                    reduction_method = "UMAP", 
                    k = 200)
plot_cells(cds, reduction_method = "UMAP", cell_size = 1)

# ==========================================================================================Marker plot
# marker_test_res = top_markers(cds, reduction_method = "UMAP", group_cells_by = "cluster")
# top_specific_markers = marker_test_res %>%
#   filter(fraction_expressing >= 0.10) %>%
#   group_by(cell_group) %>%
#   top_n(1, pseudo_R2)
# 
# top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="cluster",
#                     ordering_type="maximal_on_diag",
#                     max.size=3)

# ==========================================================================================trajectory
values <- c(brewer.pal(9, "Set1"), brewer.pal(9, "Pastel1"), "black", "green")

cds <- learn_graph(cds)
plot_cells(cds, 
           color_cells_by = "celltype", 
           cell_size = 1.5, 
           label_leaves = F,
           graph_label_size = 0, 
           label_cell_groups = F, 
           trajectory_graph_segment_size = 2,
           trajectory_graph_color = "red") +
  scale_color_manual(values = values[c(8,8,7,6,5,4,3,11,1,20,19)])

# cds <-  order_cells(cds, root_cells = as.character(colData(cds)$fUll_umi[colData(cds)$celltype == "HSC"] %>% na.omit()))
cds <-  order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = F,
           label_leaves = T,
           label_branch_points = T,
           graph_label_size = 4,
           cell_size = 2) 
scale_color_gradient2(low = "#E41A1C", high = "#377EB8", mid = "#999999", midpoint = 10)
brewer.pal(9, "Set1")

cds@clusters$UMAP$clusters
save.image("Seurat3_Monocle3_4th.rdata")
# load("monocle3_4th.rdata")
```
