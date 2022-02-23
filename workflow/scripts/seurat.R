#---Load required libraries ####
library(Seurat)
library(magrittr)
library(ggplot2)

#---Load data ####
# Load transfer data
print("Loading data")
tmp <- load(snakemake@input[["metaTrain"]]) #sample table
tmp2 <- load(snakemake@input[["expTrain"]]) #expression matrix
gse_meta <- as.data.frame(gse_meta)

transfExp <- gse_exp
transfMeta <- gse_meta

# Load Query Data
tmp <- load(snakemake@input[["metaQuery"]]) #sample table
tmp2 <- load(snakemake@input[["expQuery"]]) #expression matrix
gse_meta <- as.data.frame(gse_meta)

queryExp <- gse_exp
queryMeta <- gse_meta

rm(gse_meta)
rm(gse_exp)
rm(tmp)
rm(tmp2)

# Rename important colnames to that given by the user in config.yaml
colnames(transfMeta)[colnames(transfMeta) == snakemake@params[["Cell_type_colname"]]] <- "Cell_type"
colnames(transfMeta)[colnames(transfMeta) == snakemake@params[["Cell_ID_colname"]]] <- "Cell_ID"
colnames(queryMeta)[colnames(queryMeta) == snakemake@params[["Cell_description_colname"]]] <- "Cell_type"
colnames(queryMeta)[colnames(queryMeta) == snakemake@params[["Cell_ID_colname"]]] <- "Cell_ID"

#---Create Seurat object ####
print("Creating Seurat object")
obj <- CreateSeuratObject(counts = transfExp, meta.data = transfMeta) %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)

p1 <- DimPlot(obj, reduction = "umap", group.by = "Cell_type") + ggtitle("Transfer Data UMAP")

print("Creating transfer UMAP")
png(filename = snakemake@output[["transfdata_umap"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])
p1
dev.off()

#---Transfer Data ####
print("Predicting query data")
query <- CreateSeuratObject(counts = queryExp, meta.data = queryMeta)

query.anchors <- FindTransferAnchors(reference = obj,
                                     query = query,
                                     dims = 1:30,
                                     reference.reduction = "pca")

predictions <- TransferData(anchorset = query.anchors,
                            refdata = obj$Cell_type,
                            dims = 1:30)

query <- AddMetaData(query, metadata = predictions)

# Evaluation
# features = c("PBX1", "TH", "DDC")
features = snakemake@params[["Seurat_features"]]
print(features)
text <- "No features selected"

print("Creating feature plots")

png(filename = snakemake@output[["query_violin"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])
if (!is.null(features)) {
  VlnPlot(query, features = features, group.by = "predicted.id")
} else {
  ggplot() + 
    annotate("text", x = 4, y = 25, size=8, label = text) + 
    theme_void()
}
dev.off()

query_scale <- ScaleData(query)

png(filename = snakemake@output[["query_heatmap"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])
if (!is.null(features)) {
  DoHeatmap(query_scale, features = features, group.by = "predicted.id")
} else {
  ggplot() + 
    annotate("text", x = 4, y = 25, size=8, label = text) + 
    theme_void()
}
dev.off()

# UMAP Projection
print("Creating query UMAP")
obj <- RunUMAP(obj, dims = 1:30,
               reduction = "pca",
               return.model = TRUE)

query <- MapQuery(anchorset = query.anchors,
                  reference = obj,
                  query = query,
                  refdata = list(celltype = "Cell_type"),
                  reference.reduction = "pca",
                  reduction.model = "umap")

p2 <- DimPlot(query,
              reduction = "ref.umap",
              group.by = "predicted.celltype",
              label = TRUE,
              label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")

png(filename = snakemake@output[["refAndQuery_umap"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])
p2
dev.off()
