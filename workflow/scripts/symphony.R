#---Load required libraries and functions ####
library(symphony)
library(singlecellmethods)

source('resources/software/symphony/libs.R') # imports
source('resources/software/symphony/utils.R') # color definitions and plotting functions
source('resources/software/symphony/utils_seurat.R') # color definitions and plotting functions

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

# Filter datasets for common genes
commonGenes <- intersect(rownames(transfExp), rownames(queryExp))
transfExp <- transfExp[commonGenes,]
queryExp <- queryExp[commonGenes,]


# #---Option 1: Build Reference from Harmony ####
# var_genes = vargenes_vst(transfExp, groups = as.character(transfMeta$Cell_type), topn = 2000)
# ref_exp = transfExp[var_genes, ]
# dim(ref_exp)
#
# # Calculate and save the mean and standard deviations for each gene
# vargenes_means_sds = tibble(symbol = var_genes, mean = Matrix::rowMeans(ref_exp))
# vargenes_means_sds$stddev = singlecellmethods::rowSDs(ref_exp, vargenes_means_sds$mean)
#
# # Scale data using calculated gene means and standard deviations
# ref_exp_scaled = singlecellmethods::scaleDataWithStats(ref_exp, vargenes_means_sds$mean, vargenes_means_sds$stddev, 1)
#
# # Run SVD, save gene loadings (s$u)
# set.seed(0)
# s = irlba(ref_exp_scaled, nv = 20)
# Z_pca_ref = diag(s$d) %*% t(s$v) # [pcs by cells]
# loadings = s$u
#
# # Run Harmony integration
# set.seed(0)
# ref_harmObj = harmony::HarmonyMatrix(
#   data_mat = t(Z_pca_ref),  ## PCA embedding matrix of cells
#   meta_data = transfMeta, ## dataframe with cell labels
#   theta = c(2),             ## cluster diversity enforcement
#   vars_use = c('Cell_type'),    ## variable to integrate out
#   nclust = 100,             ## number of clusters in Harmony model
#   max.iter.harmony = 20,
#   return_object = TRUE,     ## return the full Harmony model object
#   do_pca = FALSE            ## don't recompute PCs
# )
#
# # Compress a Harmony object into a Symphony reference
# reference = symphony::buildReferenceFromHarmonyObj(
#   ref_harmObj,            # output object from HarmonyMatrix()
#   transfMeta,           # reference cell metadata
#   vargenes_means_sds,     # gene names, means, and std devs for scaling
#   loadings,               # genes x PCs matrix
#   verbose = TRUE,         # verbose output
#   do_umap = TRUE,         # Set to TRUE only when UMAP model was saved for reference
#   save_uwot_path = './testing_uwot_model_1')

#---Option 2: Build from scratch ####
print("Building model")
reference = symphony::buildReference(
  transfExp,
  transfMeta,
  vars = c('Cell_type'),     # variables to integrate over
  K = 100,                   # number of Harmony clusters
  verbose = TRUE,            # verbose output
  do_umap = TRUE,            # can set to FALSE if want to run umap separately later
  do_normalize = TRUE,      # set to TRUE if input counts are not normalized yet
  vargenes_method = 'vst',   # method for variable gene selection ('vst' or 'mvp')
  vargenes_groups = 'Cell_type', # metadata column specifying groups for variable gene selection
  topn = 2000,               # number of variable genes to choose per group
  d = 20,                    # number of PCs
  save_uwot_path = tempfile("uwot") #snakemake@output[["reference_uwottemp"]] # Must declare, otherwise next steps wont work
)

umap_labels = cbind(transfMeta, reference$umap$embedding)

fig.size(3, 5)

png(filename = snakemake@output[["reference_umap"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])
plotBasic(umap_labels,
          title = 'Reference',
          color.by = 'Cell_type')
dev.off()

#---Map query ####
print("Mapping query cells")
query = mapQuery(queryExp,              # query gene expression (genes x cells)
                 queryMeta,             # query metadata (cells x attributes)
                 reference,             # Symphony reference object
                 vars = NULL,           # Query batch variables to harmonize over (NULL treats query as one batch)
                 do_normalize = TRUE,   # perform log(CP10k) normalization on query (set to FALSE if already normalized)
                 do_umap = TRUE)        # project query cells into reference UMAP

# Predict query cell types using k-NN
query = knnPredict(query, reference, reference$meta_data$Cell_type, k = 5)

# Sync the column names for both data frames
reference$meta_data$cell_type_pred_knn = NA
reference$meta_data$cell_type_pred_knn_prob = NA
reference$meta_data$ref_query = 'reference'
query$meta_data$ref_query = 'query'

# Add the UMAP coordinates to the metadata
meta_data_combined = rbind(query$meta_data, reference$meta_data)
umap_combined = rbind(query$umap, reference$umap$embedding)
umap_combined_labels = cbind(meta_data_combined, umap_combined)

# Create query only umap dataframe
umap_query = umap_combined_labels[umap_combined_labels$ref_query=="query",]

# Plot UMAP visualizations
print("Plotting umaps")
fig.size(3,5)
png(filename = snakemake@output[["umap_query"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])
plotBasic(umap_query,
          title = 'Query cells',
          color.by = 'cell_type_pred_knn')
dev.off()

fig.size(3, 5)
png(filename = snakemake@output[["umap_mixed_queryAndReference"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])
plotBasic(umap_combined_labels,
          title = 'Reference and query cells',
          color.by = 'ref_query')
dev.off()

# fig.size(3, 7)
# plotBasic(umap_combined_labels,
#           title = 'Reference and query cells',
#           color.by = 'Cell_type',
#           facet.by = 'ref_query')
# dev.off()
