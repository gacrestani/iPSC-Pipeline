#---Load required libraries and functions ####
library(symphony)
library(singlecellmethods)
library(fusca)
library(magrittr)

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


#---Build Symphony model ####
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

#---Building query Cellrouter object ####
print("Creating CellRouter object")
cellrouter <- CreateCellRouter(rawdata=queryExp,
                               min.cells=3,
                               min.genes=0,
                               is.expr=0) %>%
  fusca::Normalize() %>%
  fusca::scaleData() %>%
  fusca::computePCA(num.pcs=50, seed=42) %>%
  fusca::computeTSNE(num.pcs=11, seed=42, max_iter=1000) %>%
  fusca::computeUMAP(num.pcs=11)
  
cellrouter <- customSpace(cellrouter, cellrouter@umap$cell.embeddings)

#---Query mapping ####
query_exp <- cellrouter@assays$RNA@ndata
query_metadata <- cellrouter@assays$RNA@sampTab
query = mapQuery(query_exp,             # query gene expression (genes x cells)
                 query_metadata,        # query metadata (cells x attributes)
                 reference,             # Symphony reference object
                 do_normalize = FALSE,  # perform log(CP10k+1) normalization on query
                 do_umap = TRUE)        # project query cells into reference UMAP

query = knnPredict(query,
                   reference,
                   reference$meta_data$Cell_type,
                   k = 25)

# Sort out colnames
query$meta_data$ref_query = 'Query'
reference$meta_data$ref_query = 'Reference'
colnames(query$meta_data)[which(names(query$meta_data) == "Cell_type")] <- "Cell_type_OLD"
colnames(query$meta_data)[which(names(query$meta_data) == "cell_type_pred_knn")] <- "Cell_type"
colnames(query$meta_data)[which(names(query$meta_data) == "sample_id")] <- "Cell_ID"

# Filter metadata to contain only important columns
keep <- c("Cell_ID", "Cell_type", "ref_query")
meta_data_combined = rbind(query$meta_data[,keep], reference$meta_data[,keep])
umap_combined = rbind(query$umap, reference$umap$embedding)

umap_combined_labels = cbind(meta_data_combined, umap_combined)

png(filename = snakemake@output[["umap_queryAndReference"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])
plotBasic(umap_combined_labels,
          title = 'Reference and Query cells',
          color.by = 'Cell_type',
          facet.by = 'ref_query')
dev.off()