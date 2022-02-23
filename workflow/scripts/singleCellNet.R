#---Load required libraries and functions####
library(singleCellNet)
library(pheatmap)

source('resources/software/singleCellNet/score_function.R')

#---Load data ####
# Load Query Data
print("Loading data")
queryMeta <- utils_loadObject(snakemake@input[["metaQuery"]]) #sample table
queryExp <- utils_loadObject(snakemake@input[["expQuery"]]) #expression matrix

# Load Training Data
transfMeta <- utils_loadObject(snakemake@input[["metaTrain"]]) #sample table
transfExp <- utils_loadObject(snakemake@input[["expTrain"]]) #expression matrix

# Remove "_" from rownames to prevent name errors.
rownames(queryExp) = gsub("_",".",rownames(queryExp))
rownames(transfExp) = gsub("_",".",rownames(transfExp))

# Rename important colnames to that given by the user in config.yaml
colnames(transfMeta)[colnames(transfMeta) == snakemake@params[["Cell_type_colname"]]] <- "Cell_type"
colnames(transfMeta)[colnames(transfMeta) == snakemake@params[["Cell_ID_colname"]]] <- "Cell_ID"
colnames(queryMeta)[colnames(queryMeta) == snakemake@params[["Cell_description_colname"]]] <- "Cell_type"
colnames(queryMeta)[colnames(queryMeta) == snakemake@params[["Cell_ID_colname"]]] <- "Cell_ID"

#---Training the classifier ####
# Get common genes and filter expTables to contain only common genes
commonGenes <- intersect(rownames(transfExp), rownames(queryExp))
transfExp <- transfExp[commonGenes,]
queryExp <- queryExp[commonGenes,]

# Split data for training and assessment, and transform training data
stList <- splitCommon(sampTab = transfMeta,
                      ncells = 100,
                      dLevel = "Cell_type")

stTrain <- stList[[1]]
expTrain <- transfExp[ ,rownames(stTrain)]

stTest <- stList[[2]]
expTest <- transfExp[ ,rownames(stTest)]

# Train the classifier
print("Creating classifier")
system.time(
  class_info <- scn_train(stTrain = stTrain,
                          expTrain = expTrain,
                          nTopGenes = 10,
                          nRand = 70,
                          nTrees = 1000,
                          nTopGenePairs = 25,
                          dLevel = "Cell_type",
                          colName_samp = "Cell_ID"))

#---Assessing the classifier ####
print("Assessing classifier")

# Train data UMAP
# Predict test data
classRes_val_all <- scn_predict(cnProc = class_info[['cnProc']],
                                expDat = expTest,
                                nrand = 50)

# Summarize data for precision recall curve and metrics
tm_heldoutassessment <- assess_comm(ct_scores = classRes_val_all,
                                    stTrain = stTrain,
                                    stQuery = stTest,
                                    dLevelSID = "Cell_ID",
                                    classTrain = "Cell_type",
                                    classQuery = "Cell_type",
                                    nRand = 50)

# Plot classifier assessment and assessment metrics
print("Creating assessment curves")
png(filename = snakemake@output[["assessment_curves"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])
plot_PRs(tm_heldoutassessment)
dev.off()

print("Creating assessment metrics")
png(filename = snakemake@output[["assessment_metrics"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])
plot_metrics(tm_heldoutassessment)
dev.off()

# Heatmap for the test data
# Create a name vector label used later in classification heatmap where the values are cell types/ clusters and names are the sample names
print("Creating training data heatmap")
nrand <- 50
sla <- as.vector(stTest$Cell_type)
names(sla) <- as.vector(stTest$Cell_ID)
slaRand <- rep("rand", nrand)
names(slaRand) <- paste("rand_", 1:nrand, sep='')
sla <- append(sla, slaRand) #include in the 'random cells profile' created


png(filename = snakemake@output[["test_data_heatmap"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])
sc_hmClass(classMat = classRes_val_all,
           grps = sla,
           max = 300,
           isBig = TRUE)
dev.off()

# Attribution plot
print("Creating training data attribution plot")
png(filename = snakemake@output[["attr_plot"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])
plot_attr(classRes = classRes_val_all,
          sampTab = stTest,
          nrand = nrand,
          dLevel = "Cell_type",
          sid = "Cell_ID")
dev.off()

# Visualize average top pairs genes expression for training data
print("Creating training data top genepairs")
gpTab <- compareGenePairs(query_exp = expTest,
                          training_exp = expTrain,
                          training_st = stTrain,
                          classCol = "Cell_type",
                          sampleCol = "Cell_ID",
                          RF_classifier = class_info$cnProc$classifier,
                          numPairs = 20,
                          trainingOnly = TRUE)

train <- findAvgLabel(gpTab = gpTab,
                      stTrain = stTrain,
                      dLevel = "Cell_type")

png(filename = snakemake@output[["average_top_genepair"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])
hm_gpa_sel(gpTab,
           genes = class_info$cnProc$xpairs,
           grps = train,
           maxPerGrp = 50,
           fontsize_row = 7)
dev.off()

#---Predict Query Data ####
print("Classifying query data")
nqRand = 50
system.time(crQuery <- scn_predict(class_info[['cnProc']],
                                   queryExp,
                                   nrand = nqRand))

# Classification annotation assignment
# This classifies a cell with  the category with the highest classification score or higher than a classification score threshold of your choosing.
# The annotation result can be found in a column named category in the query sample table.
print("Assigning newly created categories")
queryMeta <- assign_cate(classRes = crQuery[,1:(ncol(crQuery)-nqRand)],
                         sampTab = queryMeta,
                         cThresh = 0.5)

# Query visualization
print("Creating query data heatmap")
if (is.null(snakemake@params[["Cell_description_colname"]])) {
  sgrp <- as.vector(queryMeta$category)
} else {
  sgrp <- as.vector(queryMeta$Cell_type)
}
names(sgrp) <- as.vector(queryMeta$Cell_ID)
grpRand <- rep("rand", nqRand)
names(grpRand) <- paste("rand_", 1:nqRand, sep='')
sgrp <- append(sgrp, grpRand)

# Heatmap classification result
png(filename = snakemake@output[["query_classification_heatmap"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])
sc_hmClass(crQuery, sgrp, max = 5000, isBig = TRUE, cCol = F, font = 8)
dev.off()

# Classification result violin plot
print("Creating query data violin plot")
png(filename = snakemake@output[["query_violin_plot"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])

if (is.null(snakemake@params[["Cell_description_colname"]])) {
  sc_violinClass(sampTab = queryMeta,
                 classRes = crQuery,
                 sid = "Cell_ID",
                 dLevel = "category",
                 ncol = length(unique(queryMeta$category)),
                 addRand = nqRand)
} else {
  sc_violinClass(sampTab = queryMeta,
                 classRes = crQuery,
                 sid = "Cell_ID",
                 dLevel = "Cell_type",
                 ncol = 12,
                 addRand = nqRand)
  
}
dev.off()

# Zoom in one cluster
print("Creating query data subcluster violin plot")
sub_cluster <- snakemake@params[["Cell_type_to_analyze"]]

png(filename = snakemake@output[["query_subcluster_violin_plot"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])
if (!is.null(sub_cluster)) {
  sc_violinClass(sampTab = queryMeta,
                 classRes = crQuery,
                 sid = "Cell_ID",
                 dLevel = "category",
                 ncol = 12,
                 addRand = nqRand,
                 sub_cluster = sub_cluster)
}

dev.off()

# UMAP by category
print("Creating UMAP")
system.time(umPrep_HS <- prep_umap_class(crQuery,
                                         queryMeta,
                                         nrand = nqRand,
                                         dLevel = "category",
                                         sid  = "Cell_ID",
                                         topPC = 5))

png(filename = snakemake@output[["query_umap"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])
plot_umap(umPrep_HS)
dev.off()

#---Score table and violin plot ####
print("Generating score figures")
# Applies score function
scoreTable <- score(expTrain = expTrain,
                    stTrain = stTrain,
                    expQuery = queryExp,
                    stQuery = queryMeta,
                    genePairs = class_info$cnProc$xpairs,
                    numPairs = 25,
                    importanceTable = randomForest::importance(class_info$cnProc$classifier))

# Score table
png(filename = snakemake@output[["score_heatmap"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])
pheatmap(scoreTable[[1]], main = "Score Heatmap")
dev.off()

# Query Gene expression heatmap
png(filename = snakemake@output[["score_expression_heatmap"]],
    width = snakemake@params[["fig_width"]],
    height = snakemake@params[["fig_height"]],
    res = snakemake@params[["fig_res"]])
pheatmap(scoreTable[[2]], main = "Query Gene Expression (relative)")
dev.off()