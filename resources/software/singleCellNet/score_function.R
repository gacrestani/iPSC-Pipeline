#Score Function: returns a table sample x genePairs of an impact score for each gene pair for that specific classification

# # For debug only:
# expTrain = expTrain
# stTrain = stTrain
# expQuery = queryExp
# stQuery = queryMeta
# genePairs = class_info$cnProc$xpairs
# numPairs = 20
# importanceTable = randomForest::importance(class_info$cnProc$classifier)


score <- function(expTrain,
                  stTrain,
                  expQuery,
                  stQuery,
                  genePairs,
                  numPairs = 20,
                  importanceTable) {
  
  # Sort importanceTable and get top 20 genePairs
  importanceTable = sort(importanceTable[, 1], decreasing = TRUE)
  importanceTable = importanceTable[grep("_", names(importanceTable))]
  
  if (numPairs > length(importanceTable)) {
    userPairs = importanceTable
  } else {
    userPairs = importanceTable[1:numPairs]
  }
  
  # Coerce data structures and remove "_" from names
  if (class(stTrain) != "data.frame") {
    stTrain <- as.data.frame(stTrain)
  }
  
  if (class(stQuery) != "data.frame") {
    stQuery <- as.data.frame(stQuery)
  }
  
  x = grep("_", rownames(expTrain))
  if (is.integer(x) && length(x) != 0L) {
    cat("converting _ in gene names to .\n")
    rownames(expTrain) = gsub("_", ".", rownames(expTrain))
  }
  
  y = grep("_", rownames(expQuery))
  if (is.integer(y) && length(y) != 0L) {
    cat("converting _ in gene names to .\n")
    rownames(expQuery) = gsub("_", ".", rownames(expQuery))
  }
  
  # Get avgCounts tables
  avgCountsQuery <- getAvgCounts(queryExp = expQuery,
                                 queryMeta = queryMeta,
                                 Cell_type = "category")
  
  avgCountsTrain <- getAvgCounts(queryExp = expTrain,
                                 queryMeta = stTrain,
                                 Cell_type = "Cell_type")
  
  scale_avgCountsQuery <- scaleMatrix(avgCountsQuery)
  print("Query dim: ")
  print(dim(scale_avgCountsQuery))
  scale_avgCountsTrain <- scaleMatrix(avgCountsTrain)
  print("Train dim: ")
  print(dim(scale_avgCountsTrain))
  
  # Filter for genes in top genepairs
  genelist <- unique(unlist(strsplit(names(userPairs), "_")))
  genelist_dataQuery <- scale_avgCountsQuery[rownames(scale_avgCountsQuery) %in% genelist, ]
  genelist_dataTrain <- scale_avgCountsTrain[rownames(scale_avgCountsTrain) %in% genelist, ]
  
  genelist_dataQuery <- genelist_dataQuery[,order(colnames(genelist_dataQuery))]
  genelist_dataTrain <- genelist_dataTrain[,order(colnames(genelist_dataTrain))]
  
  genelist_dataQuery <- as.data.frame(genelist_dataQuery)
  genelist_dataQuery$rand.Avg <- NULL
  
  ratio_train <- genelist_dataTrain[ , colnames(genelist_dataTrain) %in% colnames(genelist_dataQuery)]
  print("Ratio train dim:")
  print(dim(ratio_train))
  
  print("Genelist DataQuery dim:")
  print(dim(genelist_dataQuery))
  
  genelist_dataRatio <- ratio_train * genelist_dataQuery
  print("Ratio dim: ")
  print(dim(genelist_dataRatio))
  
  # Creates scoreTable and asssign score values
  scoreTable <- as.data.frame(matrix(0, nrow = length(userPairs), ncol=ncol(genelist_dataRatio)))
  rownames(scoreTable) <- names(userPairs)
  colnames(scoreTable) <- colnames(genelist_dataQuery)

  for (i in 1:ncol(scoreTable)) {
    celltype <- colnames(scoreTable[i])
    for (j in 1:nrow(scoreTable)) {
      genePair <- rownames(scoreTable[j,])
      gene1 <- unlist(strsplit(genePair, "_"))[1]
      gene2 <- unlist(strsplit(genePair, "_"))[2]
      
      # expGene1 <- avgCountsQuery[gene1, celltype]
      # expGene2 <- avgCountsQuery[gene2, celltype]
      # 
      # train_expGene1 <- avgCountsTrain[gene1, celltype]
      # train_expGene2 <- avgCountsTrain[gene2, celltype]
      # 
      # if (train_expGene1 == 0) {
      #   gene1Ratio <- 0
      # } else {
      #   gene1Ratio <- expGene1/train_expGene1
      # }
      # 
      # if (train_expGene2 == 0) {
      #   gene2Ratio <- 0
      # } else {
      #   gene2Ratio <- expGene2/train_expGene2
      # }
      
      importance <- importanceTable[genePair]
      
      calculate <- log2(1 + importance * (genelist_dataRatio[gene1, celltype] + genelist_dataRatio[gene2, celltype]))
      
      if (is.nan(calculate)) {calculate = 0}
      
      scoreTable[j,i] <- calculate
      }
    }
  return(list(scoreTable, genelist_dataQuery))
}

# pheatmap(scoreTable, cluster_rows = F, cluster_cols = F, main = "Score")
# 
# pheatmap(genelist_dataQuery, cluster_rows = F, cluster_cols = F, main = "Query")
# pheatmap(temp, cluster_rows = F, cluster_cols = F, main = "Train")
# pheatmap(genelist_dataRatio, cluster_rows = F, cluster_cols = F, main = "Ratio")


# Auxiliar function - get average values of expression for each gene, for each celltype
getAvgCounts <- function(
  queryExp,
  queryMeta,
  Cell_type
) {
  colnames(queryMeta)[colnames(queryMeta) == Cell_type] <- "category"
  queryExp <- as.data.frame(as.matrix(queryExp))
  celltypes <- unique(queryMeta$"category")
  x <- as.data.frame(matrix(0, nrow = nrow(queryExp), ncol=0))
  for (i in celltypes) {
    queryMeta_aux <- queryMeta[queryMeta["category"] == i,]
    queryExp_aux <- as.data.frame(queryExp[,colnames(queryExp) %in% queryMeta_aux$Cell_ID])
    average_exp <- as.data.frame(apply(queryExp_aux, 1, mean))
    colnames(average_exp) <- paste0(i, ".Avg")
    x <- cbind(x, average_exp)
  }
  return(x)
}

scaleMatrix <- function(x) {
  data <- as.matrix(x)
  data <- t(data)
  data <- scale(data)
  data <- t(data)
  return(data)
}
