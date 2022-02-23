# Install required packages inside the Docker image. This may take a while.
print("Installing required packages for the pipeline. This may take a while.")
print("Installing remotes")
install.packages("remotes")
print("Remotes instalation completed")

# singleCellNet
print("Installing singleCellNet")
remotes::install_github("pcahan1/singleCellNet")
print("singleCellNet instalation completed")

# Seurat
print("Installing Seurat")
install.packages('Seurat')
print("Seurat instalation completed")

#Symphony
print("Installing symphony")
install.packages("symphony")
remotes::install_github("immunogenomics/singlecellmethods")
install.packages("tidyverse")
install.packages("ggthemes")
install.packages("ggpubr")
print("Symphony instalation completed")

#FUSCA
print("Installing FUSCA")
package_list <- c('reshape', 'reshape2', 'grid', 'ggplot2', 'pheatmap',
                 'igraph', 'mclust', 'Rtsne', 'cccd', 'irlba', 'ggrepel',
                 'dplyr', 'data.table', 'proxy', 'cluster', 'scales',
                 'Matrix', 'princurve', 'magrittr', 'Hmisc',
                 'RColorBrewer', 'ggnetwork', 'methods', 'umap', 'plyr',
                 'rJava', 'statmod', 'readbitmap',
                 'rjson', 'dbscan', 'tibble')
install.packages(package_list)

install.packages("BiocManager")
bioc_package_list <- c('GO.db', 'org.Hs.eg.db', 'org.Mm.eg.db', 'AnnotationDbi', 'bluster', 'BiocParallel', 'bluster')
BiocManager::install(bioc_package_list)

remotes::install_github("sctyner/geomnet")
install.packages(c('usethis', 'httr', 'rcmdcheck', 'roxygen2', 'rversions'))

remotes::install_github("edroaldo/fusca") #This needs to be tested

print("FUSCA instalation completed")
