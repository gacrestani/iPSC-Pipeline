library(GEOquery)
library(hash)

#Get data from GEO
getGEOSuppFiles('Data/GSE76381')
gse <- read.delim("Data/GSE76381/GSE76381_EmbryoMoleculeCounts.cef.txt.gz", header=FALSE)

#Create metadata and expression matrix
gse_meta <- gse[2:4, -1]
gse_meta <- as.data.frame(t(gse_meta))
colnames(gse_meta) <- gse_meta[1,]
rownames(gse_meta) <- gse_meta[,1]
gse_meta <- gse_meta[-1,]

# prototype dictionary
h <- hash()

h[['hEndo']]  <- 'Endo'
h[['hPeric']] <-  'Peric'
h[['hMgl']] <- 'Mgl'
h[['hDA1']] <- 'DA'
h[['hDA2']] <- 'DA'
h[['hDA0']] <- 'DA'
h[['hSert']] <- 'Sert'
h[['hOMTN']] <- 'OMTN'
h[['hRgl1']] <- 'Rgl'
h[['hRgl3']] <- 'Rgl'
h[['hRgl2c']] <- 'Rgl'
h[['hRgl2b']] <- 'Rgl'
h[['hRgl2a']] <- 'Rgl'
h[['hOPC']] <- 'OPC'
h[['hProgFPM']] <- 'ProgFP'
h[['hProgFPL']] <- 'ProgFP'
h[['hProgM']] <- 'ProgFP'
h[['hProgBP']] <- 'ProgBP'
h[['hNbML5']] <- 'Gaba'
h[['hGaba']] <- 'Gaba'
h[['hNbGaba']] <- 'Gaba'
h[['hNbML1']] <- 'NbML1'
h[['hNProg']] <- 'NProg'
h[['hNbM']] <- 'NbM'
h[['hRN']] <- 'RN'

for (i in 1:nrow(gse_meta)) {
  if (has.key(gse_meta$Cell_type[i], h)) {
    gse_meta$Cell_type[i] <- h[[gse_meta$Cell_type[i]]]
  }
}

gse_exp <- gse
rowNames <- gse_exp[5:nrow(gse_exp), 1]
gse_exp[5:nrow(gse_exp), 2] <- rowNames
gse_exp <- gse_exp[-1,-1]
colnames(gse_exp) <- gse_exp[1,]
rownames(gse_exp) <- gse_exp[,1]
gse_exp <- gse_exp[5:nrow(gse_exp), -1]
gse_exp <- as(as.matrix(gse_exp), 'dgCMatrix')

save(gse_exp, file='lamanno_human_exp.rda')
save(gse_meta, file='lamanno_human_meta.rda')


#Same thing for iPSC data
##Get data
gse <- read.delim("Data/GSE76381/GSE76381_iPSMoleculeCounts.cef.txt.gz", header=FALSE)

#Create metadata and expression matrix
gse_meta <- gse[2:4, -1]
gse_meta <- as.data.frame(t(gse_meta))
colnames(gse_meta) <- gse_meta[1,]
rownames(gse_meta) <- gse_meta[,1]
gse_meta <- gse_meta[-1,]

g <- hash()

g[['iDAa']] <- 'DA'
g[['iDAb']] <- 'DA'
g[['iDAc']] <- 'DA'
g[['iRgl1']] <- 'Rgl'
g[['iRgl2']] <- 'Rgl'
g[['iRN']] <- 'RN'

for (i in 1:nrow(gse_meta)) {
  if (has.key(gse_meta$Cell_type[i], g)) {
    gse_meta$Cell_type[i] <- g[[gse_meta$Cell_type[i]]]
  }
}

gse_exp <- gse
rowNames <- gse_exp[5:nrow(gse_exp), 1]
gse_exp[5:nrow(gse_exp), 2] <- rowNames
gse_exp <- gse_exp[-1,-1]
colnames(gse_exp) <- gse_exp[1,]
rownames(gse_exp) <- gse_exp[,1]
gse_exp <- gse_exp[5:nrow(gse_exp), -1]
gse_exp <- as(as.matrix(gse_exp), 'dgCMatrix')

save(gse_exp, file='lamanno_ips_exp.rda')
save(gse_meta, file='lamanno_ips_meta.rda')