# Read root trait and microbiome data, and fit model.

library(data.table)
library(mixOmics)

datadir <- 'project/data/doi_10_5061_dryad_5p414__v20190105/'
agb <- fread(file.path(datadir, 'Host_plant_drought_biomass.csv'))

seqtab.total<-readRDS(file.path(datadir, "Root_microbiome_seqtab.rds"))
# bac.tree<-read.tree(file= "16S_Root_Microbiome.tre")
# 
samples.out<-rownames(seqtab.total)
# 
meta.data<-fread(file=file.path(datadir, "meta.data.csv"))
# 
# taxa<-readRDS(file.path(datadir, "Root_microbiome_RDP_taxa.rds"))

# Perform some dimension reduction on the data
# Eliminate all the taxa that are present in <=10 samples
seqtab_gt10 <- seqtab.total[, colSums(seqtab.total>0) >= 10]

# Do CLR transformation, adding 1 to all counts. Replace original values with transformed.
seqtab_CLR <- logratio.transfo(seqtab_gt10, logratio = 'CLR', offset = 1)
seqtab_CLR <- as(seqtab_CLR, 'matrix')

# Change column names to short ID
dimnames(seqtab_CLR)[[2]] <- paste0('seq', 1:ncol(seqtab_CLR))

# Scale the columns
# seqtab_gt10_scale <- scale(seqtab_gt10)

seqtab_pca <- prcomp(seqtab_CLR)
summ_pca <- summary(seqtab_pca)

# Join the transformed sequence table with the trait data.
