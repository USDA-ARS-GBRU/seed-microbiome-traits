# Import and process data from Jin et al. 2017

library(data.table)
library(mixOmics)

# Import traits and join with metadata (doesn't seem like metadata will be that useful)

traits <- fread('project/data/jin2017/Traits - Jin et al. 2017.csv', skip = 2)
trait_names <- c('plantID', 'top_second_leaf_length', 'top_second_leaf_width', 'main_stem_height', 'main_stem_width', 'panicle_length', 'panicle_diameter', 'fringe_neck_length', 'panicle_weight', 'grain_weight', 'hundred_kernel_weight', 'spikelet_number', 'grain_number')
setnames(traits, trait_names)

metadata <- fread('project/data/jin2017/Metadata - Jin et al. 2017.csv', skip = 1)
setnames(metadata, c('plantID', 'cultivar', 'location', 'compartment', 'number'))
metadata[cultivar == '', cultivar := NA]
metadata[, cultivar := zoo::na.locf(cultivar)]

traits <- merge(traits, metadata, by = 'plantID', all.x = TRUE)

# Import OTU table
otu <- fread('project/data/jin2017/OTU table - Jin et al. 2017.csv', skip = 1)

# Subset the OTU table for only the individuals that we have trait data for
otu_use <- otu[, mget(names(otu)[names(otu) %in% c('#OTU ID', 'taxonomy', traits$plantID)])]

otu_mat <- as.matrix(otu_use[, mget(setdiff(names(otu_use), c('#OTU ID', 'taxonomy')))])
dimnames(otu_mat) <- list(otu_use[['#OTU ID']], setdiff(names(otu_use), c('#OTU ID', 'taxonomy')))

# Transpose to a plant by taxon matrix (site by species)
otu_mat <- t(otu_mat)

# Perform CLR transform on the OTU data and then do a dimension reduction with sPCA
otu_clr <- logratio.transfo(otu_mat, logratio = 'CLR', offset = 1)

otu_spca_10x50 <- spca(otu_clr, ncomp = 10, keepX = rep(50, 10))  # run the method
plotIndiv(otu_spca_10x50)  # plot the samples
plotVar(otu_spca_10x50)    # plot the variables

# extract the variables used to construct the first PC
selectVar(otu_spca_10x50, comp = 1)$name 
# depict weight assigned to each of these variables
plotLoadings(otu_spca_10x50, method = 'mean', contrib = 'max') 

# Join the OTU PCA and trait data together
otu_pca_loading_df <- data.frame(plantID = row.names(otu_spca_10x50$x), otu_spca_10x50$x)
traits_otu <- traits[otu_pca_loading_df, on = 'plantID']

fwrite(traits_otu, 'project/data/jin2017/traits_otupca.csv')

# Test relationship with traits -------------------------------------------

# Determine which traits have the highest correlation with the PCA axes
cormat <- cor(traits_otu[, mget(c(trait_names[-1], paste0('PC', 1:10)))], use = 'pairwise.complete.obs')[trait_names[-1], paste0('PC', 1:10)]
sort(apply(abs(cormat), 1, mean))

# Grain weight has one of the highest correlations

lm_test <- lm(log(grain_weight) ~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=traits_otu)

# Try regularization. We do get some nonzero values.
library(glmnet)
set.seed(333)
glmnet_test <- cv.glmnet(x = as.matrix(traits_otu[, mget(paste0('PC', 1:10))]), y = log(traits_otu$grain_weight), alpha = 1, nfolds = 10, type.measure = 'mse',
          lambda = seq(0.001, 0.1, by = 0.001))

