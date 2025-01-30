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

traits_otu <- fread('project/data/jin2017/traits_otupca.csv')

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


# Analysis grouped by cultivar --------------------------------------------

# Use only data with >1 observation per cultivar
multiple_obs_cultivars <- unique(traits_otu$cultivar[duplicated(traits_otu$cultivar)])

traits_otu_multipleobs <- traits_otu[cultivar %in% multiple_obs_cultivars]

library(lme4)
lmer_test <- lmer(log(grain_weight) ~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 + (1|cultivar), data=traits_otu_multipleobs)


# BRM model with no missing values ----------------------------------------

# Use similar horseshoe priors
library(brms)
library(rstan)
options(brms.backend = 'cmdstanr', brms.file_refit = 'on_change', mc.cores = 4)

n_PCs <- 10

# Construct priors programmatically because we have one for each taxon.
sd_X_priors <- lapply(1:n_PCs, function(i) prior_string('gamma(1, 1)', class = 'sd', resp = paste0('PC', i)))
sd_X_priors <- do.call(c, sd_X_priors)
sigma_X_priors <- lapply(1:n_PCs, function(i) prior_string('gamma(1, 1)', class = 'sigma', resp = paste0('PC', i)))
sigma_X_priors <- do.call(c, sigma_X_priors)

# Also construct formula programmatically.
X_formula <- paste0('mvbind(', paste(paste0('PC',1:n_PCs), collapse = ','), ') ~ 0 + (1||cultivar)')
y_formula <- paste0('grain_weight ~ ', paste(paste0('PC',1:n_PCs), collapse = '+'), ' + (1||cultivar)')

fit_grainweight_nomiss <- brm(
  bf(mvbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) ~ 0 + (1||cultivar)) + bf(grain_weight ~ PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 + (1||cultivar)) + set_rescor(FALSE),
  prior = c(
    sd_X_priors,
    sigma_X_priors,
    prior(gamma(1, 1), class = sd, resp = grainweight),
    prior(gamma(1, 1), class = sigma, resp = grainweight),
    prior(horseshoe(df = 1, df_global = 1, scale_slab = 20, df_slab = 10, par_ratio = 1), class = b, resp = grainweight) 
  ),
  data = traits_otu_multipleobs,
  chains = 4, iter = 7500, warmup = 5000,
  init = 0, seed = 1004,
  control = list(adapt_delta = 0.95),
  file = 'project/fits/brm_grainweight_nomiss'
)


# BRM model with missing values -------------------------------------------

# Artificially generate the missing data
# 50% of the individuals within each cultivar will have missing trait data, and 50% will have missing microbiome data
set.seed(444)

# Function to randomly sample half
samplehalf <- function(n) sample(rep(1:2, times=ceiling(n/2)))[1:n]

traits_otu_withmiss <- copy(traits_otu_multipleobs)
traits_otu_withmiss[, grp := samplehalf(.N) , by = .(cultivar)]

traits_otu_withmiss[grp == 1, grain_weight := NA]
PCcols <- paste0('PC', 1:n_PCs)
traits_otu_withmiss[ , (PCcols) := lapply(.SD, function(x) fifelse(grp == 2, NA_real_, x)), .SDcols = PCcols ]

# Construct formulas programmatically.
Xmiss_formula <- paste0('mvbind(', paste(PCcols, collapse = ','), ') | mi() ~ 0 + (1||cultivar)')
ymiss_formula <- paste0('grain_weight | mi() ~ ', paste(paste0('mi(PC', 1:n_PCs, ')'), collapse = '+'), ' + (1||cultivar)')

fit_grainweight_miss <- brm(
  bf(mvbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10) | mi() ~ 0 + (1||cultivar)) + bf(grain_weight | mi() ~ mi(PC1)+mi(PC2)+mi(PC3)+mi(PC4)+mi(PC5)+mi(PC6)+mi(PC7)+mi(PC8)+mi(PC9)+mi(PC10) + (1||cultivar)) + set_rescor(FALSE),
  prior = c(
    sd_X_priors,
    sigma_X_priors,
    prior(gamma(1, 1), class = sd, resp = grainweight),
    prior(gamma(1, 1), class = sigma, resp = grainweight),
    prior(horseshoe(df = 1, df_global = 1, scale_slab = 20, df_slab = 10, par_ratio = 1), class = b, resp = grainweight) 
  ),
  data = traits_otu_withmiss,
  chains = 4, iter = 7500, warmup = 5000,
  init = 0, seed = 1104,
  control = list(adapt_delta = 0.95),
  file = 'project/fits/brm_grainweight_miss'
)
