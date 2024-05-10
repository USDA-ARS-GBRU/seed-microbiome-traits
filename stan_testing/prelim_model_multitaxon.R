# Multitaxon single trait model in Stan

library(data.table)
library(mixOmics)
library(cmdstanr)

options(mc.cores=4)

# Import data. Corrected labels, data accessed from g drive 2023-06-06
abundance16S <- fread('data/16S_abundance.csv')
abundanceITS <- fread('data/ITS_Abundance.csv')
traits <- fread('data/seedling_data_updated.csv')

CLR_abundance <- function(abundance) {
  # Remove maternal plant J as it has no seedling traits.
  J_columns <- grep('J$', names(abundance), value = TRUE)
  abundance[, c(J_columns) := NULL]
  
  # Remove any taxa that now have all zero abundance in the remaining 9 maternal plants.
  zero_rows <- rowSums(abundance[,-1]) == 0
  abundance <- abundance[!zero_rows]
  
  # Do CLR transformation, adding 1 to all counts. Replace original values with transformed.
  abundance_CLR <- logratio.transfo(abundance[,-1], logratio = 'CLR', offset = 1)
  abundance_CLR <- cbind(abundance[, .(V1)], as(abundance_CLR, 'matrix'))
  
}

CLR_16S <- CLR_abundance(abundance16S)
CLR_16S_array <- t(CLR_16S[,-1])

# Index which maternal plants and populations go with each row of the abundance dataset.
CLR_16S_metadata <- data.table(plantID = names(CLR_16S)[-1])
CLR_16S_metadata[, maternal_plant_code := substr(plantID, 3, 3)]
CLR_16S_metadata[, country_origin := ifelse(maternal_plant_code %in% LETTERS[1:5], 'Afghanistan', 'Turkey')]

traits[, plantID := paste0('t', 1:10, maternal_plant_code)]

data_root <- traits[!is.na(rooting_depth_cm), .(plantID, maternal_plant_code, country_origin, rooting_depth_cm)]
data_root[, rooting_depth_std := as.vector(scale(rooting_depth_cm))]

maternal_plants <- unique(data_root[, .(maternal_plant_code, country_origin)])

### Here take a random sample of 10 taxa for testing.
set.seed(131)
CLR_16S_sample <- CLR_16S_array[, sample(1:ncol(CLR_16S_array), size = 10)]

stan_data_root_16S_test <- list(
  Nmicro = nrow(CLR_16S_sample),
  Ntrait = nrow(data_root),
  M = nrow(maternal_plants),
  D = ncol(CLR_16S_sample),
  y = data_root$rooting_depth_std,
  X = CLR_16S_sample,
  pop = as.numeric(maternal_plants$country_origin == 'Turkey'),
  Mmicro = as.numeric(factor(CLR_16S_metadata$maternal_plant_code)),
  Mtrait = as.numeric(factor(data_root$maternal_plant_code))
)

multitaxonmodel <- cmdstan_model(stan_file = 'stan_testing/onetraitmultitaxonmodel.stan')
multitaxonmodelopt <- cmdstan_model(stan_file = 'stan_testing/onetraitmultitaxonmodeloptimized.stan')

fit_root_16S_test <- multitaxonmodel$sample(
  data = stan_data_root_16S_test,
  seed = 27606,
  chains = 4,
  iter_warmup = 2500,
  iter_sampling = 2500,
  init = 0.1 # Attempt to deal with initial value of correlation matrix not being positive definite.
)

fit_root_16S_test2 <- multitaxonmodelopt$sample(
  data = stan_data_root_16S_test,
  seed = 27606,
  chains = 4,
  iter_warmup = 2500,
  iter_sampling = 2500
)

summ_root_16S_test <- fit_root_16S_test$summary()
summ_root_16S_test2 <- fit_root_16S_test2$summary()
