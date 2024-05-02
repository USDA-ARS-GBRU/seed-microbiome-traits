# Single trait single taxon model in Stan

library(data.table)
library(mixOmics)
library(cmdstanr)

options(mc.cores=4)

# Import data. Corrected labels, data accessed from g drive 2023-06-06
abundance16S <- fread('project/data/16S_abundance.csv')
abundanceITS <- fread('project/data/ITS_Abundance.csv')
traits <- fread('project/data/seedling_data_updated.csv')

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

# Set up data for a single taxon and trait. 

# Use row 1 of the CLR transformed matrix.
data_taxon1 <- data.table(plantID = names(CLR_16S)[-1], taxon1abundance = as.numeric(CLR_16S[1, -1]))
data_taxon1[, maternal_plant_code := substr(plantID, 3, 3)]
data_taxon1[, country_origin := ifelse(maternal_plant_code %in% LETTERS[1:5], 'Afghanistan', 'Turkey')]

traits[, plantID := paste0('t', 1:10, maternal_plant_code)]

data_root <- traits[!is.na(rooting_depth_cm), .(plantID, maternal_plant_code, country_origin, rooting_depth_cm)]
data_root[, rooting_depth_std := as.vector(scale(rooting_depth_cm))]

maternal_plants <- unique(data_root[, .(maternal_plant_code, country_origin)])

stan_data_root_taxon1 <- list(
  Nmicro = nrow(data_taxon1),
  Ntrait = nrow(data_root),
  M = length(unique(data_taxon1$maternal_plant_code)),
  y = data_root$rooting_depth_std,
  x = data_taxon1$taxon1abundance,
  pop = as.numeric(maternal_plants$country_origin == 'Turkey'),
  Mmicro = as.numeric(factor(data_taxon1$maternal_plant_code)),
  Mtrait = as.numeric(factor(data_root$maternal_plant_code))
)

onetraitonetaxonmodel <- cmdstan_model(stan_file = 'stan_testing/onetraitonetaxonmodel.stan')

fit <- onetraitonetaxonmodel$sample(
  data = stan_data_root_taxon1,
  seed = 27606,
  chains = 4,
  iter_warmup = 2500,
  iter_sampling = 2500
)
