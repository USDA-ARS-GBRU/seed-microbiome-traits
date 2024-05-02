# Ben's multilevel model

library(data.table)
library(mixOmics)
library(brms)

options(mc.cores=4, brms.backend='cmdstanr', brms.file_refit='on_change')

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
# Heterogeneous vector of data.

# Use row 1 of the CLR transformed matrix.
data_taxon1 <- data.table(plantID = names(CLR_16S)[-1], taxon1abundance = as.numeric(CLR_16S[1, -1]))
data_taxon1[, maternal_plant_code := substr(plantID, 3, 3)]
data_taxon1[, country_origin := ifelse(maternal_plant_code %in% LETTERS[1:5], 'Afghanistan', 'Turkey')]

traits[, plantID := paste0('t', 1:10, maternal_plant_code)]

root_trait_taxon1 <- rbind(
  traits[!is.na(rooting_depth_cm), .(plantID, maternal_plant_code, country_origin, rooting_depth_cm)],
  data_taxon1,
  use.names = TRUE, fill = TRUE
)

root_trait_taxon1[, rooting_depth_std := as.vector(scale(rooting_depth_cm))]
root_trait_taxon1[, y := fcoalesce(rooting_depth_std, taxon1abundance)]
root_trait_taxon1[, trt := ifelse(is.na(rooting_depth_cm), 'community', 'trait')]

# Ben's multilevel model

mlmod_root_taxon1 <- brm(
  bf(
    y ~ country_origin * trt + (1 | maternal_plant_code),
    sigma ~ trt
  ),
  family = gaussian, data = root_trait_taxon1,
  chains = 1, iter = 2000, warmup = 1000,
  file = 'project/mlmod_root_taxon1'
)
