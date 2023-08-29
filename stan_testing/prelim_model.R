# Preliminary exploration of seed microbiome and plant trait data

library(data.table)
library(cmdstanr)

# Import data. Corrected labels, data accessed from g drive 2023-06-06
abundance <- fread('data/foxx/data/16S_abundance.csv')
traits <- fread('data/foxx/data/seedling_data_updated.csv')

# Here, we will try to fit a "toy model" using a few haphazardly selected microbial taxa abundances, and one continuous plant trait (rooting depth) as the response variable
# Currently, for simplicity I will convert the abundance matrix to a presence absence matrix. 
# I will model the presence-absences of all the taxa for a given seed as a draw from a multivariate probit distribution, where each maternal plant has its own parameters for that distribution.
# I will model the (log) rooting depths of the offspring of a given maternal plant as draws from a normal distribution.
# Then, we'll estimate the effects of each taxon abundance on rooting depth, using the parameters of the above distributions.
# In both cases maternal plant identity will be a random effect, and country of origin will be a fixed effect

# Convert microbiome data to presence-absence
pa <- 1 * (abundance[,-1] > 0)

# Exclude maternal plant I from the 

# Which maternal plants do each of the microbiome samples come from?
# last character of the column names of the matrix
maternal_plants_microbiome <- substr(colnames(pa), nchar(colnames(pa)), nchar(colnames(pa)))

# Build the Stan data
# We can only use the ones that have rooting depth

traits_formodel <- traits[!is.na(rooting_depth_cm)]

maternal_plants <- unique(traits[, .(maternal_plant, maternal_plant_code, country_origin)])

stan_micro_trait_data <- list(
  Nmicro = ncol(pa),
  Ntrait = nrow(traits_formodel),
  M = nrow(maternal_plants),
  K = nrow(pa),
  y = traits_formodel$rooting_depth_cm,
  X = pa,
  pop = as.numeric(maternal_plants$country_origin == 'Turkey'),
  Mmicro = as.numeric(factor(maternal_plants_microbiome)),
  Mtrait = as.numeric(factor(traits_formodel$maternal_plant_code))
)

onetraitmodel <- cmdstan_model(stan_file = 'seed_microbiome_traits/onetraitmodel.stan')
