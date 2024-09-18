### Function to create stan data list for the missing data vectorized horseshoe Stan model
### Requires entire rows of microbiome data to be missing, not individual values in each row
### Only supports univariate trait response
make_standata_vectorized_horseshoe <- function(dt, x_names, y_name, grp_name, hs_df, hs_df_global, hs_df_slab, hs_scale_global, hs_scale_slab, prior_sd_Intercept_X, prior_sd_Intercept_Y) {
  N <- nrow(dt)
  N_taxa <- length(x_names)
  X <- as.list(dt[, x_names])
  Nmi_X <- sum(is.na(X[[1]]))
  Jmi_X <- which(is.na(X[[1]]))
  
  for (i in 1:length(X)) X[[i]][is.na(X[[i]])] <- Inf
  
  Y <- dt[[y_name]]
  Nmi_Y <- sum(is.na(Y))
  Jmi_Y <- which(is.na(Y))
  
  Y[is.na(Y)] <- Inf
  
  Ksp_Y <- N_taxa
  Kscales_Y <- N_taxa
    
  J <- as.numeric(dt[[grp_name]])
  N_grps <- length(unique(J))
  
  Z_X <- replicate(n = N_taxa, rep(1, N), simplify = FALSE)
  Z_Y <- rep(1, N)
  
  prior_only <- 0
  
  list(N = N, N_taxa = N_taxa, X = X, Nmi_X = Nmi_X, Jmi_X = Jmi_X, Y = Y, Nmi_Y = Nmi_Y, Jmi_Y = Jmi_Y, Ksp_Y = Ksp_Y, Kscales_Y = Kscales_Y, J = J, N_grps = N_grps, Z_X = Z_X, Z_Y = Z_Y, hs_df = hs_df, hs_df_global = hs_df_global, hs_df_slab = hs_df_slab, hs_scale_global = hs_scale_global, hs_scale_slab = hs_scale_slab, prior_sd_Intercept_X = prior_sd_Intercept_X, prior_sd_Intercept_Y = prior_sd_Intercept_Y, prior_only = prior_only)
  
}

library(cmdstanr)
library(mvtnorm)

# Start with manageable number of taxa.
n_mothers <- 10
n_taxa <- 5
offspring_per_mother <- 20 # half will be retained for traits, half for microbiome

set.seed(1)

X_maternal <- rmvnorm(n_mothers, mean = rep(0, n_taxa), sigma = diag(n_taxa))
sigma_maternal <- cov(X_maternal)

# Coefficients indicating which taxa predict the outcome.
# We will not include any interaction effect.
beta <- c(5, -5, 0, 0, 0)

y_maternal <- 10 + X_maternal %*% beta + rnorm(n_mothers, 0, 1)

# To get offspring microbiome, take multivariate normal draws from the mean vector for each mother (rows of X_maternal)
X_offspring <- apply(X_maternal, 1, function(Xi) rmvnorm(offspring_per_mother, mean = Xi, sigma = sigma_maternal), simplify = FALSE)

# Use regression coefficients (beta) to get value for offspring trait, plus noise
y_offspring <- lapply(X_offspring, function(Xoi) Xoi %*% beta + rnorm(offspring_per_mother, 0, 1))

# Combine together
dt <- data.frame(
  maternal_id = factor(rep(1:n_mothers, each = offspring_per_mother)),
  offspring_id = 1:offspring_per_mother,
  do.call(rbind, X_offspring),
  y = do.call(c, y_offspring)
)

# Within each mother, set half of the values to be missing for x, and the other half for y.
xmiss <- lapply(1:nrow(dt), function(i) {
  if (dt[i, 'offspring_id'] %in% 1:(offspring_per_mother/2)) {
    setNames(dt[i, paste0('X', 1:n_taxa)], paste0('Xmiss',1:n_taxa))
  } else {
    setNames(rep(NA, n_taxa), paste0('Xmiss', 1:n_taxa))
  }
})
dt <- cbind(dt, do.call(rbind, xmiss))
dt$ymiss <- ifelse(dt$offspring_id %in% 1:(offspring_per_mother/2), NA, dt$y)


standata_miss5 <- make_standata_vectorized_horseshoe(dt = dt, x_names = grep('Xmiss', names(dt), value = TRUE), y_name = 'ymiss', grp_name = 'maternal_id', hs_df = 1, hs_df_global = 1, hs_df_slab = 4, hs_scale_slab = 20, hs_scale_global = (2/3)/sqrt(nrow(dt)), prior_sd_Intercept_X = 2.5, prior_sd_Intercept_Y = 10)

horsemod <- cmdstan_model('stan_testing/miss_horseshoe_vectorized.stan')

fit_miss5 <- horsemod$sample(data = standata_miss5, seed = 333, 
                             chains = 4, parallel_chains = 4, 
                             iter_warmup = 2500, iter_sampling = 2500,
                             adapt_delta = 0.9, max_treedepth = 15)
