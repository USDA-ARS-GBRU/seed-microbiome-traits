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

standata_miss5 <- make_standata_vectorized_horseshoe(dt = dt, x_names = grep('Xmiss', names(dt), value = TRUE), y_name = 'ymiss', grp_name = 'maternal_id', hs_df = 1, hs_df_global = 1, hs_df_slab = 4, hs_scale_slab = 20, hs_scale_global = (2/3)/sqrt(nrow(dt)), prior_sd_Intercept_X = 2.5, prior_sd_Intercept_Y = 10)

horsemod <- cmdstan_model('stan_testing/miss_horseshoe_vectorized.stan')

fit_miss5 <- horsemod$sample(data = standata_miss5, seed = 111, chains = 4, parallel_chains = 4, iter_warmup = 1000, iter_sampling = 1000)
