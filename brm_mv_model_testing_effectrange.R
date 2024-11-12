# Simulate multivariate normal data in two levels: maternal taxa abundance means are MVN, 
# then offspring taxa abundance means are MVN from that.
# Update 2024-09-06: Do not estimate separate intercept for each taxon

library(mvtnorm)
library(brms)

options(mc.cores = 4, brms.backend = 'cmdstanr', brms.file_refit = 'on_change')

# Increase number of taxa
n_mothers <- 20
n_taxa <- 200
offspring_per_mother <- 10 # 5 will be retained for traits, 5 for microbiome

set.seed(2)

X_maternal <- rmvnorm(n_mothers, mean = rep(0, n_taxa), sigma = diag(n_taxa))
sigma_maternal <- cov(X_maternal)

# Coefficients indicating which taxa predict the outcome.
# We will not include any interaction effect.
# Only include a few taxa with a nonzero effect.
beta <- c(50, 20, 10, 5, 2, 1, rep(0, n_taxa - 6))

y_maternal <- 0 + X_maternal %*% beta + rnorm(n_mothers, 0, 1)

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


# Model without missing data ----------------------------------------------

# Model without missing data, including regularized horseshoe prior
# Let's see if the coefficients can be recovered.
# Interactions between taxa aren't included.

# Construct priors programmatically because we have one for each taxon.
sd_X_priors <- lapply(1:n_taxa, function(i) prior_string('gamma(1, 1)', class = 'sd', resp = paste0('X', i)))
sd_X_priors <- do.call(c, sd_X_priors)
sigma_X_priors <- lapply(1:n_taxa, function(i) prior_string('gamma(1, 1)', class = 'sigma', resp = paste0('X', i)))
sigma_X_priors <- do.call(c, sigma_X_priors)

# Also construct formula programmatically.
X_formula <- paste0('mvbind(', paste(paste0('X',1:n_taxa), collapse = ','), ') ~ 0 + (1||maternal_id)')
y_formula <- paste0('y ~ ', paste(paste0('X',1:n_taxa), collapse = '+'), ' + (1||maternal_id)')

modmv_nomiss_reghorseshoe <- brm(
  bf(X_formula) + bf(y_formula) + set_rescor(FALSE),
  prior = c(
    sd_X_priors,
    sigma_X_priors,
    prior(gamma(1, 1), class = sd, resp = y),
    prior(gamma(1, 1), class = sigma, resp = y),
    prior(horseshoe(df = 1, df_global = 1, scale_slab = 50, df_slab = 4, par_ratio = 6/(n_taxa-6)), class = b, resp = y) 
  ),
  data = dt,
  chains = 4, iter = 7500, warmup = 5000,
  init = 0, seed = 1240,
  control = list(adapt_delta = 0.95),
  file = 'project/fits/brm_mv_nomiss_reghorseshoe'
)


# Model with missing data -------------------------------------------------

# With regularized horseshoe prior on fixed effects.

# Construct sd and sigma priors programmatically.
sd_Xmiss_priors <- lapply(1:n_taxa, function(i) prior_string('gamma(1, 1)', class = 'sd', resp = paste0('Xmiss', i)))
sd_Xmiss_priors <- do.call(c, sd_Xmiss_priors)
sigma_Xmiss_priors <- lapply(1:n_taxa, function(i) prior_string('gamma(1, 1)', class = 'sigma', resp = paste0('Xmiss', i)))
sigma_Xmiss_priors <- do.call(c, sigma_Xmiss_priors)

# Also construct formula programmatically.
Xmiss_formula <- paste0('mvbind(', paste(paste0('Xmiss', 1:n_taxa), collapse = ','), ') | mi() ~ 0 + (1||maternal_id)')
ymiss_formula <- paste0('ymiss | mi() ~ ', paste(paste0('mi(Xmiss', 1:n_taxa, ')'), collapse = '+'), ' + (1||maternal_id)')

modmv_miss_reghorseshoe <- brm(
  bf(Xmiss_formula) + bf(ymiss_formula) + set_rescor(FALSE),
  prior = c(
    sd_Xmiss_priors,
    sigma_Xmiss_priors,
    prior(gamma(1, 1), class = sd, resp = ymiss),
    prior(gamma(1, 1), class = sigma, resp = ymiss),
    prior(horseshoe(df = 1, df_global = 1, scale_slab = 50, df_slab = 4, par_ratio = 6/(n_taxa-6)), class = b, resp = ymiss)
  ),
  data = dt,
  chains = 4, iter = 7500, warmup = 5000,
  init = 0, seed = 1239,
  control = list(adapt_delta = 0.95),
  file = 'project/fits/brm_mv_miss_reghorseshoe'
)

# Export summaries to download locally ------------------------------------

summ_nomiss <- summary(modmv_nomiss_reghorseshoe)
summ_miss <- summary(modmv_miss_reghorseshoe)

save(summ_nomiss, summ_miss, file = 'project/fits/brm_mv_summaries.RData')
