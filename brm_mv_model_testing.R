# Simulate multivariate normal data in two levels: maternal taxa abundance means are MVN, 
# then offspring taxa abundance means are MVN from that.

library(mvtnorm)
library(brms)

options(mc.cores = 4, brms.backend = 'cmdstanr', brms.file_refit = 'on_change')

# Start with manageable number of taxa.
n_mothers <- 10
n_taxa <- 5
offspring_per_mother <- 10 # 5 will be retained for traits, 5 for microbiome

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


# Model without missing data ----------------------------------------------

# Model without missing data, including regularized horseshoe prior
# Let's see if the coefficients can be recovered.
# Interactions between taxa aren't included.

get_prior(
  bf(mvbind(X1, X2, X3, X4, X5) ~ (1||maternal_id)) + bf(y ~ X1 + X2 + X3 + X4 + X5 + (1||maternal_id)) + set_rescor(FALSE),
  data=dt
)

# First try without regularization prior.
modmv_nomiss <- brm(
  bf(mvbind(X1, X2, X3, X4, X5) ~ (1||maternal_id)) + bf(y ~ X1 + X2 + X3 + X4 + X5 + (1||maternal_id)) + set_rescor(FALSE),
  prior = c(
    prior(gamma(1, 1), class = sd, resp = X1), 
    prior(gamma(1, 1), class = sd, resp = X2), 
    prior(gamma(1, 1), class = sd, resp = X3), 
    prior(gamma(1, 1), class = sd, resp = X4), 
    prior(gamma(1, 1), class = sd, resp = X5), 
    prior(gamma(1, 1), class = sd, resp = y),
    prior(gamma(1, 1), class = sigma, resp = X1), 
    prior(gamma(1, 1), class = sigma, resp = X2), 
    prior(gamma(1, 1), class = sigma, resp = X3), 
    prior(gamma(1, 1), class = sigma, resp = X4), 
    prior(gamma(1, 1), class = sigma, resp = X5), 
    prior(gamma(1, 1), class = sigma, resp = y),
    prior(normal(0, 5), class = b, resp = y) # Only y carries any fixed effects.
  ),
  data = dt,
  chains = 4, iter = 2000, warmup = 1000,
  file = 'project/fits/brmtest_mv_nomiss'
)

# Second try, with regularization prior
modmv_nomiss_reghorseshoe <- brm(
  bf(mvbind(X1, X2, X3, X4, X5) ~ (1||maternal_id)) + bf(y ~ X1 + X2 + X3 + X4 + X5 + (1||maternal_id)) + set_rescor(FALSE),
  prior = c(
    prior(gamma(1, 1), class = sd, resp = X1), 
    prior(gamma(1, 1), class = sd, resp = X2), 
    prior(gamma(1, 1), class = sd, resp = X3), 
    prior(gamma(1, 1), class = sd, resp = X4), 
    prior(gamma(1, 1), class = sd, resp = X5), 
    prior(gamma(1, 1), class = sd, resp = y),
    prior(gamma(1, 1), class = sigma, resp = X1), 
    prior(gamma(1, 1), class = sigma, resp = X2), 
    prior(gamma(1, 1), class = sigma, resp = X3), 
    prior(gamma(1, 1), class = sigma, resp = X4), 
    prior(gamma(1, 1), class = sigma, resp = X5), 
    prior(gamma(1, 1), class = sigma, resp = y),
    prior(horseshoe(df = 1, df_global = 1, scale_slab = 20, df_slab = 4, par_ratio = 2/3), class = b, resp = y) # Only y carries any fixed effects. 
    # We are setting par_ratio to 2/3 because we "know" that is the ratio of nonzero to zero we are going for. Otherwise default values are used.
  ),
  data = dt,
  chains = 4, iter = 2000, warmup = 1000,
  file = 'project/fits/brmtest_mv_nomiss_reghorseshoe'
)


# Model with missing data -------------------------------------------------

# First attempt. No regularization prior is used.
modmv_miss <- brm(
  bf(mvbind(Xmiss1, Xmiss2, Xmiss3, Xmiss4, Xmiss5) | mi() ~ (1|p|maternal_id)) + bf(ymiss | mi() ~ mi(Xmiss1) + mi(Xmiss2) + mi(Xmiss3) + mi(Xmiss4) + mi(Xmiss5) + (1|maternal_id)) + set_rescor(TRUE),
  prior = c(
    prior(gamma(1, 1), class = sd, resp = Xmiss1), 
    prior(gamma(1, 1), class = sd, resp = Xmiss2), 
    prior(gamma(1, 1), class = sd, resp = Xmiss3), 
    prior(gamma(1, 1), class = sd, resp = Xmiss4), 
    prior(gamma(1, 1), class = sd, resp = Xmiss5), 
    prior(gamma(1, 1), class = sd, resp = ymiss),
    prior(gamma(1, 1), class = sigma, resp = Xmiss1), 
    prior(gamma(1, 1), class = sigma, resp = Xmiss2), 
    prior(gamma(1, 1), class = sigma, resp = Xmiss3), 
    prior(gamma(1, 1), class = sigma, resp = Xmiss4), 
    prior(gamma(1, 1), class = sigma, resp = Xmiss5), 
    prior(gamma(1, 1), class = sigma, resp = ymiss),
    prior(normal(0, 5), class = b, resp = ymiss) # Only y carries any fixed effects.
  ),
  data = dt,
  chains = 4, iter = 2000, warmup = 1000,
  file = 'project/fits/brmtest_mv_miss'
)

# With regularized horseshoe prior on fixed effects.
modmv_miss_reghorseshoe <- brm(
  bf(mvbind(Xmiss1, Xmiss2, Xmiss3, Xmiss4, Xmiss5) | mi() ~ (1||maternal_id)) + bf(ymiss | mi() ~ mi(Xmiss1) + mi(Xmiss2) + mi(Xmiss3) + mi(Xmiss4) + mi(Xmiss5) + (1||maternal_id)) + set_rescor(FALSE),
  prior = c(
    prior(gamma(1, 1), class = sd, resp = Xmiss1), 
    prior(gamma(1, 1), class = sd, resp = Xmiss2), 
    prior(gamma(1, 1), class = sd, resp = Xmiss3), 
    prior(gamma(1, 1), class = sd, resp = Xmiss4), 
    prior(gamma(1, 1), class = sd, resp = Xmiss5), 
    prior(gamma(1, 1), class = sd, resp = ymiss),
    prior(gamma(1, 1), class = sigma, resp = Xmiss1), 
    prior(gamma(1, 1), class = sigma, resp = Xmiss2), 
    prior(gamma(1, 1), class = sigma, resp = Xmiss3), 
    prior(gamma(1, 1), class = sigma, resp = Xmiss4), 
    prior(gamma(1, 1), class = sigma, resp = Xmiss5), 
    prior(gamma(1, 1), class = sigma, resp = ymiss),
    prior(horseshoe(df = 1, df_global = 1, scale_slab = 20, df_slab = 4, par_ratio = 2/3), class = b, resp = ymiss)
  ),
  data = dt,
  chains = 4, iter = 2000, warmup = 1000,
  file = 'project/fits/brmtest_mv_miss_reghorseshoe'
)
