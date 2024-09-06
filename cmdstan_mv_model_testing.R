library(mvtnorm)
library(cmdstanr)

# Increase number of taxa
n_mothers <- 20
n_taxa <- 200
offspring_per_mother <- 10 # 5 will be retained for traits, 5 for microbiome

set.seed(1)

X_maternal <- rmvnorm(n_mothers, mean = rep(0, n_taxa), sigma = diag(n_taxa))
sigma_maternal <- cov(X_maternal)

# Coefficients indicating which taxa predict the outcome.
# We will not include any interaction effect.
# Only include a few taxa with a nonzero effect.
beta <- c(50, 10, -50, rep(0, n_taxa-3))

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