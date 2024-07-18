# Go back to simple form. We have a relationship between x and y, and some offspring from each mother with those means
# Generate data for all offspring, but delete the x data for half and the y data for the other half

set.seed(201)
n_mothers <- 10
maternal_mean_x <- rnorm(n_mothers, mean = 10, sd = 2)
maternal_mean_y <- 5 + 2 * maternal_mean_x + rnorm(n_mothers, mean = 0, sd = 1)

# 5 offspring are taken from each mother for x, and 5 for y. We will model the within-mother variability as being the same for each.
offspring_per_mother <- 10
offspring_x <- rnorm(offspring_per_mother * n_mothers, mean = rep(maternal_mean_x, offspring_per_mother), sd = 2)
offspring_y <- rnorm(offspring_per_mother * n_mothers, mean = rep(maternal_mean_y, offspring_per_mother), sd = 1)

dt <- data.frame(
  maternal_id = factor(rep(1:n_mothers, each = offspring_per_mother)),
  offspring_id = 1:offspring_per_mother,
  x = offspring_x,
  y = offspring_y
)

# Within each mother, set half of the values to be missing for x, and the other half for y.
dt$xmiss <- ifelse(dt$offspring_id %in% 1:(offspring_per_mother/2), dt$x, NA)
dt$ymiss <- ifelse(dt$offspring_id %in% 1:(offspring_per_mother/2), NA, dt$y)

library(brms)

options(mc.cores = 4, brms.backend = 'cmdstanr', brms.file_refit = 'on_change')

get_prior(
  bf(x | mi() ~ (1|maternal_id)) + bf(y | mi() ~ mi(x) + (1|maternal_id)) + set_rescor(FALSE),
  data = dt
)

## Model fit with no missing data
# We should recover xint 9, yint 15, slope ~0.9
modrescorfalse_nomiss <- brm(
  bf(x ~ (1|maternal_id)) + bf(y ~ x + (1|maternal_id)) + set_rescor(FALSE),
  prior = c(
    prior(gamma(1, 1), class = sd, resp = x),
    prior(gamma(1, 1), class = sd, resp = y),
    prior(gamma(1, 1), class = sigma, resp = x),
    prior(gamma(1, 1), class = sigma, resp = y),
    prior(normal(0, 5), class = b, resp = y)
  ),
  data = dt,
  chains = 4, iter = 5000, warmup = 2500,
  control = list(adapt_delta = 0.95),
  file = 'project/fits/brmtest_rescorfalse_nomiss'
)

# This is getting close. The regression coefficients and sigma y are mixing poorly. Not sure how to address.
# Maybe because there's, like, no data
modrescorfalse <- brm(
  bf(xmiss | mi() ~ (1|maternal_id)) + bf(ymiss | mi() ~ mi(xmiss) + (1|maternal_id)) + set_rescor(FALSE),
  prior = c(
    prior(gamma(1, 1), class = sd, resp = xmiss),
    prior(gamma(1, 1), class = sd, resp = ymiss),
    prior(gamma(1, 1), class = sigma, resp = xmiss),
    prior(gamma(1, 1), class = sigma, resp = ymiss),
    prior(normal(0, 5), class = b, resp = ymiss)
  ),
  data = dt,
  chains = 4, iter = 5000, warmup = 2500,
  control = list(adapt_delta = 0.95),
  seed = 1506,
  file = 'project/fits/brmtest_rescorfalse'
)

modrescortrue <- brm(
  bf(xmiss | mi() ~ (1|maternal_id)) + bf(ymiss | mi() ~ mi(xmiss) + (1|maternal_id)) + set_rescor(TRUE),
  prior = c(
    prior(gamma(1, 1), class = sd, resp = xmiss),
    prior(gamma(1, 1), class = sd, resp = ymiss),
    prior(gamma(1, 1), class = sigma, resp = xmiss),
    prior(gamma(1, 1), class = sigma, resp = ymiss),
    prior(normal(0, 5), class = b, resp = ymiss)
  ),
  data = dt,
  chains = 4, iter = 5000, warmup = 2500,
  control = list(adapt_delta = 0.95),
  seed = 1522,
  file = 'project/fits/brmtest_rescortrue'
)
