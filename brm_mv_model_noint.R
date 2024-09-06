# Try to make model, within brms, have a shared intercept 
modmv_miss_reghorseshoe_noint <- brm(
  bf(mvbind(Xmiss1, Xmiss2, Xmiss3, Xmiss4, Xmiss5) | mi() ~ 0 + (1||maternal_id)) + bf(ymiss | mi() ~ mi(Xmiss1) + mi(Xmiss2) + mi(Xmiss3) + mi(Xmiss4) + mi(Xmiss5) + (1||maternal_id)) + set_rescor(FALSE),
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
  file = 'project/fits/brmtest_mv_miss_reghorseshoe_noint'
)
