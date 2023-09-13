aldexsamples <- aldex.clr(abundance[, -1], abundance_wide$maternal_plant_code, mc.samples = 1000)
aldexeff <- aldex.effect(aldexsamples, include.sample.summary = TRUE)

# Create model matrix from population and maternal plant

covariates <- abundance_wide[, .(maternal_plant_code)]
covariates[traits, on = .(maternal_plant_code), country := i.country_origin]
mm <- model.matrix(~ maternal_plant_code + country, covariates)

aldexsamples_mm <- aldex.clr(abundance[,-1], mm, mc.samples = 10, denom = 'all')
glm_test <- aldex.glm(aldexsamples_mm, mm)
glm_eff <- aldex.glm.effect(aldexsamples_mm)

###
conds2 <- c(rep("NS", 4), rep("S", 4), rep("thing",6))
x.all <- aldex(selex.sub, conds2, mc.samples=16, test="kw", effect=TRUE,
               include.sample.summary=TRUE, denom="all", verbose=FALSE, paired.test=FALSE)
