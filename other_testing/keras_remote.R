# Run keras model remotely

library(brms)
library(keras)
library(data.table)

options(brms.file_refit = 'never')

abundance <- fread('project/data/16S_abundance.csv')
traits <- fread('project/data/seedling_data_updated.csv')

mv_traits <- c("days_to_germinate", "leaf_height_cm", "rooting_depth_cm", "root_tip_count", "no_primary_roots", "no_leaves")

traits_std <- traits[apply(traits[, lapply(.SD, \(x) !is.na(x)), .SDcols = mv_traits], 1, sum) > 1, 
                     mget(c('maternal_plant_code', 'country_origin', mv_traits))]

traits_std[, days_to_germinate := ifelse(is.na(days_to_germinate), median(days_to_germinate, na.rm = TRUE), days_to_germinate),
           by = .(maternal_plant_code)]

traits_std[, (mv_traits) := lapply(.SD, scale), .SDcols = mv_traits]

mv_trait_fit <- brm(
  mvbind(days_to_germinate, leaf_height_cm, rooting_depth_cm, root_tip_count, no_primary_roots, no_leaves) ~ country_origin + (1 | maternal_plant_code),
  data = traits_std, 
  chains = 4, iter = 2000, warmup = 1000, seed = 825,
  file = 'project/fits/mv_trait_fit'
)

pred_grid <- unique(traits_std[, .(maternal_plant_code, country_origin)])

trait_post_blups <- predict(mv_trait_fit, newdata = pred_grid, summary = FALSE)

abundance_long <- melt(abundance, id.vars = 'V1') 
abundance_long[, maternal_plant_code := substr(variable, 3, 3)]
abundance_long[, sampleID := substr(variable, 1, 2)]

abundance_wide <- dcast(abundance_long,  maternal_plant_code + sampleID ~ V1)

abundance_agg <- abundance_long[, .(value = sum(value)), by = .(V1, maternal_plant_code)]
abundance_agg_wide <- dcast(abundance_agg,  maternal_plant_code ~ V1)

abundance_agg_forfitting <- as.matrix(abundance_agg_wide[maternal_plant_code %in% pred_grid$maternal_plant_code, -1])

set.seed(1104)
n_samples <- 100
trait_post_blups_subsample <- trait_post_blups[sample(1:dim(trait_post_blups)[1], n_samples), , ]


model <- keras_model_sequential() |>
  layer_dense(units = 100, activation = 'relu', input_shape = dim(abundance_agg_forfitting)[2]) |>
  layer_dense(units = 32, activation = 'relu') |>
  layer_dense(units = dim(trait_post_blups)[3], activation = 'linear')

compile(model, loss = 'mse', optimizer = 'adam')

# For each of 100 MC samples, do leave-one-out cross-validation.
for (i in 1:n_samples) {
  for (j in 1:nrow(pred_grid)) {
    ab_train <- abundance_agg_forfitting[-j, ]
    ab_test <- abundance_agg_forfitting[j, , drop = FALSE]
    trait_train <- trait_post_blups_subsample[i, -j, ]
    trait_test <- trait_post_blups_subsample[i, j, , drop = FALSE]
    
    model |> fit(
      x = ab_train, y = trait_train, epochs = 100
    )
    
    pred_ij <- predict(model, ab_test)
    
  }
}

fit(model, abundance_agg_forfitting, tr1, epochs = 100)
