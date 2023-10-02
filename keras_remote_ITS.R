# Run keras model remotely

library(brms)
library(keras)
library(data.table)

options(brms.file_refit = 'on_change', brms.backend = 'cmdstanr', mc.cores = 4)

abundanceITS <- fread('project/data/ITS_Abundance.csv')
traits <- fread('project/data/seedling_data_updated.csv')

mv_traits <- c("days_to_germinate", "leaf_height_cm", "rooting_depth_cm", "root_tip_count", "no_primary_roots", "no_leaves")

traits_std <- traits[apply(traits[, lapply(.SD, \(x) !is.na(x)), .SDcols = mv_traits], 1, sum) > 1, 
                     mget(c('maternal_plant_code', 'country_origin', mv_traits))]

traits_std[, days_to_germinate := ifelse(is.na(days_to_germinate), median(days_to_germinate, na.rm = TRUE), days_to_germinate),
           by = .(maternal_plant_code)]

traits_std[, leaf_height_cm := leaf_height_cm + 1]
traits_std[, no_leaves := no_leaves + 1]
traits_std[, (mv_traits) := lapply(.SD, log), .SDcols = mv_traits]

traits_std[, (mv_traits) := lapply(.SD, scale), .SDcols = mv_traits]

mv_trait_fit <- brm(
  mvbind(days_to_germinate, leaf_height_cm, rooting_depth_cm, root_tip_count, no_primary_roots, no_leaves) ~ country_origin + (1 | maternal_plant_code),
  data = traits_std, 
  chains = 4, iter = 2000, warmup = 1000, seed = 825,
  file = 'project/fits/mv_trait_fit'
)

pred_grid <- unique(traits_std[, .(maternal_plant_code, country_origin)])

trait_post_blups <- predict(mv_trait_fit, newdata = pred_grid, summary = FALSE)

abundanceITS <- abundanceITS[rowSums(abundanceITS[,-1]) > 0]

abundanceITS_long <- melt(abundanceITS, id.vars = 'V1') 
abundanceITS_long[, maternal_plant_code := substr(variable, 3, 3)]
abundanceITS_long[, sampleID := substr(variable, 1, 2)]

abundanceITS_wide <- dcast(abundanceITS_long,  maternal_plant_code + sampleID ~ V1)

abundanceITS_agg <- abundanceITS_long[, .(value = sum(value)), by = .(V1, maternal_plant_code)]
abundanceITS_agg_wide <- dcast(abundanceITS_agg,  maternal_plant_code ~ V1)

abundanceITS_agg_forfitting <- as.matrix(abundanceITS_agg_wide[maternal_plant_code %in% pred_grid$maternal_plant_code, -1])

set.seed(1104)
n_samples <- 100
trait_post_blups_subsample <- trait_post_blups[sample(1:dim(trait_post_blups)[1], n_samples), , ]


# For each of 100 MC samples, do leave-one-out cross-validation.
# Model must be defined and compiled within each fold so that it does not save info from previous ones.
results <- list()
for (i in 1:n_samples) {
  for (j in 1:nrow(pred_grid)) {
    # Define and compile model.
    nnetmodel <- keras_model_sequential() |>
      layer_dense(units = 100, activation = 'relu', input_shape = dim(abundanceITS_agg_forfitting)[2]) |>
      layer_dropout(rate = 0.3) |>
      layer_dense(units = 32, activation = 'relu') |>
      layer_dropout(rate = 0.3) |>
      layer_dense(units = dim(trait_post_blups)[3], activation = 'linear')
    
    nnetmodel |> compile(loss = 'mse', optimizer = 'adam')
    
    # Create train and test split with just a single holdout row and the rest used for fitting.
    ab_train <- abundanceITS_agg_forfitting[-j, ]
    ab_test <- abundanceITS_agg_forfitting[j, , drop = FALSE]
    trait_train <- trait_post_blups_subsample[i, -j, ]
    trait_test <- t(as.matrix(trait_post_blups_subsample[i, j, ]))
    
    # Fit model and get predictions for the holdout row.
    nnetmodel |> fit(
      x = ab_train, y = trait_train, epochs = 1000, verbose = 0
    )
    
    results[[length(results) + 1]] <- data.frame(sample = i, row = j, predict(nnetmodel, ab_test))
    
  }
  message('Sample ', i, ' of ', n_samples, ' complete.')
}

results_df <- do.call(rbind, results)
fwrite(results_df, 'project/fits/nnet_mv_trait_results_ITS.csv')