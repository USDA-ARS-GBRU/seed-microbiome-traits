# Testing
tr1 <- trait_post_blups_subsample[1,,]

model <- keras_model_sequential() |>
  layer_dense(units = 100, activation = 'relu', input_shape = dim(abundance_agg_forfitting)[2]) |>
  layer_dense(units = 32, activation = 'relu') |>
  layer_dense(units = dim(trait_post_blups)[3], activation = 'linear')

compile(model, loss = 'mse', optimizer = 'adam')