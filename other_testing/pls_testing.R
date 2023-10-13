pls_16s_test <- pls(X = abundance16S_agg_forfitting, Y = trait_post_blups[101,,], ncomp = 3, mode = 'regression')
pls_16s_testpred_insample <- predict(pls_16s_test, newdata = abundance16S_agg_forfitting)

sqrt(apply((trait_post_blups[101,,] - pls_16s_testpred_insample$predict[,,1])^2, 2, mean))
sqrt(apply((trait_post_blups[101,,] - pls_16s_testpred_insample$predict[,,2])^2, 2, mean))
sqrt(apply((trait_post_blups[101,,] - pls_16s_testpred_insample$predict[,,3])^2, 2, mean))

# Example with LOO cross validation
pls_16s_loo_test <- lapply(1:9, \(i) pls(X = abundance16S_agg_forfitting[-i,], Y = trait_post_blups[101,-i,], ncomp = 3, mode = 'regression'))
pls_16s_loo_predoutsample <- lapply(1:9, \(i) predict(pls_16s_loo_test[[i]], newdata = abundance16S_agg_forfitting[i, , drop = FALSE]))

# Reshape out of sample prediction into three matrices
pls_16s_loo_predmats <- list()
for (i in 1:3) {
  pls_16s_loo_predmats[[i]] <- do.call(rbind, lapply(pls_16s_loo_predoutsample, \(p) p$predict[,,i]))
}

pls_16_loo_rmses <- lapply(1:3, \(i) sqrt(apply((trait_post_blups[101,,] - pls_16s_loo_predmats[[i]])^2, 2, mean)))