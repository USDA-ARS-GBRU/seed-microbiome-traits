# Dirichlet regression

abmat <- abundance_wide[,-c(1:2)]

colsumsab <- colSums(abmat > 0)
bigones <- tail(sort(colsumsab))
abmatreduced <- abmat[, mget(names(bigones))]
propmatreduced <- sweep(abmatreduced, 1, rowSums(abmatreduced), `/`)

yuse <- as.matrix(propmatreduced[is.finite(rowSums(propmatreduced)),])

datfit <- as.data.frame(abundance_wide[is.finite(rowSums(propmatreduced)), 1:2])
datfit$y <- yuse

fit <- brm(y ~ maternal_plant_code, data = datfit, family = dirichlet)
