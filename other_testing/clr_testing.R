# CLR transformation of abundance16S and ITS. logratio(x+1) transformation is used.
# FIXME first remove maternal plant J, then remove all taxa that are all zero, then do the transformation.

aggregate_abundance <- function(abundance) {
  # Remove maternal plant J as it has no seedling traits.
  J_columns <- grep('J$', names(abundance), value = TRUE)
  abundance[, c(J_columns) := NULL]
  
  # Remove any taxa that now have all zero abundance in the remaining 9 maternal plants.
  zero_rows <- rowSums(abundance[,-1]) == 0
  abundance <- abundance[!zero_rows]
  
  # Do CLR transformation, adding 1 to all counts. Replace original values with transformed.
  abundance_CLR <- logratio.transfo(abundance[,-1], logratio = 'CLR', offset = 1)
  abundance_CLR <- cbind(abundance[, .(V1)], as(abundance_CLR, 'matrix'))
  
  # Reshape to longform
  abundance_long <- melt(abundance_CLR, id.vars = 'V1') 
  abundance_long[, maternal_plant_code := substr(variable, 3, 3)]
  abundance_long[, sampleID := substr(variable, 1, 2)]
  
  # Aggregate by taking the average log ratio of each taxon grouped by maternal plant
  abundance_agg <- abundance_long[, .(value = mean(value)), by = .(V1, maternal_plant_code)]
  
  # Return to wideform
  dcast(abundance_agg, maternal_plant_code ~ V1)
}

abundance16S_CLR <- logratio.transfo(abundance16S[,-1], logratio = 'CLR', offset = 1)
abundance16S_CLR <- cbind(abundance16S[, .(V1)], as(abundance16S_CLR, 'matrix'))

abundance16S_long <- melt(abundance16S_CLR, id.vars = 'V1') 
abundance16S_long[, maternal_plant_code := substr(variable, 3, 3)]
abundance16S_long[, sampleID := substr(variable, 1, 2)]

abundance16S_wide <- dcast(abundance16S_long,  maternal_plant_code + sampleID ~ V1)

abundance16S_agg <- abundance16S_long[, .(value = mean(value)), by = .(V1, maternal_plant_code)]
abundance16S_agg_wide <- dcast(abundance16S_agg,  maternal_plant_code ~ V1)


