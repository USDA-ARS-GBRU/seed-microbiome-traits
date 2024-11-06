# Make graph of model summaries from the missing and nonmissing models with 200 taxa.

library(brms)
library(tidyverse)

load('project/fits/brmtest_mv_summaries.RData')
#load('project/fits/brm_mv_summaries.RData')

# Combine coefficients into one data table and then parse the parameter name to get only the slopes
coef_summ <- bind_rows(
  tibble(model = 'mismatch', parameter = rownames(summ_miss$fixed), summ_miss$fixed),
  tibble(model = 'complete data', parameter = rownames(summ_nomiss$fixed), summ_nomiss$fixed)
) %>%
  mutate(type = if_else(grepl('Intercept', parameter), 'intercept', 'slope'),
         taxon = as.numeric(str_extract(parameter, '[0-9]+')))

pd = position_dodge(width = 0.5)
coef_summ %>%
  filter(type == 'slope', taxon %in% 1:10) %>%
  ggplot(aes(x = taxon, y = Estimate, ymin = `l-95% CI`, ymax = `u-95% CI`, group = model)) +
    geom_errorbar(width = 0, position = pd) +
    geom_point(aes(color = model), position = pd, size = 2) +
    coord_flip() + theme_bw() +
    scale_x_reverse(breaks = 1:10) +
    theme(legend.position = 'inside', legend.position.inside = c(0, 0), legend.justification = c(0, 0), legend.background = element_blank()) +
    see::scale_color_okabeito()

ggsave('project/examplefigtaxa1to10.png', height = 5, width = 6)
