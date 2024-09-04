// Original model generated with brms 2.21.0
// Modified by QDR to constrain each of the X variables to have the same standard deviation and the same intercept
// Where brms had separate lines for each X variable, convert to arrays and loops so the code is scalable to any number of taxa
functions {
  /* Efficient computation of the horseshoe scale parameters
   * see Appendix C.1 in https://projecteuclid.org/euclid.ejs/1513306866
   * Args:
   *   lambda: local shrinkage parameters
   *   tau: global shrinkage parameter
   *   c2: slab regularization parameter
   * Returns:
   *   scale parameter vector of the horseshoe prior
   */
  vector scales_horseshoe(vector lambda, real tau, real c2) {
    int K = rows(lambda);
    vector[K] lambda2 = square(lambda);
    vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
    return lambda_tilde * tau;
  }
  /* compute scale parameters of the R2D2 prior
   * Args:
   *   phi: local weight parameters
   *   tau2: global scale parameter
   * Returns:
   *   scale parameter vector of the R2D2 prior
   */
  vector scales_R2D2(vector phi, real tau2) {
    return sqrt(phi * tau2);
  }

}
data {
  int<lower=1> N;                       // total number of observations
  int<lower=1> N_taxa;                  // Number of taxa (predictors)
  
  int<lower=0> Nmi_X;                   // number of missing X (CONSTRAINED TO BE THE SAME FOR ALL X)
  array[Nmi_X] int<lower=1> Jmi_X;      // positions of missings (CONSTRAINED TO BE THE SAME FOR ALL X)
  
  array[N_taxa] vector[N] X;                   // X (microbiome array), which may contain missing values.
  
  vector[N] Y;                          // Y (trait), response variable
  int<lower=0> Nmi_Y;                   // number of missing Y
  array[Nmi_Y] int<lower=1> Jmi_Y;      // positions of missing Y
  
  int<lower=1> Ksp_Y;               // number of special effects terms
  int<lower=1> Kscales_Y;           // number of local scale parameters
  
  // data for the horseshoe prior
  real<lower=0> hs_df;            // local degrees of freedom
  real<lower=0> hs_df_global;     // global degrees of freedom
  real<lower=0> hs_df_slab;       // slab degrees of freedom
  real<lower=0> hs_scale_global;  // global prior scale
  real<lower=0> hs_scale_slab;    // slab prior scale
  
  real<lower=0> prior_sd_Intercept_X;    // SD of Student t prior on X intercept
  real<lower=0> prior_sd_Intercept_Y;    // SD of Student t prior on Y intercept
  
  int<lower=1> N_grps;            // number of grouping levels (CONSTRAINED TO BE THE SAME FOR ALL X AND Y)
  array[N] int<lower=1> J;        // grouping indicator per observation (CONSTRAINED TO BE THE SAME FOR ALL X AND Y)
  
  // group-level predictor values for each X
  array[N_taxa] vector[N] Z_X;
  
  // data for group-level effects of y
  // group-level predictor values
  vector[N] Z_Y;
  
  int prior_only;  // should the likelihood be ignored?
}

transformed data {
}

parameters {
  array[N_taxa] vector[Nmi_X] Xmi;  // Array of estimated missing values for microbiome

  real Intercept_X;  // temporary intercept for centered predictors (CONSTRAINED TO BE SAME FOR ALL X)
  real<lower=0> sigma_X;  // dispersion parameter (CONSTRAINED TO BE SAME FOR ALL X)

  vector[Nmi_Y] Ymi;  // estimated missing values for trait
  real Intercept_Y;  // temporary intercept for centered predictors
  vector[Ksp_Y] zbsp_Y;  // unscaled coefficients
  
  // horseshoe shrinkage parameters
  real<lower=0> hs_global;  // global shrinkage parameter
  real<lower=0> hs_slab;  // slab regularization parameter
  vector<lower=0>[Kscales_Y] hs_local;  // local parameters for the horseshoe prior
  
  real<lower=0> sigma_Y;    // dispersion parameter for response variable y
  
  real<lower=0> sd_X;  // group-level standard deviations (CONSTRAINED TO BE SAME FOR ALL X)
  array[N_taxa] vector[N_grps] z; // Standardized group-level effects (random effects for microbiome array, X)

  vector<lower=0>[1] sd_Y;  // group-level standard deviations (Y)
  vector[N_grps] z_y;  // standardized group-level effects (random effects for trait, Y)
}

transformed parameters {
  vector[Ksp_Y] bsp_Y;  // scaled coefficients
  vector<lower=0>[Ksp_Y] sdbsp_Y;  // SDs of the coefficients
  vector<lower=0>[Kscales_Y] scales_Y;  // local horseshoe scale parameters
  
  array[N_taxa] vector[N_grps] r_X; // actual group-level effects for microbiome array (X)
  vector[N_grps] r_Y;  // actual group-level effects for trait (Y)
  
  real lprior = 0;  // prior contributions to the log posterior
  
  // compute horseshoe scale parameters
  scales_Y = scales_horseshoe(hs_local, hs_global, hs_scale_slab^2 * hs_slab);
  sdbsp_Y = scales_Y[(1):(Ksp_Y)];
  bsp_Y = zbsp_Y .* sdbsp_Y;  // scale coefficients
  
  for (taxon in 1:N_taxa) {
	r_X[taxon] = sd_X .* z[taxon];
  }
  r_Y = sd_X .* z_y;
  
  // Update lprior (CAN BE DONE AS A LOOP BECAUSE ALL X ARE THE SAME)
  for (taxon in 1:N_taxa) {
	lprior += student_t_lpdf(Intercept_X | 3, 0, prior_sd_Intercept_X);
	lprior += gamma_lpdf(sigma_X | 1, 1);
  }

  lprior += student_t_lpdf(Intercept_Y | 3, 0, prior_sd_Intercept_Y);
  lprior += student_t_lpdf(hs_global | hs_df_global, 0, hs_scale_global * sigma_Y)
    - 1 * log(0.5);
  lprior += inv_gamma_lpdf(hs_slab | 0.5 * hs_df_slab, 0.5 * hs_df_slab);
  lprior += gamma_lpdf(sigma_Y | 1, 1);
  
  for (taxon in 1:N_taxa) {
	lprior += gamma_lpdf(sd_X | 1, 1); // All taxa constrained to have same SD
  }

  lprior += gamma_lpdf(sd_Y | 1, 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // vectors combining observed and missing responses
	array[N_taxa] vector[N] Xl = X;
    vector[N] Yl = Y;
	
    // initialize linear predictor terms
	array[N_taxa] vector[N] mu_X;
	vector[N] mu_Y = rep_vector(0.0, N);
	
	for (taxon in 1:N_taxa) {
		mu_X[taxon] = rep_vector(0.0, N);
		Xl[taxon][Jmi_X] = Xmi[taxon];
		mu_X[taxon] += Intercept_X;
	}	
	
    Yl[Jmi_Y] = Ymi;
    mu_Y += Intercept_Y;
	
    for (n in 1:N) {
      // add more terms to the linear predictor
	  for (taxon in 1:N_taxa) {
		mu_X[taxon][n] += r_X[taxon][J[n]] * Z_X[taxon][n];
		mu_Y[n] += (bsp_Y[taxon]) * Yl[n];
	  }

      mu_Y[n] += r_Y[J[n]] * Z_Y[n];
    }
	
	for (taxon in 1:N_taxa) {
		target += normal_lpdf(Xl[taxon] | mu_X[taxon], sigma_X);
	}
	
    target += normal_lpdf(Yl | mu_Y, sigma_Y);
  }
  
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(zbsp_Y);
  target += student_t_lpdf(hs_local | hs_df, 0, 1) - rows(hs_local) * log(0.5);
  
  for (taxon in 1:N_taxa) {
	target += std_normal_lpdf(Z_X[taxon]);
  }
  
  target += std_normal_lpdf(Z_Y);
}
