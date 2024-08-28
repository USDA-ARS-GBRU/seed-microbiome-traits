// Original model generated with brms 2.21.0
// Modified by QDR to constrain each of the X variables to have the same standard deviation
// Also if possible vectorize the code which currently has separate lines for each X variable
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
  
  array[N] vector[N_taxa] X;                   // X (microbiome array), which may contain missing values.
  
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
  array[Nmi_X] vector[N_taxa] Xmi;  // Array of estimated missing values for microbiome

  real Intercept_X;  // temporary intercept for centered predictors (CONSTRAINED TO BE SAME FOR ALL X)
  real<lower=0> sigma_X;  // dispersion parameter (CONSTRAINED TO BE SAME FOR ALL X)

  vector[Nmi_ymiss] Ymi;  // estimated missing values for trait
  real Intercept_Y;  // temporary intercept for centered predictors
  vector[Ksp_Y] zbsp_Y;  // unscaled coefficients
  
  // horseshoe shrinkage parameters
  real<lower=0> hs_global;  // global shrinkage parameter
  real<lower=0> hs_slab;  // slab regularization parameter
  vector<lower=0>[Kscales_Y] hs_local;  // local parameters for the horseshoe prior
  
  real<lower=0> sigma_Y;    // dispersion parameter for response variable y
  
  vector<lower=0>[1] sd_X;  // group-level standard deviations (CONSTRAINED TO BE SAME FOR ALL X)
    array[N_taxa] vector[N_grps] z; // Standardized group-level effects (random effects for microbiome array, X)

  vector<lower=0>[1] sd_Y;  // group-level standard deviations (Y)
  array[1] vector[N_grps] z_y;  // standardized group-level effects (random effects for trait, Y)
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
	r_X[taxon] = (sd_X[1] * (z[taxon][1]));
  }
  r_Y = (sd_X[1] * (z_y[1]));
  
  // Update lprior (CAN BE DONE AS A LOOP BECAUSE ALL X ARE THE SAME)
  for (taxon in 1:N_taxa) {
	lprior += student_t_lpdf(Intercept_X | 3, 0, 2.5); // FIXME change SD on this
	lprior += gamma_lpdf(sigma_X | 1, 1);
  }

  lprior += student_t_lpdf(Intercept_Y | 3, 0, 8.9); // FIXME change SD on this
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
	// FIXME this is as far as I got.
    vector[N] mu_Xmiss1 = rep_vector(0.0, N);
    vector[N] mu_Xmiss2 = rep_vector(0.0, N);
    vector[N] mu_Xmiss3 = rep_vector(0.0, N);
    vector[N] mu_Xmiss4 = rep_vector(0.0, N);
    vector[N] mu_Xmiss5 = rep_vector(0.0, N);
    vector[N] mu_ymiss = rep_vector(0.0, N);
    Yl_Xmiss1[Jmi_X] = Ymi_Xmiss1;
    Yl_Xmiss2[Jmi_X] = Ymi_Xmiss2;
    Yl_Xmiss3[Jmi_X] = Ymi_Xmiss3;
    Yl_Xmiss4[Jmi_X] = Ymi_Xmiss4;
    Yl_Xmiss5[Jmi_X] = Ymi_Xmiss5;
    Yl_ymiss[Jmi_y] = Ymi_ymiss;
    mu_Xmiss1 += Intercept_Xmiss;
    mu_Xmiss2 += Intercept_Xmiss;
    mu_Xmiss3 += Intercept_Xmiss;
    mu_Xmiss4 += Intercept_Xmiss;
    mu_Xmiss5 += Intercept_Xmiss;
    mu_ymiss += Intercept_ymiss;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu_Xmiss1[n] += r_1_Xmiss1_1[J[n]] * Z_1_Xmiss1_1[n];
	  mu_Xmiss2[n] += r_2_Xmiss2_1[J[n]] * Z_2_Xmiss2_1[n];
	  mu_Xmiss3[n] += r_3_Xmiss3_1[J[n]] * Z_3_Xmiss3_1[n];
	  mu_Xmiss4[n] += r_4_Xmiss4_1[J[n]] * Z_4_Xmiss4_1[n];
	  mu_Xmiss5[n] += r_5_Xmiss5_1[J[n]] * Z_5_Xmiss5_1[n];

      mu_ymiss[n] += (bsp_ymiss[1]) * Yl_Xmiss1[n] + (bsp_ymiss[2]) * Yl_Xmiss2[n] + (bsp_ymiss[3]) * Yl_Xmiss3[n] + (bsp_ymiss[4]) * Yl_Xmiss4[n] + (bsp_ymiss[5]) * Yl_Xmiss5[n] + r_6_ymiss_1[J_6_ymiss[n]] * Z_6_ymiss_1[n];
    }
    target += normal_lpdf(Yl_Xmiss1 | mu_Xmiss1, sigma_Xmiss);
    target += normal_lpdf(Yl_Xmiss2 | mu_Xmiss2, sigma_Xmiss);
    target += normal_lpdf(Yl_Xmiss3 | mu_Xmiss3, sigma_Xmiss);
    target += normal_lpdf(Yl_Xmiss4 | mu_Xmiss4, sigma_Xmiss);
    target += normal_lpdf(Yl_Xmiss5 | mu_Xmiss5, sigma_Xmiss);
    target += normal_lpdf(Yl_ymiss | mu_ymiss, sigma_ymiss);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(zbsp_ymiss);
  target += student_t_lpdf(hs_local_ymiss | hs_df_ymiss, 0, 1)
    - rows(hs_local_ymiss) * log(0.5);
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(z_2[1]);
  target += std_normal_lpdf(z_3[1]);
  target += std_normal_lpdf(z_4[1]);
  target += std_normal_lpdf(z_5[1]);
  target += std_normal_lpdf(z_y[1]);
}
