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
  int<lower=1> N;  // total number of observations
  int<lower=1> N_taxa; // Number of taxa (predictors)
  
  int<lower=0> Nmi_X;  // number of missing X (CONSTRAINED TO BE THE SAME FOR ALL X)
  array[Nmi_X] int<lower=1> Jmi_X;  // positions of missings (CONSTRAINED TO BE THE SAME FOR ALL X)
  
  vector[N] Y_Xmiss1;  // response variable
  vector[N] Y_Xmiss2;  // response variable
  vector[N] Y_Xmiss3;  // response variable
  vector[N] Y_Xmiss4;  // response variable
  vector[N] Y_Xmiss5;  // response variable
  
  vector[N] Y_ymiss;  // response variable
  int<lower=0> Nmi_ymiss;  // number of missings
  array[Nmi_ymiss] int<lower=1> Jmi_y;  // positions of missing Y
  int<lower=1> Ksp_ymiss;  // number of special effects terms
  int<lower=1> Kscales_ymiss;  // number of local scale parameters
  
  // data for the horseshoe prior
  real<lower=0> hs_df_ymiss;  // local degrees of freedom
  real<lower=0> hs_df_global_ymiss;  // global degrees of freedom
  real<lower=0> hs_df_slab_ymiss;  // slab degrees of freedom
  real<lower=0> hs_scale_global_ymiss;  // global prior scale
  real<lower=0> hs_scale_slab_ymiss;  // slab prior scale
  
  int<lower=1> N_grps;  // number of grouping levels (CONSTRAINED TO BE THE SAME FOR ALL X AND Y)
  array[N] int<lower=1> J;  // grouping indicator per observation (CONSTRAINED TO BE THE SAME FOR ALL X AND Y)
  
  // group-level predictor values for each X
  vector[N] Z_1_Xmiss1_1;
  vector[N] Z_2_Xmiss2_1;
  vector[N] Z_3_Xmiss3_1;
  vector[N] Z_4_Xmiss4_1;
  vector[N] Z_5_Xmiss5_1;
  
  // data for group-level effects of y
  // group-level predictor values
  vector[N] Z_6_ymiss_1;
  
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  vector[Nmi_X] Ymi_Xmiss1;  // estimated missings
  vector[Nmi_X] Ymi_Xmiss2;  // estimated missings
  vector[Nmi_X] Ymi_Xmiss3;  // estimated missings
  vector[Nmi_X] Ymi_Xmiss4;  // estimated missings
  vector[Nmi_X] Ymi_Xmiss5;  // estimated missings
  
  real Intercept_Xmiss;  // temporary intercept for centered predictors (CONSTRAINED TO BE SAME FOR ALL X)
  real<lower=0> sigma_Xmiss;  // dispersion parameter (CONSTRAINED TO BE SAME FOR ALL X)

  vector[Nmi_ymiss] Ymi_ymiss;  // estimated missings
  real Intercept_ymiss;  // temporary intercept for centered predictors
  vector[Ksp_ymiss] zbsp_ymiss;  // unscaled coefficients
  
  // horseshoe shrinkage parameters
  real<lower=0> hs_global_ymiss;  // global shrinkage parameter
  real<lower=0> hs_slab_ymiss;  // slab regularization parameter
  vector<lower=0>[Kscales_ymiss] hs_local_ymiss;  // local parameters for the horseshoe prior
  
  real<lower=0> sigma_ymiss;  // dispersion parameter
  vector<lower=0>[1] sd_x;  // group-level standard deviations (CONSTRAINED TO BE SAME FOR ALL X)
  
  array[1] vector[N_grps] z_1;  // standardized group-level effects
  array[1] vector[N_grps] z_2;  // standardized group-level effects
  array[1] vector[N_grps] z_3;  // standardized group-level effects
  array[1] vector[N_grps] z_4;  // standardized group-level effects
  array[1] vector[N_grps] z_5;  // standardized group-level effects
  vector<lower=0>[1] sd_y;  // group-level standard deviations (Y)
  array[1] vector[N_grps] z_y;  // standardized group-level effects
}
transformed parameters {
  vector[Ksp_ymiss] bsp_ymiss;  // scaled coefficients
  vector<lower=0>[Ksp_ymiss] sdbsp_ymiss;  // SDs of the coefficients
  vector<lower=0>[Kscales_ymiss] scales_ymiss;  // local horseshoe scale parameters
  vector[N_grps] r_1_Xmiss1_1;  // actual group-level effects
  vector[N_grps] r_2_Xmiss2_1;  // actual group-level effects
  vector[N_grps] r_3_Xmiss3_1;  // actual group-level effects
  vector[N_grps] r_4_Xmiss4_1;  // actual group-level effects
  vector[N_grps] r_5_Xmiss5_1;  // actual group-level effects
  vector[N_grps] r_6_ymiss_1;  // actual group-level effects
  real lprior = 0;  // prior contributions to the log posterior
  // compute horseshoe scale parameters
  scales_ymiss = scales_horseshoe(hs_local_ymiss, hs_global_ymiss, hs_scale_slab_ymiss^2 * hs_slab_ymiss);
  sdbsp_ymiss = scales_ymiss[(1):(Ksp_ymiss)];
  bsp_ymiss = zbsp_ymiss .* sdbsp_ymiss;  // scale coefficients
  r_1_Xmiss1_1 = (sd_x[1] * (z_1[1]));
  r_2_Xmiss2_1 = (sd_x[1] * (z_2[1]));
  r_3_Xmiss3_1 = (sd_x[1] * (z_3[1]));
  r_4_Xmiss4_1 = (sd_x[1] * (z_4[1]));
  r_5_Xmiss5_1 = (sd_x[1] * (z_5[1]));
  r_6_ymiss_1 = (sd_y[1] * (z_y[1]));
  
  // Update lprior (CAN BE DONE AS A LOOP BECAUSE ALL X ARE THE SAME)
  for (taxon in 1:N_taxa) {
	lprior += student_t_lpdf(Intercept_Xmiss | 3, 0, 2.5); // FIXME change SD on this
	lprior += gamma_lpdf(sigma_Xmiss | 1, 1);
  }

  lprior += student_t_lpdf(Intercept_ymiss | 3, 0, 8.9); // FIXME change SD on this
  lprior += student_t_lpdf(hs_global_ymiss | hs_df_global_ymiss, 0, hs_scale_global_ymiss * sigma_ymiss)
    - 1 * log(0.5);
  lprior += inv_gamma_lpdf(hs_slab_ymiss | 0.5 * hs_df_slab_ymiss, 0.5 * hs_df_slab_ymiss);
  lprior += gamma_lpdf(sigma_ymiss | 1, 1);
  
  for (taxon in 1:N_taxa) {
	lprior += gamma_lpdf(sd_x | 1, 1); // All taxa constrained to have same SD
  }

  lprior += gamma_lpdf(sd_y | 1, 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // vectors combining observed and missing responses
    vector[N] Yl_Xmiss1 = Y_Xmiss1;
    vector[N] Yl_Xmiss2 = Y_Xmiss2;
    vector[N] Yl_Xmiss3 = Y_Xmiss3;
    vector[N] Yl_Xmiss4 = Y_Xmiss4;
    vector[N] Yl_Xmiss5 = Y_Xmiss5;
    vector[N] Yl_ymiss = Y_ymiss;
	
    // initialize linear predictor terms
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
