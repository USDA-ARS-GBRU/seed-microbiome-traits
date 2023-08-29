// Try to implement multivariate probit regression for the microbiome part.
// Use https://mc-stan.org/docs/2_22/stan-users-guide/multivariate-outcomes.html
// Newer: https://mc-stan.org/docs/stan-users-guide/multivariate-outcomes.html


functions {
  int sum2d(array[,] int a) {
    int s = 0;
    for (i in 1:size(a)) {
      s += sum(a[i]);
    }
    return s;
  }
}

data {
	int<lower=1> Nmicro;				// Number of offspring with microbiome data
	int<lower=1> Ntrait;				// Number of offspring with trait data
	int<lower=1> M;						  // Number of maternal plants
	int<lower=1> D;						  // Number of taxa in the microbiome data
	int<lower=1> K;             // 
	
	vector<lower=0>[Ntrait] y;			          // So far we are using a single trait as response variable. Later it'll be made more complex. 
	array[Nmicro, D] int<lower=0, upper=1> X;	// Matrix of presence-absences of the different taxa
	int<lower=0, upper=1> pop[M];		    // Population of origin of the maternal plants (0=Afghanistan, 1=Turkey)

	int<lower=1,upper=M> Mmicro[Nmicro];		  // Mapping to maternal plants 1-M for microbiome data
	int<lower=1,upper=M> Mtrait[Ntrait];		  // Mapping to maternal plants 1-M for trait data
}

transformed data {
	vector[Ntrait] logy;
	
	array[sum2d(X)] int<lower=1, upper=Nmicro> n_pos;
  array[size(n_pos)] int<lower=1, upper=D> d_pos;
  int<lower=0> N_neg;
  array[(Nmicro * D) - size(n_pos)] int<lower=1, upper=Nmicro> n_neg;
  array[size(n_neg)] int<lower=1, upper=D> d_neg;
  
  logy = log(y);
  
  N_pos = size(n_pos);
  N_neg = size(n_neg);
  {
    int i;
    int j;
    i = 1;
    j = 1;
    for (n in 1:Nmicro) {
      for (d in 1:D) {
        if (X[n,d] == 1) {
          n_pos[i] = n;
          d_pos[i] = d;
          i += 1;
        } else {
          n_neg[j] = n;
          d_neg[j] = d;
          j += 1;
        }
      }
    }
  }
}

parameters {
	real alpha_trait;					      // Global mean of the trait
	real beta_trait;					      // Effect of population on the trait
	vector[M] u_trait;					    // Intercept for each maternal plant for the trait
	real<lower=0> sigma_alpha;			// Variation in intercepts

  // Parameters for multivariate probit:
  matrix[D, K] beta;
  cholesky_factor_corr[D] L_Omega;
  vector<lower=0>[N_pos] z_pos;
  vector<upper=0>[N_neg] z_neg;
	
}

transformed parameters {
  vector[D] z[N];
  for (n in 1:N_pos)
    z[n_pos[n], d_pos[n]] = z_pos[n];
  for (n in 1:N_neg)
    z[n_neg[n], d_neg[n]] = z_neg[n];
}

model {
  // Priors
  b0 ~ normal(0, 10);
  b1 ~ normal(0, 10);
  u ~ normal(0, sigma_u);
  sigma_u ~ exponential(1);                           // SD of random effects
  sigma ~ exponential(1);                             // SD of model residuals
  
  // Likelihood
  logy ~ normal(b0 + b1 * pop + u[Mtrait], sigma);		// Trait of each individual is drawn from a normal distribution
  
  L_Omega ~ lkj_corr_cholesky(4);
  to_vector(beta) ~ normal(0, 5);
  {
    vector[D] beta_x[Nmicro];
    for (n in 1:Nmicro)
    beta_x[n] = beta * x[n];
    z ~ multi_normal_cholesky(beta_x, L_Omega);
  }
  
  // Next, here we need to put the likelihood for the part of the model where we use the params of the multivariate probit to predict the maternal plant level predicted values of the trait.
  
}

generated quantities {
  corr_matrix[D] Omega;
  Omega = multiply_lower_tri_self_transpose(L_Omega);
}
