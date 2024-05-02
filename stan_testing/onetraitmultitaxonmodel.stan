// Multivariate microbiome by single trait model
// QDR 2024-05-02

data {
	int<lower=1> Nmicro;				                    // Number of offspring with microbiome data
	int<lower=1> Ntrait;				                    // Number of offspring with trait data
	int<lower=1> M;						                      // Number of maternal plants
	int<lower=1> D;						                      // Number of taxa in the microbiome data

	vector[Ntrait] y;			                          // Vector of standardized values of a trait for each offspring with trait data
	array[Nmicro] vector[D] X;                      // Array of CLR transformed abundance vectors for each offspring with microbiome data
	vector<lower=0,upper=1>[M] pop;   		          // Population of origin of the maternal plants (0=Afghanistan, 1=Turkey)

	array[Nmicro] int<lower=1,upper=M> Mmicro;		  // Mapping to maternal plants 1-M for microbiome data
	array[Ntrait] int<lower=1,upper=M> Mtrait;		  // Mapping to maternal plants 1-M for trait data
}

parameters {
  real b0;                                        // Global intercept of trait across all maternal plants
  real b_pop;                                     // Fixed effects at maternal plant level of population, taxon abundance, and their interaction
  row_vector[D] b_tax;
  row_vector[D] b_tax_pop;
	vector[M] y_M;					                        // Predicted trait value for each maternal plant
	array[M] vector[D] X_M;                         // Array of predicted abundance vectors for each maternal plant
	real<lower=0> sigma_yM;		                      // Variation (SD) among maternal plants in traits
	corr_matrix[D] Omega_xM;                        // Correlation matrix of taxon abundances
  vector<lower=0>[D] tau;                         // Parameter to scale correlation matrix to variance-covariance matrix
	real<lower=0> sigma;                            // Variation (SD) of model residuals
}

transformed parameters {
  
}

model {
  // Priors
  b0 ~ normal(0, 10);                             // Prior on intercept 
  b_pop ~ normal(0, 10);                          // Fixed effects of population, taxon abundance, and their interactions
  b_tax ~ normal(0, 10);
  b_tax_pop ~ normal(0, 10);
  y_M ~ normal(0, 10);                            // Mean traits at maternal plant level
  for (i in 1:M) {
    X_M[i] ~ normal(0, 10);                       // Mean taxon abundances at maternal plant level
  }
  sigma_yM ~ exponential(1);                      // SD of trait means of maternal plants
  Omega_xM ~ lkj_corr(2);                         // Correlation matrix of taxon abundances for each maternal plant
  tau ~ cauchy(0, 2.5);                           // Scale correlation matrix to v-cv matrix
  sigma ~ exponential(1);                         // SD of model residuals
  
  // Likelihood
  // Trait of each offspring is drawn from a normal distribution with mean being the mean trait of its mother
  for (i in 1:Ntrait) {
    y[i] ~ normal(y_M[Mtrait[i]], sigma_yM);           
  }
  // Taxon abundance of each offspring is drawn from a multivariate normal distribution with mean being the mean abundance of its mother
  for (i in 1:Nmicro) {
    X[i] ~ multi_normal(X_M[Mmicro[i]], quad_form_diag(Omega_xM, tau));    
  }
  // Trait of each maternal plant is modeled as a function of its population of origin, its microbiome, and their interactions
  for (i in 1:M) {
    y_M[i] ~ normal(b0 + b_pop * pop[i] + b_tax * X_M[i] + b_tax_pop * (pop[i] .* X_M[i]), sigma);
  }

}

