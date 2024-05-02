// Single taxon by single trait model

data {
	int<lower=1> Nmicro;				                    // Number of offspring with microbiome data
	int<lower=1> Ntrait;				                    // Number of offspring with trait data
	int<lower=1> M;						                      // Number of maternal plants

	vector[Ntrait] y;			                          // Vector of standardized values of a trait
	vector[Nmicro] x;	                              // Vector of CLR transformed abundance of a single taxon (for now)
	vector<lower=0,upper=1>[M] pop;   		          // Population of origin of the maternal plants (0=Afghanistan, 1=Turkey)

	array[Nmicro] int<lower=1,upper=M> Mmicro;		  // Mapping to maternal plants 1-M for microbiome data
	array[Ntrait] int<lower=1,upper=M> Mtrait;		  // Mapping to maternal plants 1-M for trait data
}

parameters {
  real b0;                                        // Global intercept of trait across all maternal plants
  vector[3] b;                                    // Fixed effects at maternal plant level of population, taxon abundance, and their interaction
	vector[M] y_M;					                        // Predicted value for each maternal plant for the trait
	vector[M] x_M;                                  // Predicted value for each maternal plant for taxon abundance
	real<lower=0> sigma_yM;		                      // Variation (SD) among maternal plants in traits
	real<lower=0> sigma_xM;                         // Variation (SD) among maternal plants in taxon abundance
	real<lower=0> sigma;                            // Variation (SD) of model residuals
}

model {
  // Priors
  b0 ~ normal(0, 10);                             // Prior on intercept 
  b ~ normal(0, 10);                              // Fixed effects of population, taxon abundance, and their interactions
  sigma_yM ~ exponential(1);                      // SD of trait intercepts of each maternal plant
  sigma_xM ~ exponential(1);                      // SD of taxon abundance intercepts for each maternal plant
  sigma ~ exponential(1);                         // SD of model residuals
  
  // Likelihood
  y ~ normal(y_M[Mtrait], sigma_yM);		// Trait of each individual is drawn from a normal distribution with mean being the mean trait of its mother
  x ~ normal(x_M[Mmicro], sigma_xM);    // Taxon abundance of each individual is drawn from a normal distribution with mean being the mean abundance of its mother
  
  for (i in 1:M) {
    y_M[i] ~ normal(b0 + b[1] * pop[i] + b[2] * x_M[i] + b[3] * pop[i] * x_M[i], sigma);
  }

}

