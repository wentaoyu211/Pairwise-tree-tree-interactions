data{
  int ND; // number of total data points
  int NS; // number of species
  int NP; // number of plots in total
  int NY; // number of years
  int N_spy; // number of unique species * year
  int N_spp; // number of unique species * plot
  real focal_bm[ND];  // biomass of focal tree
  real next_bm[ND];   // the next year biomass of corresponding focal tree
  int fsp_ID[ND]; // identify species of each data point
  int year_ID[ND];   // identify the collecting year for each data point
  int spy[ND];    // species * year indices
  int fplot_ID[ND]; // identify plot for each data point
  int spp[ND];      // species * plot indices
}

parameters{
  real<lower=0> beta[NS]; // species specific growth parameter
  real theta;             // test if the growth rate scale with power law (0.75)  
  vector[NP] eps_p_raw;   // model plots as random effect
  real<lower=0> sigma_eps_p;
  vector[N_spp] eps_spp_raw;   // model species * plots
  real<lower=0> sigma_eps_spp;   
  vector[NY] eps_y_raw;   // model year as random effect
  real<lower=0> sigma_eps_y;
  vector[N_spy] eps_spy_raw;   // model species * year 
  real<lower=0> sigma_eps_spy;
  real<lower=0> sigma;    // sigma is observation error   
}

transformed parameters{
  vector[NP] eps_p;
  vector[NY] eps_y; 
  vector[N_spy] eps_spy;
  vector[N_spp] eps_spp;
  
  //non-centred parameterization to improve efficiency and convergence of hierachical model
  eps_p = sigma_eps_p * eps_p_raw;
  eps_y = sigma_eps_y * eps_y_raw;
  eps_spy = sigma_eps_spy * eps_spy_raw;
  eps_spp = sigma_eps_spp * eps_spp_raw;
}

model{
  // for model simplification
  real predBM;
  
  // specify priors
  for (i in 1:NS){
    beta[i] ~ exponential(1.0);
  	  // based on 3/4 exponent
	}
  theta ~ normal(0.75, 1);
  
  // random effect structure
  sigma_eps_p ~ normal(0, 0.1);
  sigma_eps_y ~ normal(0, 0.1);
  sigma_eps_spy ~ normal(0, 0.1);  
  sigma_eps_spp ~ normal(0, 0.1); 
  eps_p_raw ~ std_normal();
  eps_y_raw ~ std_normal();
  eps_spy_raw ~ std_normal();
  eps_spp_raw ~ std_normal();
   
  for(i in 1:ND){    // loop over all data points
        // prediction
        predBM = focal_bm[i] + beta[fsp_ID[i]]*focal_bm[i]^theta * max({1 + eps_p[fplot_ID[i]] + eps_spp[spp[i]] + eps_y[year_ID[i]] + eps_spy[spy[i]], 0}); 
		
        // likelihood		
	next_bm[i] ~ lognormal( log(predBM), sigma );
  }
}

generated quantities {
  real predBM; 
  vector[ND] log_lik;  // number of log likelihood for loo
  
  for (i in 1:ND){    // loop through focal individual tree
      // prediction
      predBM = focal_bm[i] + beta[fsp_ID[i]]*focal_bm[i]^theta * max({1 + eps_p[fplot_ID[i]] + eps_spp[spp[i]] + eps_y[year_ID[i]] + eps_spy[spy[i]], 0}); 
	  // pointwise log likelihood
	  log_lik[i] = lognormal_lpdf( next_bm[i] | log(predBM), sigma);
  }
}

