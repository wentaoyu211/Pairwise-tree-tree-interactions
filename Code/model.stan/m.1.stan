data{
  int ND; // number of total data points
  int NS; // number of species
  real focal_bm[ND];  // biomass of focal tree
  real next_bm[ND];   // the next year biomass of corresponding focal tree
  int fsp_ID[ND]; // identify species of each data point
}

parameters{
  real<lower=0> beta[NS]; // species specific growth parameter
  real theta;             // test if the growth rate scale with power law (0.75) 
  real<lower=0> sigma;    // sigma is observation error 
}

model{
  // for model simplification
  real predBM;
  
  // specify priors
  for (i in 1:NS){
    beta[i] ~ exponential(1.0); 	  
	}
  theta ~ normal(0.75, 1); // based on 3/4 exponent
  sigma ~ normal(0,1);
   
  for(i in 1:ND){    // loop over all data points
        // prediction
        predBM = focal_bm[i] + beta[fsp_ID[i]]*focal_bm[i]^theta; 
		
        // likelihood		
       next_bm[i] ~ lognormal( log(predBM), sigma );
  }
}

generated quantities {
  real predBM; 
  vector[ND] log_lik;  // number of log likelihood for loo
  
  for (i in 1:ND){    // loop through focal individual tree
      // prediction
      predBM = focal_bm[i] + beta[fsp_ID[i]]*focal_bm[i]^theta; 
	  // pointwise log likelihood
	  log_lik[i] = lognormal_lpdf( next_bm[i] | log(predBM), sigma);
  }
}

