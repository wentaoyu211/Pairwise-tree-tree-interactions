# Pairwise-tree-tree-interactions 
(Codes will become available upon the acceptance of this manuscript)

by Wentao Yu, Georg Albert, Benjamin Rosenbaum, Florian Schnabel, Helge Bruelheide, John Connolly, Werner Härdtle, Goddert von Oheimb, Stefan Trogisch, Nadja Rüger, Ulrich Brose

This repository contains codes for data processing, statistical analysis, simulations that perturb the interaction matrix, identifiability analysis. It modeled individual tree growth using empirical data (forest inventory data) from BEF-China spanning from 2009 to 2016. It partitioned the individual tree growth into intrinsic growth rate described by metabolic theory (MTE) and identity specific pairwise tree-tree interactions with direct neighbours. 

## Code in this repository
  1. **code/data_processing.r:** Code to select 8 species included in the study from a full dataset which contains all 40 tree species in BEF-China. Subsequently, 
  suitable data structures were created for fitting the models in rstan.
  2. **Code/model_array.r:** Code to fit the model candidates with specs (i.e. iteration, acceptance).
  3. **Code/model.stan:** A folder contains all the model formulations tested in the study.
  4. **Code/params_file.csv:** Parameter files for running the array job on HPC. Files of 2,3,4 are need to run the models.
  5. **Code/analysis.r:** Code to extract alpha (interaction coefficients) matrices and visualize the results.
  6. **Code/alpha_sim.r:** Code to perturb the interaction matrices and comptuing corresponding biomass.
  7. **Code/identifiability.r** Code to perform indentifiability analysis given the number of parameters in the pairwise model.
