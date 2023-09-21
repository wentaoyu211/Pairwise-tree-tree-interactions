#------------------------------------------------------------------------------#
# creating input file ID for submitting array jobs on HPC EVE
# require data A_8sp.RData, models coded in Stan and a parameter files 
# every job is a run on a model with specific structure 
# null/neutral/pairwise  * two random effect scenario
# in total 20 models were tested 
#------------------------------------------------------------------------------#

rm(list=ls())

# load required libraries
library(data.table)
library(docopt)


library(rstan)
library(coda)
rstan_options(auto_write = TRUE)
options(mc.cores = 12)

load("A_8sp.RData")

# provide initial values for stan model
inits <- rep(list(list(sigma = 1.0)), 3)

doc <- "usage: array_job.R <params>"
opts <- docopt(doc)

process_single <- function(params) {
  ## actual processing goes here
  ## you can access the parameters from the parameter file
  ## like you would normally with readCSV or data.table:
  
  # to recognize data.1 as data not character
  dat <- eval(parse(text = params$data)) 
  stan_model <- params$stan_model
  
  # compile model and sampling
  stan_fit <- sampling(  stan_model(file = stan_model),
                         data = dat,
                         iter = 8000,
                         warmup = 4000,
                         thin = 1,
                         chains = 3, 
                         init = inits,
                         control = list(max_treedepth = 10,
                                        adapt_delta = 0.9))
  
  # save the output stan_fit object
  # saveRDS(stan_fit, file = res)
  return (stan_fit)

}


## read parameter file
params <- fread(opts$params)

## try to get task id
task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## check if task id exists
if (is.na(task)) {
  ## process everything
  result <- by(params, 1:nrow(params), process_single)
  attributes(result) <- NULL
  
  saveRDS(result, file.path(opts$output_dir, "everything.rds"))
} else {
  ## process single item
  result <- process_single(params[task])
  
  saveRDS(result, file.path(opts$output_dir, paste0("chunk-", task, ".rds")))
}


