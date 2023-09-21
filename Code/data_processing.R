#------------------------------------------------------------------------------#
# original data: BEF-China inventory data of site A and B from 2009-2016  
# select data from mono culture to eight species mixture from site A
# because site B has very high mortality rate, thus less data points
# those chosen specific two sets of eight species mixture has the best data quality
# see the supplementary for selected species 
# this script is for one set data (from mono to eight species mixture used in the main text)
# for the other set of the data shown in the supplementary
# change dat_set_1.csv to dat_set_2.csv 
# R version 4.2.1
#------------------------------------------------------------------------------#

rm(list=ls())

################################################################################
# prepare the correct data structure for Stan model 
################################################################################

# read the data of set 1
dat <- read.csv("dat_set_1.csv", header = T)

# create the data structure for stan_model
ND <- nrow(dat)    # number of total data points

NS <- length(unique(dat$sp)) # number of species

NN <- 8 # number of neighboring tree (without mortality)

NP <- length(unique(dat$plotID_2)) # number of plots

NR <- length(unique(BM_data_filter$plot_richness)) # number of richness level

focal_bm <- dat$biomass  # biomass of focal tree

next_bm <- dat$biomass_next # biomass of the next year

neib_bm <- as.matrix(apply(dat[ , grep("^NB*", colnames(dat))], 2, as.numeric)) # corresponding biomass of neighbour trees

fsp_ID <- as.numeric(as.factor(dat$sp)) # the species identity of focal tree

neib_sp_ID <- as.matrix(apply(dat[ , grep("^sp_NB*", colnames(dat))], 2, as.numeric)) # species ID for neighbour trees ND*8

fplot_ID <- as.numeric(droplevels.factor(dat$plotID_2)) # plot ID for focal trees

fdiv_ID <- dat$plot_richness # the richness of the plot where the focal tree is located

NY <- 6 # total years

N_spy <- length(unique(dat$spy)) # total number of species-year combination

N_spp <- length(unique(dat$spp)) # total number of species-plot-identifier combination

year_ID <- dat$year # year for each data point

spy <- as.numeric(as.factor(dat$spy)) # species-year-identifier

spp <- as.numeric(as.factor(dat$spp)) # species-plot-identifier


# data for the null model
data.1 = list(ND = ND,  
                NS = NS,  
                NN = NN,
                focal_bm = focal_bm,
                next_bm = next_bm,
                fsp_ID = fsp_ID)


# data for the neutral model
data.5 = list(ND = ND,  
              NS = NS,  
              NP = NP,
              NY = NY,
              N_spy = N_spy,
              N_spp =N_spp,
              focal_bm = focal_bm,
              next_bm = next_bm,
              fsp_ID = fsp_ID,
              year_ID = year_ID,
              spy = spy,
              spp = spp,
              fplot_ID= fplot_ID)  


# data for the pairwise interaction model
data.6 = list(ND = ND,  
               NS = NS,  
               NN = NN,
               NP = NP,
               NY = NY,
               N_spy = N_spy,
               N_spp =N_spp,
               focal_bm = focal_bm,
               next_bm = next_bm,
               neib_bm = neib_bm,
               neib_sp_ID = neib_sp_ID,
               fsp_ID = fsp_ID,
               year_ID = year_ID,
               spy = spy,
               spp = spp,
               fplot_ID= fplot_ID)


save(data.1,  data.5, data.6, file="A_8sp.RData")

# Note that A_8sp.RData are input data files for Stan models and it should run
# together with params_file.csv, model.array.R, and the stan file in model.stan folder 
# the study ran the model on HPC using 12 cores for each model
# the time used for each model ranged from around 4 hours to 14 days.
