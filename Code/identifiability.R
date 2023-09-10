library(tidyverse)


# 0.1 functions <--- ####


# identify direct neighbours; removed calculating distance, as we don't use it anymore !
neigh.wrap_function <- function(df_cell, ind, indID = "ind", wrap = TRUE){
  
  
  
  ## define which are edge trees and need wrap around (if we wrap!)
  if(wrap){
    
    edge <- df_cell[df_cell$cell_x %in% c(min(df_cell$cell_x), max(df_cell$cell_x)) | 
                    df_cell$cell_y %in% c(min(df_cell$cell_y), max(df_cell$cell_y)), ]
    edge$cell_x[edge$cell_x == min(df_cell$cell_x)] <- max(df_cell$cell_x) + 1
    edge$cell_x[edge$cell_x == max(df_cell$cell_x)] <- min(df_cell$cell_x) - 1
    edge$cell_y[edge$cell_y == min(df_cell$cell_y)] <- max(df_cell$cell_y) + 1
    edge$cell_y[edge$cell_y == max(df_cell$cell_y)] <- min(df_cell$cell_y) - 1
    
    edge_xfix <- edge_yfix <- edge[paste(edge$cell_x, edge$cell_y) %in% paste(rep(range(edge$cell_x), 2), 
                                                                              rep(range(edge$cell_y), each = 2)), ]
    edge_xfix$cell_x <- abs(edge_xfix$cell_x - max(df_cell$cell_x))
    edge_yfix$cell_y <- abs(edge_yfix$cell_y - max(df_cell$cell_y))
    
    
    edge <- rbind(edge, edge_xfix, edge_yfix)
    
    
    tmp.df_cell <- rbind(cbind(df_cell, core = rep("core", nrow(df_cell))), cbind(edge, core = rep("edge", nrow(edge))))
    
  } else {
    
    tmp.df_cell <- cbind(df_cell, core = rep("core", nrow(df_cell)))
    
  }
  
  
  ## specify neighbouring tree IDs !
  tmp <- tmp.df_cell %>%
    {.[.[, indID] == ind & .$core == "core", ]} %>%
    {tmp.df_cell[tmp.df_cell$cell_x %in% (.$cell_x-1):(.$cell_x+1) & 
                   tmp.df_cell$cell_y %in% (.$cell_y-1):(.$cell_y+1) &
                   paste(tmp.df_cell$cell_x, tmp.df_cell$cell_y) != paste(.$cell_x, .$cell_y), ]} %>%
    {cbind(.[, indID], 1)}
  
  
  
  ## prepare output
  out <- suppressMessages(left_join(as.data.frame(df_cell[, indID]) %>% rename(!!indID := !!names(.)), 
                                    as.data.frame(tmp) %>% rename(!!indID := V1, dist = V2))[, 2])
  
  
  return(as.numeric(out))
  
  
}





# 0.2 data - overview for parameterization ####

# model_pars <- read.csv("C:/Users/ga41quru/Nextcloud/own/PhD/Project/02 local diversity effects - 2nd order diversity effects/identifyability/pars.csv")
# input_data <- read.csv("C:/Users/ga41quru/Nextcloud/own/PhD/Project/02 local diversity effects - 2nd order diversity effects/identifyability/8sp_data_filtered.csv")
# 
# 
# ## alphas:
# alphas_overview <- model_pars[grepl("alpha", model_pars$X), ]
# alphas_overview <- alphas_overview[order(alphas_overview$mean, decreasing = TRUE), ]
# alphas_overview$rank = 1:nrow(alphas_overview)
# 
# ggplot(alphas_overview, aes(mean, rank)) +
#   geom_point() +
#   geom_errorbar(aes(xmin=X2.5., xmax=X97.5.)) +
#   geom_vline(xintercept = 0, lty = 2)
# 
# hist(alphas_overview$mean)
# 
# hist(rnorm(100))
# # I think we can just go with alphas drawn from normal distribution and then setting 10% of them to zero !
# 
# 
# 
# ## betas:
# betas_overview <- model_pars[grepl("beta", model_pars$X), ]
# betas_overview <- betas_overview[order(betas_overview$mean, decreasing = TRUE), ]
# betas_overview$rank = 1:nrow(betas_overview)
# 
# ggplot(betas_overview, aes(log(mean), rank)) +
#   geom_point() +
#   geom_errorbar(aes(xmin=log(X2.5.), xmax=log(X97.5.))) +
#   geom_vline(xintercept = 0, lty = 2)
# 
# hist(betas_overview$mean)
# hist(rnorm(8, mean = 4))
# 
# 
# 
# ## starting densities:
# dens_year1 <- input_data$biomass[input_data$year == 1]
# 
# hist(dens_year1)
# hist(log(dens_year1))
# hist(rlnorm(length(dens_year1), meanlog = 4, sdlog = 1.1)) # looks alright I'd say? (sd = 1 is a bit to short tailed?!)
# 
# 
# grow_year1 <- input_data$growth[input_data$year == 1]
# 
# hist(grow_year1)
# hist(log(grow_year1))
# hist(rlnorm(length(grow_year1), meanlog = 5, sdlog = 1)) # looks alright I'd say?
# 
# # growth also lognormal ! 
# 
# # in the model parameters, sigma (on log already) is 0.4, but that produces very different distributions !
# 
# 
# 
# 
# ## estimate mortality rates:
# sum(input_data[, 5:12] == 0) / sum(input_data[, 5:12] >= 0) / 7
# # - that's a bit bullshit, but maybe we just assume 1% mortality rates per year ! (i.e. by year 7, 6% of trees died)







# 1. simulate data ####
# + 1.1 parameters ####


## list of species compositions:
sp_comp <- list(1, 2, 3, 4, 5, 6, 7, 8, c(1,2), c(3,4), c(5,6), c(7,8), c(1, 2, 3, 4), c(5, 6, 7, 8), c(1, 2, 3, 4, 5, 6, 7, 8))


## number of individuals:
# - 12x12 for plot-richness > 2, 6x6 otherwise (same as inventory assessments BEF-China)
nind <- sapply(sp_comp, function(x) ifelse(length(x) > 2, 12*12, 6*6))


## total number of species considered:
nsp <- unlist(sp_comp) %>% unique() %>% length()


## years of simulation (we have 5-10 years):
years <- 7


## number of plots per composition (i.e. repetitions):
# - varies depending on diversity level and species composition; max. 4
reps <- 3 # currently, only one value for all compositions !


## mortality rate per year
mort <- 0.01 # 1%


## process and measurement errors:
# - both errors are sigmas in log-normal distributions
sigma_measure <- 0.05 # !!!!!!!!! according to the model parameters, this (?) is ~0.4 !!!!!!!!!!!!!
sigma_process <- 0.05 







## !!!!! parameters we want to estimate !!!!!

# 1) alpha matrix (sp x sp):
df_alpha_sp <- matrix(rnorm(nsp*nsp), ncol = nsp) # using normal distribution should be close enough (see 0.2)
df_alpha_sp[sample(1:(nsp*nsp), size = round(nsp*nsp/10), replace = FALSE)] = 0.0 # set ~10% (i.e. 6) to true zero !

colnames(df_alpha_sp) <- paste0("sp", 1:nsp)
rownames(df_alpha_sp) <- paste0("sp", 1:nsp)
# - cols: focal i; rows: competitors j --> alpha_ji
# - we assume that effect of sp1 on sp2 are independent from effects of sp2 on sp1



# 2) beta vector (sp):
vec_beta_sp  <- rnorm(nsp, mean = 4) # we use mean of 4 because it kinda looks like the one we see in the data! 
while(sum(vec_beta_sp < 0) != 0){ vec_beta_sp  <- rnorm(nsp, mean = 4) }



# 3) combined:
par <- list(
  
  alpha_spsp = df_alpha_sp,
  beta_sp = vec_beta_sp,
  theta = 0.8,
  b = 0.2
  
)







## specify which interaction type we use

interaction_type = 1
# 1: sum(alpha_ij * B_j^b)
# 2: sum(alpha_ij * (B_j * B_i)^b)
# 3: sum(alpha_ij * (B_j / B_i)^b)


## specify how interactions affect intrinsic growth

interaction_effect = 1
# 1: additive, i.e. intrinsic + interaction
# 2: multiplicative, i.e. intrinsic * (1 + interaction)









# + 1.2 simulation ####

out <- data.frame()
indID_start <- 1

for(i in 1:length(sp_comp)){
  for(j in 1:reps){
    
    ## define unique IDs for each tree:
    indID_end <- indID_start + nind[i] - 1
    indID_range <- indID_start:indID_end
    
    # - update first value (for next plot!)
    indID_start <- indID_end + 1
    
    
    ## species identifier
    if(length(sp_comp[[i]]) == 1){
      spID <- rep(sp_comp[[i]], nind[i])
    } else {
      spID <- sample(sp_comp[[i]], nind[i], replace = TRUE)
    }
    
    
    ## neighbor identifier
    df_cell <- data.frame(
      ind = paste0("ind", indID_range),
      cell_x = rep(1:sqrt(nind[i]), each = sqrt(nind[i])),
      cell_y = rep(1:sqrt(nind[i]), sqrt(nind[i]))
    ) # !!! only works with square tree layout !!!
    
    df_neigh <- sapply(paste0("ind", indID_range), function(x) neigh.wrap_function(df_cell, x))
    colnames(df_neigh) <- paste0("ind", indID_range)
    rownames(df_neigh) <- paste0("ind", indID_range)
    
    
    ## beta vector (ind)
    vec_beta_ind <- vec_beta_sp[spID]
    
    
    ## alpha matrix (ind x ind)
    df_alpha_ind <- matrix(ncol = nind[i], nrow = nind[i])
    for(k in 1:ncol(df_alpha_ind)){
      df_alpha_ind[, k] <- df_alpha_sp[, spID[k]][spID]
    }
    colnames(df_alpha_ind) <- paste0("ind", indID_range)
    rownames(df_alpha_ind) <- paste0("ind", indID_range)

    
    ## starting densities:
    biomass_df <- as.data.frame(matrix(rlnorm(nind[i], meanlog = 4, sdlog = 1.1), nrow = 1))
    colnames(biomass_df) <- paste0("ind", indID_range)

    
    ## initialize output df for neighbours
    neigh_out <- data.frame() 
    
    
    ## simulate growth over years 
    for(k in 1:(years-1)){
      
      biomass <- as.numeric(biomass_df[k, ])
      
      ## calculate intrinsic growth term
      intrinsic_growth <- vec_beta_ind * biomass^par$theta
      
  

      ## calculate interaction term
      if(interaction_type == 1){
        interaction_sum <- (df_alpha_ind * df_neigh) %>% apply(1, function(x){x * biomass^par$b}) %>% colSums(na.rm = TRUE)
          # use colSums because of apply!
      } else if(interaction_type == 2){
        interaction_sum <- (df_alpha_ind * df_neigh * sapply(biomass, function(x) x * biomass)^par$b) %>% rowSums(na.rm = TRUE)
      } else if(interaction_type == 3){
        interaction_sum <- (df_alpha_ind * df_neigh * sapply(biomass, function(x) x / biomass)^par$b) %>% rowSums(na.rm = TRUE)
      }
      
      
      ## calculate new biomass
      if(interaction_effect == 1){
        biomass_df[k+1, ] <- biomass + intrinsic_growth + interaction_sum
      } else if(interaction_effect == 2){
        biomass_df[k+1, ] <- biomass + intrinsic_growth * (1 + interaction_sum)
      }
      
      ## set negative biomasses to zero (rare!) & and assure dead trees stay dead
      biomass_df[k+1, biomass_df[k+1, ] <= 0] <- 0.0
      biomass_df[k+1, biomass_df[k, ] <= 0] <- 0.0
      
      
      ## introduce process error
      biomass_df[k+1, ] <- sapply(biomass_df[k+1, ], function(x){
          rlnorm(1, meanlog = log(x), sdlog = sigma_process)
        })
      
      
      ## extinctions
      # - negative growth is something we can find in the data (measurement error?!), so I prefer to not use this a criteria
      #   --> we exclude those cases when fitting the models on real data ! ( excluded later )
      # - each tree has a chance to die (using random binomial with mortality rate mort)
      biomass_df[k+1, which(rbinom(nind[i], size = 1, prob = mort) == 1)] <- 0.0

      
      ## save neighborhood information in output formatting:
      neigh_out <- rbind(neigh_out, cbind(
        apply(df_neigh, 1, function(x){x * biomass}) %>% apply(2, function(x){x[!is.na(x)]}) %>% t(), # neighbor's biomass
        apply(df_neigh, 1, function(x){x * spID}) %>% apply(2, function(x){x[!is.na(x)]}) %>% t() # neighbor's species ID
            ))
    }
    
    ## finish neighborhood information for last year !:
    biomass <- as.numeric(biomass_df[years, ])
    neigh_out <- rbind(neigh_out, cbind(
      apply(df_neigh, 1, function(x){x * biomass}) %>% apply(2, function(x){x[!is.na(x)]}) %>% t(), # neighbor's biomass
      apply(df_neigh, 1, function(x){x * spID}) %>% apply(2, function(x){x[!is.na(x)]}) %>% t() # neighbor's species ID
    ))
    
    
    ## create output
    biomass_df$year <- 1:years
    biomass_df$plotID_1 <- rep("A", years)
    biomass_df$plotID_2 <- rep(paste0(i, "_", j), years)
    biomass_df$plot_richness <- rep(length(sp_comp[[i]]), years)
    
    out_prep <- biomass_df %>% 
      pivot_longer(cols = starts_with("ind"), names_to = "ID_focal", values_to = "biomass") %>%
      as.data.frame()
    
    # out_prep$biomass_next <- biomass_df[2:years, ] %>% 
    #   pivot_longer(cols = starts_with("ind"), names_to = "ID_focal", values_to = "biomass") %>%
    #   {as.data.frame(.)[, "biomass"]}
    
    colnames(neigh_out) <- c(paste0("NB", 1:8), paste0("sp_NB", 1:8))
    out_prep <- cbind(out_prep, neigh_out)
    
    out_prep$sp <- rep(spID, years)
    
    
    
    
    
    out <- rbind(out, out_prep)
    
    
  }
}


# + 1.3 finalize output ####

## add measurement error

out$biomass <- sapply(out$biomass, function(x){
  rlnorm(1, meanlog = log(x), sdlog = sigma_measure)
})


## add and sort columns

# - biomass_next
out$biomass_next <- NA
for(i in 1:nrow(out)){
  if(out$year[i] != years){
    out$biomass_next[i] <- out$biomass[out$ID_focal == out$ID_focal[i] &
                                       out$year == (out$year[i] + 1)]
  }
}
  

# - transform ID_focal to numeric: 
out$ID_focal <- substring(out$ID_focal, 4) %>% as.numeric()

# - calculate growth
out$growth <- out$biomass_next - out$biomass

# - spp
out$spp <- paste(out$sp, out$plotID_2)

# - spy
out$spy <- paste(out$sp, out$year)

# - sort columns
out <- out[, c('ID_focal','biomass','biomass_next','growth','NB1','NB2','NB3','NB4','NB5','NB6','NB7','NB8','sp_NB1','sp_NB2','sp_NB3','sp_NB4','sp_NB5','sp_NB6','sp_NB7','sp_NB8','plot_richness','sp','plotID_1','plotID_2','year','spp','spy')]

# - remove misleading rownames
rownames(out) <- NULL




## remove values

# - remove edge trees
baseID_36 <- c(1:6, 7, 12, 13, 18, 19, 24, 30, 31:36)
baseID_144 <- c(1:12, 13, 24, 25, 36, 37, 48, 49, 60, 61, 72, 73, 84, 85, 96, 97, 108, 109, 120, 121, 132, 133:144)

edgeID <- c()

for(i in 1:(reps*length(sp_comp))){
  if(i <= (reps*sum(sapply(sp_comp, length) < 4))){
    edgeID <- c(edgeID, baseID_36 + (i-1)*36)
  } else {
    edgeID <- c(edgeID, (baseID_144 + 
                  (sum(sapply(sp_comp, length) < 4) * reps * 36) + 
                  ((i-(sum(sapply(sp_comp, length) < 4) * reps)-1) * 144)))
  }
}

out <- out[!(out$ID_focal %in% edgeID), ]


# - remove zero value focal trees
out <- out[out$biomass > 0, ]


# - remove negative growth or NA growth (i.e. last year)
out <- out[out$growth >= 0 & !is.na(out$growth), ]




## prepare parameter output file:
# model_pars # for reference
out_pars <- data.frame()

# - beta
for(i in 1:length(vec_beta_sp)){
  out_pars[nrow(out_pars)+1, 1] <- paste0("beta[", i, "]")
  out_pars[nrow(out_pars), 2] <- vec_beta_sp[i]
}

# - theta
out_pars[nrow(out_pars)+1, 1] <- "theta"
out_pars[nrow(out_pars), 2] <- par$theta

# - alpha
for(i in 1:nrow(df_alpha_sp)){
  for(j in 1:ncol(df_alpha_sp)){
    out_pars[nrow(out_pars)+1, 1] <- paste0("alpha[", i, ",", j, "]")
    out_pars[nrow(out_pars), 2] <- df_alpha_sp[i, j]
  }
}

# - b
out_pars[nrow(out_pars)+1, 1] <- "b"
out_pars[nrow(out_pars), 2] <- par$b

# - sigma
out_pars[nrow(out_pars)+1, 1] <- "sigma_process"
out_pars[nrow(out_pars), 2] <- sigma_process
out_pars[nrow(out_pars)+1, 1] <- "sigma_measure"
out_pars[nrow(out_pars), 2] <- sigma_measure





# + 1.4 output ####


write.csv(out, "simulated_input.csv")
write.csv(out_pars, "simulated_input_pars.csv")

# we fitted this simulated data for which we know the parameter values to the 
# pairwise interactio models and retrieved the estimated parameters
# we then compared the estimated parameter values with the values we used 
# for generating the data. Results are shown in supplementary Fig.2







