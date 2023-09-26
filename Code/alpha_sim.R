#------------------------------------------------------------------------------#
################ simulation study based on alpha matrix ########################
#------------------------------------------------------------------------------#
rm(list=ls())

library(tidyverse)
library(plyr)
library(ggh4x)
library(gridExtra)

################################################################################
####################### perturbing alpha matrices ##############################
################## for alpha matrix estimated using dataset 1 ##################
#####################   for dataset 2 change input file  #######################
###   to test run the code, set the n_shuffle to e.g.10 instead of 1000      ###
###              the 1000 run takes almost one day on HPC                    ###
################################################################################

# load parameter estimates 
alpha_matrix <- as.matrix(read.csv("alpha_matrix_1.csv", row.names = 1, header = T))

pars <- as.matrix(read.csv("pars_1.csv", row.names = 1, header = T))


## function to summarize data
# - data as data frame
# - varname to specify the column names of the variables to summarize
# - groupnames to specify how the data should be grouped
data_summary <- function(data, varnames, groupnames){
  # modified from http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      
      median = median(x[[col]], na.rm=TRUE),
      quantile25 = quantile(x[[col]], probs = 0.25, na.rm=TRUE),
      quantile75 = quantile(x[[col]], probs = 0.75, na.rm=TRUE),
      
      n = length(x[[col]][!is.na(x[[col]])]))
  }
  
  data_sum <- c()
  for(i in 1:length(varnames)){
    data_sum <- rbind(data_sum, 
                      ddply(data, groupnames, .fun=summary_func, varnames[i]))
  }
  
  data_sum <- cbind(varname = rep(varnames, each = nrow(data_sum) / length(varnames)), data_sum)
  data_sum <- rename(data_sum, c("quantile25.25%" = "quantile25",
                                 "quantile75.75%" = "quantile75"))
  return(data_sum)
}


## identify direct neighbours
# - df_cell is a dataframe that has treeIDs, x, and y coordinates - created in simulation loop
# - ind is identifier for individual trees in plot
# - indID can be used if values in ind don't use the format "ind1"
# - wrap specifies if tree interactions should be considered with periodic boundary conditions (TRUE) or not (FALSE)
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
    {cbind(.[, indID], 1)} %>%
    as.data.frame()
  colnames(tmp) <- c(indID, "dist")
  
  
  ## prepare output
  out <- suppressMessages(left_join(df_cell[indID], tmp)[, 2])
  
  
  return(as.numeric(out))
}


# sample values making sure each value is sampled an equal number of times
# - vals are the values to sample from
# - n is the total number of samples
sample_equal <- function(vals, n){
  out <- rep(NA, n)
  
  for(i in vals){
    out[sample((1:n)[is.na(out)], n/length(vals), replace = FALSE)] <- i
  }
  return(out)
}




# _ 0.1 data ####
vec_beta_sp <- c(pars[grep("^beta", row.names(pars)), "mean"])
theta <- pars["theta", "mean"]
b <- pars["b", "mean"]



# 1 shuffle alpha matrix ####
# _ 1.0 setup ####

# specify number of reshuffling
# for trying out the code, change the number to 10
# otherwise it will take too long
n_shuff <- 1000




# _ 1.1 unconstrained ####

# initialize output
sim_alpha_uc <- list()

# shuffle alphas
for (i in 1:n_shuff){
  sim_alpha_uc[[i]] <- matrix(sample(c(alpha_matrix)), ncol=8, nrow=8)
}



# _ 1.2 constrained ####

# initialize output
sim_alpha_c <- list()

# get indices for non diagonal values
j = row(alpha_matrix) != col(alpha_matrix) 

# shuffle alphas
for (i in 1:n_shuff){
  # creat empty matrices
  sim_alpha <- diag(0, ncol = 8, nrow = 8)
  # fill in the non-diagonal with non-diagonal values
  sim_alpha[j] <- ave(alpha_matrix[j], row(alpha_matrix)[j], FUN = sample)
  # fill in the diagonal with diagonal values
  diag(sim_alpha) <- sample(diag(alpha_matrix))
  
  sim_alpha_c[[i]] <- sim_alpha
}




# _ 1.3 combine with no interaction and real alpha matrix ####

# create no interaction matrix:
no_inter <- matrix(0, 8, 8)

# combine matrices in list and create reference variable for the shuffling scenarios
sim_alpha_list <- c(list(no_inter, alpha_matrix), sim_alpha_uc, sim_alpha_c)
scenario <- c("no_inter", "real_inter",
              rep("unconstrained", length(sim_alpha_uc)),
              rep("constrained", length(sim_alpha_c)))







# 2. simulate tree growth ####
# _ 2.0 setup ####

# - list of species compositions:
sp_comp <- list(1, 2, 3, 4, 5, 6, 7, 8, c(1,2), c(3,4), c(5,6), c(7,8), c(1, 2, 3, 4), c(5, 6, 7, 8), c(1, 2, 3, 4, 5, 6, 7, 8))

# - number of individuals (12x12 for plot-richness > 2, 6x6 otherwise):
nind <- sapply(sp_comp, function(x) ifelse(length(x) > 2, 12*12, 6*6))

# - total number of species considered:
nsp <- unlist(sp_comp) %>% unique() %>% length()

# - years of simulation:
years <- 7

# - number of plots per composition (i.e. repetitions):
reps <- 3


# - specify interaction type
interaction_type = 1
# 0: /
# 1: sum(alpha_ij * B_j^b)
# 2: sum(alpha_ij * (B_j * B_i)^b)
# 3: sum(alpha_ij * (B_j / B_i)^b)


# - specify how interactions affect intrinsic growth
interaction_effect = 1
# 1: additive, i.e. intrinsic + interaction
# 2: multiplicative, i.e. intrinsic * (1 + interaction)





# _ 2.1 simulation ####

out_long <- c()

for(alpha_sel in 1:length(scenario)){
  
  
  
  df_alpha_sp <- sim_alpha_list[[alpha_sel]]
  
  
  
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
      spID <- sample_equal(sp_comp[[i]], nind[i])
      
      
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
      # biomass_df <- as.data.frame(matrix(rlnorm(nind[i], meanlog = 4, sdlog = 1.1), nrow = 1))
      biomass_df <- as.data.frame(matrix(rep(100, nind[i]), nrow = 1))
      colnames(biomass_df) <- paste0("ind", indID_range)
      
      
      intrinsic_df <- c()
      inter_df <- c()
      growth_df <- c()
      
      ## initialize output df for neighbours
      neigh_out <- data.frame() 
      
      
      ## simulate growth over years 
      for(k in 1:(years-1)){
        
        biomass <- as.numeric(biomass_df[k, ])
        
        ## calculate intrinsic growth term
        intrinsic_growth <- vec_beta_ind * biomass^theta
        
        intrinsic_df <- rbind(intrinsic_df, intrinsic_growth)
        
        
        
        ## calculate interaction term
        if(interaction_type == 0){
          interaction_sum <- rep(0, nind[i])
        } else if(interaction_type == 1){
          interaction_sum <- (df_alpha_ind * df_neigh) %>% apply(1, function(x){x * biomass^b}) %>% colSums(na.rm = TRUE)
          # use colSums because of apply!
        } else if(interaction_type == 2){
          interaction_sum <- (df_alpha_ind * df_neigh * sapply(biomass, function(x) x * biomass)^b) %>% rowSums(na.rm = TRUE)
        } else if(interaction_type == 3){
          interaction_sum <- (df_alpha_ind * df_neigh * sapply(biomass, function(x) x / biomass)^b) %>% rowSums(na.rm = TRUE)
        }
        
        inter_df <- rbind(inter_df, interaction_sum)
        
        
        
        ## calculate new biomass
        if(interaction_effect == 1){
          biomass_df[k+1, ] <- biomass + intrinsic_growth + interaction_sum
          growth_df <- rbind(growth_df, intrinsic_growth + interaction_sum)
        } else if(interaction_effect == 2){
          biomass_df[k+1, ] <- biomass + intrinsic_growth * (1 + interaction_sum)
          growth_df <- rbind(growth_df, intrinsic_growth * (1 + interaction_sum))
        }
        
        ## set negative biomasses to zero (rare!) & and assure dead trees stay dead
        biomass_df[k+1, biomass_df[k+1, ] <= 0] <- 0.0
        biomass_df[k+1, biomass_df[k, ] <= 0] <- 0.0
        
        
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
      # - biomass data
      biomass_df$year <- 1:years
      biomass_df$plotID_1 <- rep("A", years)
      biomass_df$plotID_2 <- rep(paste0(i, "_", j), years)
      biomass_df$plot <- paste0(biomass_df$plotID_1, biomass_df$plotID_2)
      biomass_df$plot_richness <- rep(length(sp_comp[[i]]), years)
      
      out_prep <- biomass_df %>% 
        pivot_longer(cols = starts_with("ind"), names_to = "ID_focal", values_to = "biomass") %>%
        as.data.frame()
      
      
      # - add interactions
      inter_df <- as.data.frame(rbind(inter_df, NA))
      
      inter_df$year <- 1:years
      inter_df$plotID_1 <- rep("A", years)
      inter_df$plotID_2 <- rep(paste0(i, "_", j), years)
      inter_df$plot <- paste0(inter_df$plotID_1, inter_df$plotID_2)
      inter_df$plot_richness <- rep(length(sp_comp[[i]]), years)
      
      colnames(inter_df) <- colnames(biomass_df)
      
      out_prep <- inter_df %>% 
        pivot_longer(cols = starts_with("ind"), names_to = "ID_focal", values_to = "inter") %>%
        as.data.frame() %>%
        {cbind(out_prep, .["inter"])}
      
      
      # - add growth
      growth_df <- as.data.frame(rbind(growth_df, NA))
      
      growth_df$year <- 1:years
      growth_df$plotID_1 <- rep("A", years)
      growth_df$plotID_2 <- rep(paste0(i, "_", j), years)
      growth_df$plot <- paste0(growth_df$plotID_1, growth_df$plotID_2)
      growth_df$plot_richness <- rep(length(sp_comp[[i]]), years)
      
      colnames(growth_df) <- colnames(biomass_df)
      
      out_prep <- growth_df %>% 
        pivot_longer(cols = starts_with("ind"), names_to = "ID_focal", values_to = "growth") %>%
        as.data.frame() %>%
        {cbind(out_prep, .["growth"])}
      
      
      # - add neighborhood information
      colnames(neigh_out) <- c(paste0("NB", 1:8), paste0("sp_NB", 1:8))
      out_prep <- cbind(out_prep, neigh_out)
      
      out_prep$sp <- rep(spID, years)
      
      
      
      
      
      out <- rbind(out, out_prep)
      
      
    }
  }
  
  
  out2 <- data_summary(out, c("biomass", "inter", "growth"), c("plot", "plot_richness", "year"))
  out2$log.div <- log2(out2$plot_richness)
  out2$scenario <- rep(scenario[alpha_sel], nrow(out2))
  out2$alpha_sel <- rep(alpha_sel, nrow(out2))
  
  out_long <- rbind(out_long, out2)
  
}



# _ 2.2 transform data for figures ####
## (a) - main NEW !!
# - select interaction and growth data (i.e. not biomass!)
# - focus on last year
dat_main <- out_long[out_long$varname != "biomass" & out_long$year == 6, ]


# - make sure real and no interaction can be shown in constrained and unconstrained plots 
tmp1 <- tmp2 <- dat_main[dat_main$scenario %in% c("real_inter", "no_inter"), ]
tmp1$scenario <- "constrained"
tmp2$scenario <- "unconstrained"

dat_main <- rbind(tmp1, tmp2, 
                  dat_main[dat_main$scenario %in% c("constrained", "unconstrained"), ])


# - transform faceting variables to factors and rename them
dat_main$scenario <- factor(dat_main$scenario, 
                            levels = c("unconstrained", "constrained"),
                            labels = c("Unconstrained", "Constrained"))
dat_main$varname <- factor(dat_main$varname, 
                           levels = c("inter", "growth"), 
                           labels = c("Net Interactions", "Productivity"))





## (b) - summary 

# - take mean alpha values for diagonals and off-diagonals 
offdiag <- sim_alpha_list[2:(n_shuff+2)] %>%
  sapply(function(x) mean(x[row(x) != col(x)]))

diag <- sim_alpha_list[2:(n_shuff+2)] %>%
  sapply(function(x) mean(x[row(x) == col(x)]))


# - calculate diversity~effect slopes
# slope_biomass <- c()
slope_inter <- c()
slope_growth <- c()

for(i in unique(out_long$alpha_sel)){
  # slope_biomass <- c(slope_biomass, lm(mean ~ log.div, data = out_long[out_long$varname == "biomass" & out_long$alpha_sel == i, ])$coefficients[2])
  slope_inter <- c(slope_inter, lm(mean ~ log.div, data = out_long[out_long$varname == "inter" & out_long$alpha_sel == i, ])$coefficients[2])
  slope_growth <- c(slope_growth, lm(mean ~ log.div, data = out_long[out_long$varname == "growth" & out_long$alpha_sel == i, ])$coefficients[2])
  
}


dat_summary <- rbind(
  data.frame(alpha_sel = 2:(n_shuff+2),
             diff = offdiag - diag,
             slope = slope_growth[2:(n_shuff+2)],
             response = rep("productivity", n_shuff+1),
             col = c(2, rep(1, n_shuff))),
  data.frame(alpha_sel = 2:(n_shuff+2),
             diff = offdiag - diag,
             slope = slope_inter[2:(n_shuff+2)],
             response = rep("net interactions", n_shuff+1),
             col = c(2, rep(1, n_shuff)))
)



# _ 2.3 save data

write.csv(dat_main, file = "dat_main_1.csv")
write.csv(dat_summary, file =  "dat_summary_1.csv")




# 3. plotting ####
# _ 3.1 (a) diversity~mean effect ####
p1 <- ggplot(dat_main, aes(log.div, mean)) +
  # geom_hline(yintercept = 0, lty = 2) +
  geom_smooth(data = dat_main[dat_main$alpha_sel > 2, ], aes(group = as.factor(alpha_sel)),
              method = "lm", se = FALSE, linewidth = 0.2, color = "darkgrey") +
  
  # no inter:
  geom_smooth(data = dat_main[dat_main$alpha_sel == 1, ], aes(group = as.factor(alpha_sel)),
              method = "lm", se = FALSE, color = "black", lty = 2) +
  # mean shuffled:
  geom_smooth(data = dat_main[dat_main$alpha_sel > 2, ],
              method = "lm", se = FALSE, linewidth = 1.8, color = "white") +
  geom_smooth(data = dat_main[dat_main$alpha_sel > 2, ], 
              method = "lm", se = FALSE, linewidth = 1.2, color = "#023858") +
  # real - growth & interaction:
  geom_smooth(data = dat_main[dat_main$alpha_sel == 2, ], aes(group = as.factor(alpha_sel)), 
              method = "lm", se = FALSE, linewidth = 1.8, color = "white") +
  geom_smooth(data = dat_main[dat_main$alpha_sel == 2, ], aes(group = as.factor(alpha_sel)), 
              method = "lm", se = FALSE, linewidth = 1.2, color = "#8A9A5B") + 
  
  scale_x_continuous(#breaks = c(0, 1, 2, 3), 
    name = "Species Richness (log2)") +
  scale_y_continuous(#breaks = c(-10, 0, 10, 20, 30, 2600, 2700, 2800, 2900, 3000),
    name = "Community Mean Effect") +
  
  facet_grid2(varname~scenario, scale = "free", independent = "all", switch = "y") +
  labs(tag = "(a)") +
  theme_classic() +
  theme(strip.background = element_blank(), 
        strip.text = element_text(face = "bold"),
        axis.title = element_text(size = 9),
        plot.tag = element_text(size=12),
        strip.placement = "outside")





# _ 3.2 (b) alphadiff~slope ####


p2 <- ggplot(dat_summary, aes(diff, slope, col = as.factor(col), fill = factor(col))) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 0, lty = 3) +
  geom_point(size = 1, alpha = 0.5, col = "darkgrey") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(data = dat_summary[dat_summary$alpha_sel == 2, ], 
             size = 3, shape = 22, col = "white") +
  scale_color_manual(values = c("#023858", "#8A9A5B")) + # "#93C572", "#96DED1", 
  scale_fill_manual(values = c("#023858", "#8A9A5B")) + # "", "#96DED1", 
  scale_x_continuous(#limits = c(-0.8, 0.8), breaks = c(-0.8, 0, 0.8), 
    name = "Difference of Mean Inter- and Intraspecific Interactions") +
  scale_y_continuous(#limits = c(-45, 40), breaks = c(-40, -20, 0, 20, 40),
    name = "Slope of Diversity Relationship") +
  facet_grid2(response~., scale = "free", independent = "all") +
  labs(tag = "(b)") +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white", colour = "white"),
        strip.text.y = element_blank(),
        axis.title = element_text(size = 9),
        plot.tag = element_text(size=12),
        plot.margin = margin(7, 2, 2, 0, "mm"))





# _ 3.3 merge and save plots ####
p0 <- ggplot() + theme_void()

p_out <- gridExtra::grid.arrange(grobs = list(p1, p0, p2,p0), nrow = 1, widths=c(1, 0.2,0.6, 0.2))

ggsave(filename = "set1_Fig3_1000.png", plot = p_out, width = 173, height = 108, unit = "mm", dpi = 600)

