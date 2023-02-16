###############################################################################
# Code for the article '....'
# 
# Copyright (C) 2023 Wentao Yu 
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY.  See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################


library(reshape2)
library(matrixStats)
library(plyr)
library(splitstackshape)

load("BEF_data_idx.RData")

# select data based on plot ID from site A 
# compile 8 species data from 1~8 diversity gradient from site A
# plotID_1 represent site while plotID_2 denote site ID
BEF_data_A <- BEF_data[which(BEF_data$plotID_1 == "A"), ]

# check  8 species name (exclude misplanted species)
sort(unique(BEF_data_A[BEF_data_A$plotID_2 == "T10", ]$sp))


# specified species to be selected according to the eight species mixture
sp_name <-  c("Castanea henryi" , "Castanopsis sclerophylla", "Choerospondias axillaris",
              "Liquidambar formosana",  "Nyssa sinensis", "Quercus serrata",         
              "Sapindus mukorossi",  "Sapium sebiferum"  )

richness <- c(1,2,4,8)


# compile 8 species data from 1~8 diversity gradient from site A
# this method will select other compositions than only targeted mixtures
list_8 <- list()
for (i in 1:4){
  
  sub <- list()
  
  for (j in 1:length(sp_name)){
    
    sub[[j]] <- BEF_data_A[which(BEF_data_A$plot_richness == richness[i] & BEF_data_A$sp == sp_name[j]), ]
    
  }
  list_8[[i]] <- do.call(rbind, sub)
}

BEF_8 <- do.call(rbind, list_8)


# check how many plots selected
length(unique(BEF_8$plotID_2))  
A_Plots_full <- as.character(unique(BEF_8$plotID_2))


# select plots that only consist of specified species (won't select plots with misplanted species)
compo_A_list <- split(BEF_data_A, BEF_data_A$plotID_2, drop = TRUE)
plot_include <- c()
for (i in 1:length(compo_A_list)){
  
  if(all(unique(compo_A_list[[i]]$sp) %in% sp_name) == TRUE){
    plot_include[i] = as.character(unique(compo_A_list[[i]]$plotID_2))
  } else {
    plot_include[i] = NA
  }
}

plot_include <- plot_include[!is.na(plot_include)]


# check mismatched plots, potentially contain plots with misplanted species
mismatch_plot <- A_Plots_full[ ! A_Plots_full %in% plot_include]


# creating a table with species and corresponding tree numbers
# to identify the plots with misplants which should also be included
data_mismatch <- BEF_data_A[(BEF_data_A$plotID_2 %in% mismatch_plot),  ]
mism_list <- split(data_mismatch, data_mismatch$plotID_2, drop = TRUE)
A_sp_treenum <- list()
for (i in 1:length(mism_list)){
  
  tb <- count(mism_list[[i]], "sp")
  tb_order <- tb[order(tb$sp), ]
  A_sp_treenum[[i]] = data.frame(unique(mism_list[[i]]$plotID_2),
                                 unique(mism_list[[i]]$plot_richness),
                                 sort(unique(mism_list[[i]]$sp)), tb_order$freq)
  
}

A_sp_num <- do.call(rbind, A_sp_treenum) 
colnames(A_sp_num)[1]<-"plotID_2"
colnames(A_sp_num)[2]<-"richness"
colnames(A_sp_num)[3]<-"sp"
colnames(A_sp_num)[4]<-"tree_num"
write.csv(A_sp_num, file="mismatch_plot_sp_num.csv")


# visually inspect the plots should be included due to misplanting
# it is not ideal but this ensures we don't lose plots
misplant_plot <- c("F27", "G26", "P19", "Q20", "R16")

# compile the final version of plots that should be included in the study
plots <- c(plot_include, misplant_plot)
BEF_A_8 <- BEF_data_A[(BEF_data_A$plotID_2 %in% plots), ]


# identified misplanted species
sp_all <- unique(BEF_A_8$sp)
sp_misplant <- sp_all[! sp_all %in% sp_name]

# now to deal with misplanted species, set their biomass to -1
# so that later those data points will be excluded
BEF_exp <- BEF_A_8[(BEF_A_8$sp %in% sp_name), ]
BEF_mis <- BEF_A_8[!(BEF_A_8$sp %in% sp_name), ]  
BEF_mis[ , 8:14] <- -1
BEF_A_8 <- rbind(BEF_exp, BEF_mis)







#####################################################################################
## create new structure to exclude invalid data point and aid the log_lik computation
#####################################################################################
BEF_focal <- BEF_A_8[BEF_A_8$focal, ]
BEF_focal <- BEF_focal[ , c(2:3, 6,8:14, 29:36, 45:47)]
NB_info <- as.matrix(BEF_focal[, 11:18])  ## select neighbor tree information
BM_focal <- as.matrix(apply(BEF_focal[, 4:10], 2, as.numeric))

BM_melt <- melt(t(BM_focal[ , -7]), value.name = "biomass")
BM_melt_next <- melt(t(BM_focal[ , -1]), value.name = "biomass_next")
colnames(BM_melt)[2] <- "ID_focal"  ## correspond to the row of BM_focal
BM_melt$biomass_next <- BM_melt_next$biomass_next  ## add biomass of the next year

GR <- rowDiffs(BM_focal) ## matrix of obs growth rate
GR_melt <- melt(t(GR), value.name = "growth")
BM_melt$growth <- GR_melt$growth


# create neighbour biomass matrix (exclude the last year)
# for each dataframe contain neighbour info for one focal tree, 8 row (neighbour) * 6 biomass
NB_BM_list <- list()
for (i in 1:nrow(NB_info)){ 
  NB <- list()
  for (j in 1:ncol(NB_info)){
    NB[[j]]  = BEF_A_8[ BEF_A_8$treeID == NB_info[i, j] , 8:13]
  } 
  NB_BM_list[[i]] <- do.call(rbind, NB)
}


# switch to 6 biomass *8 row (neighbour) in the order of focal tree
NB_list_t <- list()
for (i in 1:length(NB_BM_list)){
  NB_list_t[[i]] <- t(NB_BM_list[[i]])
}

NB_BM <- do.call(rbind, NB_list_t)
colnames(NB_BM) <- c("NB1", "NB2", "NB3", "NB4", "NB5", "NB6", "NB7", "NB8")


# create neighbour tree species ID matrix focal tree number * 8
NB_spID_list <- list()

for (i in 1:nrow(NB_info)){ 
  NB_sp <- list()
  for (j in 1:8){
    NB_sp[[j]]  = BEF_A_8[ BEF_A_8$treeID == NB_info[i, j], "sp" ]
  } 
  NB_spID_list[[i]] <- do.call(cbind, NB_sp)
}
NB_spID <- do.call(rbind, NB_spID_list)

colnames(NB_spID) <-  c("sp_NB1","sp_NB2","sp_NB3","sp_NB4","sp_NB5", "sp_NB6","sp_NB7","sp_NB8")


# attach richness, species, plot info to BM_data
info <- BEF_focal[ ,c("plot_richness", "sp", "plotID_1", "plotID_2")]


# replicate the info 6 times for merging
info_rep <- expandRows(info, count=6, count.is.col = FALSE)


# replicate spID matrix as well
NB_spID_rep <- expandRows(NB_spID, count = 6, count.is.col = FALSE)

# concatenate focal biomass with their neighbour biomass, neigbour spID and info
BM_A_8_data <- cbind(BM_melt, NB_BM, NB_spID_rep, info_rep)
save(BM_A_8_data, NB_BM, NB_spID, file = "BM_A_8_data.RData") # no filtering






##############################################################################
############ select positive growth and create stan data structure ###########
##############################################################################
# select the data points should be included
# biomass of the previous year and the growth should be > 0
BM_data_filter <- subset(BM_A_8_data, subset = biomass>0.0 & growth>0.0)


# biomass of neighbours should be >= 0 
BM_data_filter <- subset(BM_data_filter, subset = NB1>=0 & NB2>=0 &
                           NB3>=0 & NB4>=0 & NB5>=0 & NB6>=0 & NB7>=0 & NB8>=0 )


# replace species name with number
levels(as.factor(BM_data_filter$sp))
# levels: 1 ~ 8
sp_level <- c( "Castanea henryi" ,"Castanopsis sclerophylla", "Choerospondias axillaris",
               "Liquidambar formosana" , "Nyssa sinensis", "Quercus serrata",         
               "Sapindus mukorossi", "Sapium sebiferum")

match("Nyssa sinensis", sp_level ) ## check if it works

# replace the neighbour species with levels (numeric)
for (i in 1:nrow(BM_data_filter)){
  for (j in 1:8){ # match to find first occurance (here only occurance)
    BM_data_filter[ i, 13+j ] <- match(BM_data_filter[ i, 13+j ], sp_level) 
  }
}


# split the year identifier
BM_data_filter <- cSplit(BM_data_filter, splitCols = "Var1", sep = "_" )
BM_data_filter <- BM_data_filter[ , -25]
colnames(BM_data_filter)[25] <- "year"

# add plot*species index to each data point
BM_data_filter$spp <-  with(BM_data_filter, paste0(sp, " ",  plotID_2))

# add year*species index to each data point
BM_data_filter$spy <- with(BM_data_filter, paste0(sp, " ", year))


# create the data structure for stan_model
ND <- nrow(BM_data_filter)    ## number of total data points
NS <- length(unique(BM_data_filter$sp)) ## number of species
NN <- 8   ## number of neighboring tree
NP <- length(unique(BM_data_filter$plotID_2))
NR <- length(unique(BM_data_filter$plot_richness))

focal_bm <- BM_data_filter$biomass  ## biomass of focal tree
next_bm <- BM_data_filter$biomass_next
neib_bm <- as.matrix(apply(BM_data_filter[ , 6:13], 2, as.numeric)) ## corresponding biomass of neighbour trees
fsp_ID <- as.numeric(as.factor(BM_data_filter$sp))
neib_sp_ID <- as.matrix(apply(BM_data_filter[ , 14:21], 2, as.numeric)) ## species ID for neighbour trees ND*8
fplot_ID <- as.numeric(droplevels.factor(BM_data_filter$plotID_2)) ## plot ID for focal trees
fdiv_ID <- BM_data_filter$plot_richness
NY <- 6
N_spy <- length(unique(BM_data_filter$spy))
N_spp <- length(unique(BM_data_filter$spp))
year_ID <- BM_data_filter$year
spy <- as.numeric(as.factor(BM_data_filter$spy))
spp <- as.numeric(as.factor(BM_data_filter$spp))


# data structure for model without interactions
data.1 = list(ND = ND,  
                NS = NS,  
                NN = NN,
                focal_bm = focal_bm,
                next_bm = next_bm,
                fsp_ID = fsp_ID)


# data structure for model with average interactions
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


# data structure for model with pairwise interactions
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