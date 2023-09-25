#------------------------------------------------------------------------------#
########################### Figures of alpha matrix ############################
################ correspond to fig.1 (a,b,c,d) in the manuscript ###############
#### this script is for dataset 1, for the other the analysis is the same ###
################## but with input objects of dataset 2 ###################
#------------------------------------------------------------------------------#

rm(list=ls())

library(ggplot2)
library(reshape2)
library(scales) # for change the transparency of colors
library(qgraph) # for the network figure
library(patchwork) # for arrange the plot
library(gridExtra) # annotation below x axis


################################### plotting ###################################
# the results of model fitting are stan_fit object
# the pairwise model selected using loo is the m.7.stan in the model.stan folder
# in the list of model formulation PDF file, it is the m.7 model depicted in equation 

# the main part of model output is a 12000*637 (number of sampling * number of parameters) matrix
# which consists of sampled distribution of all parameters in the model
# for detail of stan_fit output see https://mc-stan.org/rstan/reference/stanfit-class.html

# the alpha_matrix.csv and pars.csv which consists of other estimated parameters in the pairwise model
# were mean estimates extracted from the matrix in the stan_fit object
################################################################################


# Figure 2(a) plots alpha matrix 
# load interaction strength along with their posterior distributions
# load parameter estimates from the model
alpha_matrix <- as.matrix(read.csv("alpha_matrix_1.csv", row.names = 1, header = T))

alpha_post <- as.matrix(read.csv("alpha_1_posterior.csv", row.names = 1, header = T))

pars <- as.matrix(read.csv("pars_1.csv", row.names = 1, header = T))


# change to dataframe required by ggplot
alpha_df <- melt(alpha_matrix)

p1 <- ggplot(alpha_df, aes(x = Var2, y = Var1)) + 
  geom_raster(aes(fill=value)) + 
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(alpha_df$Var1))) +
  scale_colour_gradient2(
    low = "#cb181d",     
    high = "#005e7d",             
    limits = c(-1.5,3),
    n.breaks = 9,
    space = "Lab",
    guide = "legend",
    name = "alpha",
    aesthetics = "fill") +
  labs(x="", y="", title = "Dataset 1", tag = "(a)") +
  geom_text(aes(label = round(value, 2))) +
  theme_light() + 
  theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11, hjust=0.5, face="bold")) 

ggsave("alpha_matrix.png", p1, units="cm", width=12, height=11, dpi=600)
dev.off()


# Figure 2(b) interaction network 
# arrow pointed species are species being affected, hence i in alpha_i_j

# color for positive interaction 
col.pos <- alpha("#023858", alpha = 0.8)

# color for negative interaction 
col.neg <- alpha("#cb181d", alpha = 0.7)

# output figure is called qgraph.png
p2 <- qgraph(t(alpha_matrix), layout = "circle", 
       weighted=TRUE,
       diag = TRUE,
       directed = TRUE, 
       arrows = TRUE,
       labels=colnames(alpha_matrix), 
       label.cex=1.3,
       label.color = "#636363",
       label.prop=0.8,  
       border.color = "#636363",
       border.width = 3,
       posCol= col.pos,       
       negCol= col.neg, 
       # the figure created is two figures overlapped each other
       # one "fade = T", the other "fade = F"
       # to accentuate the faded color
       fade = T, 
       esize = 6,
       edge.width = 4,
       curve = 1.3, 
       curveAll = TRUE,
       loop = 1,
       asize = 4,
       lty = "solid",
       title = "(b)",
       title.cex = 3,
       filetype = "png", width = 12, height = 11)
       
dev.off()


# Figure 2(c) + 2(d): plot alpha_i_j versus alpha_j_i
i_j <- as.vector(alpha_matrix[upper.tri(alpha_matrix)])
j_i <- as.vector(alpha_matrix[lower.tri(alpha_matrix)])
df_i_j <- data.frame(x = i_j, y = j_i)

p3 <- ggplot(data=df_i_j, aes(x=i_j, y=j_i))+
  geom_hline(yintercept=0, lty="dashed")+
  geom_vline(xintercept=0, lty="dashed") +
  geom_point(size=4, col="grey", alpha=0.7) +
  geom_point(size=3.5, col="white") +
  geom_point(size= 3, col="#368B85", alpha=0.9)+   
  theme_classic()+
  labs(x= ~ paste(alpha[i][","][j]),
       y = ~ paste(alpha[j][","][i]),
       tag = "(c)",
       title="Non-reciprocal Inter-specific Interactions") +
  theme(plot.title = element_text(hjust=0.5,face="bold",size=10),
        axis.text=element_text(size=10),
        axis.title=element_text(size=14))


# use alpha posterior distribution to display the difference
# between intra and interspecific interaction coefficients
alpha_intra <- alpha_post[, seq(1, 64, by=9)]
alpha_inter <- alpha_post[ ,!colnames(alpha_post) %in% colnames(alpha_intra)]
intra_mean <- apply(alpha_intra, 1, function(x) mean(x))
inter_mean <- apply(alpha_inter, 1, function(x) mean(x))

col.non <- alpha("#78938A", alpha=0.8)
col.diag <- alpha("#B270A2", alpha=0.8)

alpha_mean <- data.frame(mean=c(intra_mean, inter_mean), 
                         group=rep(1:2, each=length(intra_mean)))

alpha_mean$group <- as.factor(alpha_mean$group)

p4 <- ggplot()+
  geom_density(alpha_mean, mapping=aes(x=mean,col=group),linewidth=1.5)+
  labs(x="Interaction coefficients",
       y= "Density",
       title = "Posterior Distribution of Inter- and Intra-specific Interactions",
       tag="(d)",
  ) +
  scale_color_manual(values=c(col.diag, col.non),    
                     name = "",
                     labels=c("Intra-specific", "Inter-specific")) +
  geom_vline(xintercept=mean(inter_mean), color= col.non, lty=1, linewidth=1.2) +
  geom_vline(xintercept=mean(intra_mean), color= col.diag, lty=1, linewidth=1.2) +
  theme_classic() +
  theme(legend.position = c(0.85, 0.9),
        legend.key.size = unit(0.8,"cm"),
        plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
        axis.title = element_text(size=14))

p0 <- ggplot() + theme_void()
p <- p3 +p0 + p4 + plot_layout(widths = c(0.9,0.03, 1.07))

ggsave("set1_3_4.png", p, units="cm", width=24, height=11.5, dpi=600)

dev.off() 


#------------------ for the supplementary Fig.7 -------------------------------#
## difference between inter- and intra-specific interaction
dif <- inter_mean - intra_mean

## creating a global x axis label
global_title <- ggplot() +
  annotate(geom = "text", x=1, y=1, 
           label="Statistical Test for the Difference") +
  theme_void()

## p5 is for dataset 1
p5 <- ggplot() + 
     geom_density(aes(x=dif),linewidth=1.5, colour="#E96479")+
  labs(x="Difference between inter- and intra-specific interactions",
       y= "Density",
       title = "Dataset 1",
       tag="(a)",
  ) +
  geom_vline(xintercept=0, color= "#E96479", lty=2, linewidth=1) +
  annotate("text", x = mean(dif), y = 2, label = "100%", size = 3) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size=9),
        axis.title.x = element_text(size = 9))

ggsave("Sup_fig_7.png", p5, units="cm", width=12, height=11, dpi=600)
dev.off()




