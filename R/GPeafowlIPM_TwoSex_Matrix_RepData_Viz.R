
# Test Run 

library(coda)
library(ggplot2)
library(reshape2)
library(viridis)
library(tidyr)
library(plyr)


# Set up Age and time 

Amax <- 4
Tmax <- 8

# Load MCMC Results

rho_norm <- readRDS('PeafowlIPM_TwoSex_Matrix_TestRun_ComRep_Clutch_Origin.rds')
rho_rep <- readRDS('PeafowlIPM_TwoSex_Matrix_TestRun_BreedProb.rds')

str(rho_norm)
str(rho_rep)


rho_norm <- as.matrix(rho_norm)
rho_norm.data <- melt(rho_norm)
colnames(rho_norm.data) <- c('Sample', 'Parameter', 'Value')
rho_norm.data$Models <- 'Add clutch'


rho_rep <- as.matrix(rho_rep)
rho_rep.data <- melt(rho_rep)
colnames(rho_rep.data) <- c('Sample', 'Parameter', 'Value')
rho_rep.data$Models <- 'Add clutch + breed prob'

mains <- c(paste('rho[',c(1:8), ']', sep = ''), 'mean.rho')
# mains <- c(paste('Fec[',c(1:8), ']', sep = ''))
# mains <- c(paste('NBreedM[', c(1),',',' ', c(1:Tmax), ']', sep = ''))


rhoNormal <- subset(rho_norm.data, Parameter %in% mains)
rhoRepAdd <- subset(rho_rep.data, Parameter %in% mains)

all.d <- list(rhoNormal, rhoRepAdd)
all.da <- dplyr::bind_rows(rhoNormal, rhoRepAdd)

str(all.da)
rm(all.d, rho_norm, rho_norm.data, rho_rep, rho_rep.data)
ggplot(all.da) + 
  geom_density(aes(group = Models, x = Value, 
                   color = Models, fill = Models), alpha = 0.5) +
  facet_wrap(~Parameter, scales = 'free') +
  scale_fill_viridis(discrete = T) + 
  scale_color_viridis(discrete = T) + 
  theme_bw() + theme(panel.grid = element_blank())



