
# Test Run 

library(coda)
library(ggplot2)
library(reshape2)
library(viridis)
library(tidyr)
library(plyr)


# Set up Age and time 

Amax <- 4
Tmax <- 23

# Load MCMC Results

rho_norm <- readRDS('R/PeafowlIPM_TwoSex_Matrix_TestRun_ComRep_Origin.rds')
rho_rep <- readRDS('R/PeafowlIPM_TwoSex_Matrix_TestRun_ComRep_Reproduction.rds')

str(rho_norm)
str(rho_rep)


# out.norm <- rho_norm
# out.rep <- rho_rep
# 
# comb <- list(c(out.norm, out.rep))
# 
# data.list <- list()
# 
# for(n in 1:2) {
#   out.sam <- as.matrix(comb)
#   out.data <- reshape2::melt(out.sam)
#   colnames(out.data) <- c('Sample', 'Parameter', 'Value')
# }


rho_norm <- as.matrix(rho_norm)
rho_norm.data <- melt(rho_norm)
colnames(rho_norm.data) <- c('Sample', 'Parameter', 'Value')
rho_norm.data$Models <- 'Original'


rho_rep <- as.matrix(rho_rep)
rho_rep.data <- melt(rho_rep)
colnames(rho_rep.data) <- c('Sample', 'Parameter', 'Value')
rho_rep.data$Models <- 'BroodSize'

mains <- c(paste('rho[',c(1:23), ']', sep = ''), 'mean.rho')

rhoNormal <- subset(rho_norm.data, Parameter %in% mains)
rhoRepAdd <- subset(rho_rep.data, Parameter %in% mains)

all.d <- list(rhoNormal, rhoRepAdd)
all.da <- dplyr::bind_rows(rhoNormal, rhoRepAdd)

str(all.da)
rm(all.d, rho_norm, rho_norm.data, rho_rep, rho_rep.data)
ggplot(all.da) + 
  geom_density(aes(group = Models, x = Value, 
                   color = Models, fill = Models), alpha = 0.5) +
  xlim(1,10) +
  facet_wrap(~Parameter, scales = 'free') +
  scale_fill_viridis(discrete = T) + 
  scale_color_viridis(discrete = T) + 
  theme_bw() + theme(panel.grid = element_blank())



