
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

rho_norm <- readRDS('R/PeafowlIPM_TwoSex_Matrix_TestRun_ComRep_Normal.rds')
rho_rep <- readRDS('R/PeafowlIPM_TwoSex_Matrix_TestRun_ComRep_AddRep.rds')

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

mains <- c('rho', paste('rho[',c(1:23), ']', sep = ''))

rhoNormal <- subset(rho_norm.data, Parameter %in% mains)
rhoRepAdd <- subset(rho_rep.data, Parameter %in% mains)

all.d <- list(rhoNormal, rhoRepAdd)
all.da <- dplyr::bind_rows(rhoNormal, rhoRepAdd)

str(all.da)
rm(all.d, rho_norm, rho_norm.data, rho_rep, rho_rep.data)
ggplot(all.da) + 
  geom_density(aes(group = Models, x = Value, 
                   color = Models, fill = Models), alpha = 0.3) +
  xlim(1,10) +
  facet_wrap(~Parameter, scales = 'free') +
  scale_fill_viridis(discrete = T) + 
  scale_color_viridis(discrete = T) + 
  theme_bw() + theme(panel.grid = element_blank())

















# Normal Reproduction 
rho_norm <- as.matrix(rho_norm)
rho_norm.data <- reshape2::melt(rho_norm)
colnames(rho_norm.data) <- c('Sample', 'Parameter', 'Value')
rho_norm.sum <- ddply(rho_norm.data, .(Parameter), summarise, median = median(Value, na.rm = T),
                    lCI_90 = quantile(Value, probs = 0.05, na.rm = T),
                    uCI_90 = quantile(Value, probs = 0.95, na.rm = T),
                    lCI_50 = quantile(Value, probs = 0.25, na.rm = T),
                    uCI_50 = quantile(Value, probs = 0.75, na.rm = T))


# Normal Reproduction
rho_rep <- as.matrix(rho_rep)
rho_rep.data <- reshape2::melt(rho_rep)
colnames(rho_rep.data) <- c('Sample', 'Parameter', 'Value')
rho_rep.sum <- ddply(rho_rep.data, .(Parameter), summarise, median = median(Value, na.rm = T),
                      lCI_90 = quantile(Value, probs = 0.05, na.rm = T),
                      uCI_90 = quantile(Value, probs = 0.95, na.rm = T),
                      lCI_50 = quantile(Value, probs = 0.25, na.rm = T),
                      uCI_50 = quantile(Value, probs = 0.75, na.rm = T))


data.list <- dplyr::bind_rows(data.list(list(c(rho_norm.dat, rho_rep.data))))



rhoNorm <- paste('rho[',c(1:Tmax), ']', sep = '')
data.rhoNorm <- subset(rho_norm.sum, Parameter%in% c(rhoNorm))
data.rhoNorm$Model <- c(1)
data.rhoNorm$Year <- c(1:Tmax)

rhoRep <- paste('rho[',c(1:Tmax), ']', sep = '')
data.rhoRep <- subset(rho_rep.sum, Parameter%in% c(rhoRep))
data.rhoRep$Model <- c(2)
data.rhoRep$Year <- c(1:Tmax)


plotNorm <- ggplot(data.rhoNorm) +
              geom_density(aes(x = median, fill = 'red'), alpha = 0.5) + xlim(1,4) +
              scale_fill_viridis(discrete = T) +
              scale_color_viridis(discrete = T) +
              theme_bw() + theme(panel.grid = element_blank())



plotRep <- ggplot(data.rhoRep) +
              geom_density(aes(x = median, fill = 'red'), alpha = 0.5) + xlim(1,5)
              scale_fill_viridis(discrete = T) +
              scale_color_viridis(discrete = T) +
              theme_bw() + theme(panel.grid = element_blank())







data.list <- list(c(Rho_Norm, Rho_Rep))
data.all <- dplyr::bind_cols(data.list)


str(data.all)

mains <- c('rho')

ggplot(subset(data.all, mains %in% median)) +
  geom_density(aes(x = median, color = Model, fill = Model), alpha = 0.5)


Plot <- ggplot(
  geom_density(aes(data = Rho_Norm, x = Year, y = median), col = 'red', alpha = 0.5)) 
  # geom_density(data = Rho_Rep, col = 'blue', alpha = 0.5)
  
Plot <- ggplot(Rho_Norm, aes(x = Year, y = median)) +
  geom_point(aes(col = 'red')) +
  geom_point(data = Rho_Rep, aes(x = Year, y = median, col = 'blue'))

  
  
PlotNorm <- ggplot(data = Rho_Norm) +
  geom_density(aes(x = median, color = Model, fill = Model), alpha = 0.5) 

PlotRep <- ggplot(data = Rho_Rep) +
  geom_density(aes(x = median, color = Model, fill = Model), alpha = 0.5)







  
  
