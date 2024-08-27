# Visualization

library(coda)
library(ggplot2)
library(reshape2)
library(viridis)
library(tidyr)
library(plyr)


# Set up Age and time 

Amax <- 4
Tmax <- 20


## Load posterior data

gpipm <- readRDS('R/PeafowlIPM_TestRun.rds')
gpipm

## Re-arrange data

out.mat <- as.matrix(gpipm)
data <- melt(out.mat)
colnames(data) <- c('index', 'parameter', 'value')

## Summarise posterior into median and 50% and 90% CI

data.sum <- ddply(data, .(parameter), summarise, median = median(value, na.rm = T),
                  lCI_90 = quantile(value, probs = 0.05, na.rm = T),
                  uCI_90 = quantile(value, probs = 0.95, na.rm = T),
                  lCI_50 = quantile(value, probs = 0.25, na.rm = T),
                  uCI_50 = quantile(value, probs = 0.75, na.rm = T))

data.sum


# Other parameter
# popN <- paste('s_BN[', c(1:16), ']', sep = '')


## Set data for Breeding

NBr1 <- paste('NBreed[', c(1),',',' ', c(1:Tmax), ']', sep = '') # Juvenile
NBr2 <- paste('NBreed[', c(2),',',' ', c(1:Tmax), ']', sep = '') # Yearling
NBr3 <- paste('NBreed[', c(3),',',' ', c(1:Tmax), ']', sep = '') # 2 Years
NBr4 <- paste('NBreed[', c(4),',',' ', c(1:Tmax), ']', sep = '') # 3 Years

## Subset data 

data.NBR <- subset(data.sum, parameter%in% c(NBr1,NBr2,NBr3,NBr4))

## Add indexT ans age class time

data.NBR$indexT <- (1:Amax)
data.NBR$Year <- rep(1:Tmax, each = Amax)
data.NBR <- data.NBR[order(data.NBR$Year),]


## Plot - Estimate with 95% CI

plot.NBR <- ggplot(data.NBR, aes(x = Year, y = median, group = indexT, color = indexT)) + 
  geom_line(size = 1) + geom_point(size = 3) + 
  geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90), alpha = 0.2) + 
  facet_wrap(~indexT) 
plot.NBR


## Set data for Non-Breeding

NNo1 <- paste('NNon[', c(1),',',' ', c(1:Tmax), ']', sep = '') # Chick
NNo2 <- paste('NNon[', c(2),',',' ', c(1:Tmax), ']', sep = '') # Yearling
NNo3 <- paste('NNon[', c(3),',',' ', c(1:Tmax), ']', sep = '') # 2 Years
NNo4 <- paste('NNon[', c(4),',',' ', c(1:Tmax), ']', sep = '') # 3 Years

## Subset data 

data.NNO <- subset(data.sum, parameter%in% c(NNo1,NNo2,NNo3,NNo4))

## Add indexT ans age class time

data.NNO$indexT <- (1:Amax)
data.NNO$Year <- rep(1:Tmax, each = Amax)
data.NNO <- data.NNO[order(data.NNO$Year),]

## Plot - Estimate with 95% CI

plot.NNO <- ggplot(data.NNO, aes(x = Year, y = median, group = indexT, color = indexT)) + 
  geom_line(size = 1) + geom_point(size = 3) + 
  geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90), alpha = 0.2) + 
  facet_wrap(~indexT) 
plot.NNO





