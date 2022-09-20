library(coda)
library(ggplot2)
library(reshape2)
library(viridis)
library(tidyr)
library(plyr)



gpipm <- readRDS('PeafowlIPM_TestRun.rds')
gpipm

str(gpipm)

out.mat <- as.matrix(gpipm)
out.mat
data <- melt(out.mat)
data

# vr.params <- c('NBreed[1,1]','NBreed[1,2]','NBreed[1,3]','NBreed[1,4]',
#                'NNon[1,1]','NNon[1,2]','NNon[1,3]','NNon[1,4]',
#                'Fec[1]','Fec[2]','Fec[3]','Fec[4]',
#                'p[1]','p[2]','p[3]','p[4]',
#                's_NB[1]','s_NB[2]','s_NB[3]','s_NB[4]',
#                's_BN[1]','s_BN[2]','s_BN[3]','s_BN[4]')


colnames(data) <- c('index', 'parameter', 'value')
data

str(data)

data.sum <- ddply(data, .(parameter), summarise, median = median(value, na.rm = T),
                  lCI_90 = quantile(value, probs = 0.05, na.rm = T),
                  uCI_90 = quantile(value, probs = 0.95, na.rm = T),
                  lCI_50 = quantile(value, probs = 0.25, na.rm = T),
                  uCI_50 = quantile(value, probs = 0.75, na.rm = T))

data.sum


plot.data <- data.frame(data.sum)
plot.data

str(plot.data)


ggplot(plot.data, aes(x = median, y = parameter))  +
  geom_point() +
  geom_point(data = plot.data, aes(x = median), colour = 'blue', size = 1)
