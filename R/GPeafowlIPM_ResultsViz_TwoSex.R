# Visualization

library(coda)
library(ggplot2)
library(reshape2)
library(viridis)
library(tidyr)
library(plyr)
library(patchwork)

# Set up Age and time 

Amax <- 4
Tmax <- 20


# Set up
## Load posterior data

gpipm <- readRDS('R/PeafowlIPM_TwoSex_TestRun.rds')
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

NBrF1 <- paste('NBreedF[', c(1),',',' ', c(1:Tmax), ']', sep = '') # Juvenile
NBrF2 <- paste('NBreedF[', c(2),',',' ', c(1:Tmax), ']', sep = '') # Yearling
NBrF3 <- paste('NBreedF[', c(3),',',' ', c(1:Tmax), ']', sep = '') # 2 Years
NBrF4 <- paste('NBreedF[', c(4),',',' ', c(1:Tmax), ']', sep = '') # 3 Years

NBrM1 <- paste('NBreedM[', c(1),',',' ', c(1:Tmax), ']', sep = '') # Juvenile
NBrM2 <- paste('NBreedM[', c(2),',',' ', c(1:Tmax), ']', sep = '') # Yearling
NBrM3 <- paste('NBreedM[', c(3),',',' ', c(1:Tmax), ']', sep = '') # 2 Years
NBrM4 <- paste('NBreedM[', c(4),',',' ', c(1:Tmax), ']', sep = '') # 3 Years

## Subset data 

data.NBRF <- subset(data.sum, parameter%in% c(NBrF1,NBrF2,NBrF3,NBrF4))

data.NBRM <- subset(data.sum, parameter%in% c(NBrM1,NBrM2,NBrM3,NBrM4))
  
## Add indexT ans age class time

data.NBRF$indexT <- (1:Amax)
data.NBRF$gen <- rep("Female")
data.NBRF$Year <- rep(1:Tmax, each = Amax)
data.NBRF <- data.NBRF[order(data.NBRF$Year),]

data.NBRM$indexT <- (1:Amax)
data.NBRM$gen <- rep("Male")
data.NBRM$Year <- rep(1:Tmax, each = Amax)
data.NBRM <- data.NBRM[order(data.NBRM$Year),]

## Label age class 

BrF <- within(data.NBRF, {
  Age <- factor(indexT,labels = c("Chick","Yearling","2 Years","3 Years"))
})
BrM <- within(data.NBRM, {
  Age <- factor(indexT,labels = c("Chick","Yearling","2 Years","3 Years"))
})         

## Plot - Estimate with 95% CI
  
PBRF <- ggplot(BrF, aes(x = Year, y = median)) +
  geom_line(aes(color = Age),size = 1) +
  geom_point(aes(color = Age), size = 3) +
  geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90, fill = Age), alpha = 0.3) +
  facet_grid(rows = vars(Age), scales = "free_y") +
  theme_bw() +
  labs(
    title = "Prediction number of 4 female age classes in 20 years",
    subtitle = "Breeding Season",
    # caption =
    x = "Year", 
    y = "Number",
    tag = "Figure 1",
    # colour = 
  )
PBRM <- ggplot(BrM, aes(x = Year, y = median)) +
  geom_line(aes(color = Age),size = 1) +
  geom_point(aes(color = Age), size = 3) +
  geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90, fill = Age), alpha = 0.3) +
  facet_grid(rows = vars(Age), scales = "free_y") +
  theme_bw() +
  labs(
    title = "Prediction number of 4 male age classes in 20 years",
    subtitle = "Non-breeding Season",
    # caption =
    x = "Year", 
    y = "Number",
    tag = "Figure 2",
    # colour = 
  )

PBRF
PBRM

# Aggregate !!!!!! Not function yet

# Plote.Br <- PBRF + PBRM + guide_area() + plot_layout(ncol = 2)




# ## Set data for Non-Breeding

NNoF1 <- paste('NNonF[', c(1),',',' ', c(1:Tmax), ']', sep = '') # Chick
NNoF2 <- paste('NNonF[', c(2),',',' ', c(1:Tmax), ']', sep = '') # Yearling
NNoF3 <- paste('NNonF[', c(3),',',' ', c(1:Tmax), ']', sep = '') # 2 Years
NNoF4 <- paste('NNonF[', c(4),',',' ', c(1:Tmax), ']', sep = '') # 3 Years

NNoM1 <- paste('NNonM[', c(1),',',' ', c(1:Tmax), ']', sep = '') # Chick
NNoM2 <- paste('NNonM[', c(2),',',' ', c(1:Tmax), ']', sep = '') # Yearling
NNoM3 <- paste('NNonM[', c(3),',',' ', c(1:Tmax), ']', sep = '') # 2 Years
NNoM4 <- paste('NNonM[', c(4),',',' ', c(1:Tmax), ']', sep = '') # 3 Years

# ## Subset data 

data.NNOF <- subset(data.sum, parameter%in% c(NNoF1,NNoF2,NNoF3,NNoF4))

data.NNOM <- subset(data.sum, parameter%in% c(NNoM1,NNoM2,NNoM3,NNoM4))

# ## Add indexT ans age class time

data.NNOF$indexT <- (1:Amax)
data.NNOF$gen <- rep("Female")
data.NNOF$Year <- rep(1:Tmax, each = Amax)
data.NNOF <- data.NNOF[order(data.NNOF$Year),]

data.NNOM$indexT <- (1:Amax)
data.NNOM$gen <- rep("Male")
data.NNOM$Year <- rep(1:Tmax, each = Amax)
data.NNOM <- data.NNOM[order(data.NNOM$Year),]


## Label age class 

NoF <- within(data.NNOF, {
  Age <- factor(indexT,labels = c("Chick","Yearling","2 Years","3 Years"))
})
NoM <- within(data.NNOM, {
  Age <- factor(indexT,labels = c("Chick","Yearling","2 Years","3 Years"))
})    

## Plot - Estimate with 95% CI

PNOF <- ggplot(NoF, aes(x = Year, y = median)) +
  geom_line(aes(color = Age),size = 1) +
  geom_point(aes(color = Age), size = 3) +
  geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90, fill = Age), alpha = 0.3) +
  facet_grid(rows = vars(Age), scales = "free_y") +
  theme_bw() +
  labs(
    title = "Prediction number of 4 female age classes in 20 years",
    subtitle = "Breeding Season",
    # caption =
    x = "Year", 
    y = "Number",
    tag = "Figure 1",
    # colour = 
  )
PNOM <- ggplot(NoM, aes(x = Year, y = median)) +
  geom_line(aes(color = Age),size = 1) +
  geom_point(aes(color = Age), size = 3) +
  geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90, fill = Age), alpha = 0.3) +
  facet_grid(rows = vars(Age), scales = "free_y") +
  theme_bw() +
  labs(
    title = "Prediction number of 4 male age classes in 20 years",
    subtitle = "Non-breeding Season",
    # caption =
    x = "Year", 
    y = "Number",
    tag = "Figure 2",
    # colour = 
  )

PNOF
PNOM 

