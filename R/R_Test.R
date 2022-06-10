library(wiqid) # alway run to open wiqid
#GPE = Green Peafowl egg clusth all data from 2015-2022
#BPE = Blue Peafowl egg clusth


GPE<-c(3,7,4,4,4,11,3,5,5,10)
BPE<-c(9,4,6,2,6,9,7,3,6,3,6)


###Comparison of Peafowl cultch size (Green VS. Blue)
# ===================================================

#Data for clutch sizes of Green Peafowl (GPE) and Blue Peafowl (BPE) 
GPE<-c(3,7,4,4,4,11,3,5,5,10)
BPE<-c(9,4,6,2,6,9,7,3,6,3,6)

mean(GPE)
mean(BPE)

GPE.bar <- mean(GPE)
BPE.bar <- mean(BPE)
sd(GPE)
sd(BPE)

# 1. Calculate the sample means and the difference in sample means
mean(GPE)
mean(BPE)
mean(GPE)-mean(BPE)

# 2. Calculate the posterior distribution of the mean for each group
#      using the default priors and plot them.
library(wiqid)
( gp <- Bnormal(GPE) )
plot(gp)
( bp <- Bnormal(BPE) )
plot(bp)

par(mfrow=c(2,1)) # m=multiframe, to row in one collum # to display graphic in two rows and one colimm
plot(gp, xlim=c(0, 12), main="Green Peafowl")
plot(bp, xlim=c(0, 12), main="Blue Peafowl")
par(mfrow=c(1,1)) # to make before display to go back to deafault display setting


# 3. Calculate the posterior distribution for the difference in means
#      and its 95% CrI

head(gp,10)
dif<- gp$mu-bp$mu
plotPost(dif, compVal=2) #comp comparision with 0 or any number in this case more than 2 mm
mean(dif>2)






###Compare informative with uninformative priors
# ==============================================

library(wiqid) # alway run to open wiqid
#GPE = Green Peafowl egg clusth
#Gpr prior: data from 2015 (3 eggs), 2020 (7 eggs), 2021 (4,4,4,11,3)
#GPE 2022: 3,5,5,10

Gpr<-c(3,7,4,4,4,11)
GPE<-c(3,5,5,10)

mean(Gpr)
Gpr.bar <- mean(GPE)
( Gpr.bar <- mean(Gpr) ) # use () cover x.bar <- mean(x) to show them againe 
( sd.Gpr <- sd(Gpr) )   # Should be the same values that you got in the spreads

mean(GPE)
GPE.bar <- mean(GPE)
( GPE.bar <- mean(GPE) ) # use () cover x.bar <- mean(x) to show them againe 
( sd.GPE <- sd(GPE) )   # Should be the same values that you got in the spreadsheet

# Bayesian credible interval
# Uninformative priors i.e. flat priors
( gp1 <- Bnormal(GPE) )
plot(gp1)
densityPlot(gp1)
plot(gp1, "sigma")

# What has Bnormal produced? # data frame of MU (estimate population mean) and sigma (spread of squirrel weith)
head(gp1)
nrow(gp1) # nrow = n = number + row 
colMeans(gp1) # expected value for MU and sigma (calculate the mean for every colum

plot(gp1$mu[1:200], type='l') # trace plot display MU the first 200
plot(gp1$mu[1:3000], type='l') # plot first MU 3000
plot(gp1$mu[1:20], type='l')
plot(gp1$mu[1:50], type='l')
head(gp1$mu[1:20])
abline(h=6,col='green') #population mean
abline(h=GPE.bar,col='blue') #sample mean
( gp1 <- Bnormal(GPE) ) 


# Put in YOUR prior for 'mu':
?Bnormal  ## Bayesian modelling of a normal distribution 
( gp2 <- Bnormal(Gpr)) # priors from our data above
plot(gp2)
plot(gp2, "sigma")


# Compare uninformative vs informative priors: 

par(mfrow=c(2,1)) # m=multiframe, to row in one collum # to display graphic in two rows and one colimm

plot(gp1, xlim=c(1, 12), main="Flat priors")
plot(gp2, xlim=c(1, 12), main="Prior for mu ~ Normal(5.596, 1.0236)")
par(mfrow=c(1,1)) # to make before display to go back to deafault display setting