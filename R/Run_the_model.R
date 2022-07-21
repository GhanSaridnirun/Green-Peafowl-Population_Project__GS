library(bayesplot)
library(rjags)
library(bayestestR)
library(tidyverse)


month <- factor(c("Nov19","Dec19",
                  "Jan20","Feb20","Mar20","Apr20","May20","Jun20","Jul20","Aug20",
                  "Sep20","Oct20","Nov20","Dec20","Jan21","Feb21","Mar21","Apr21",
                  "May21","Jun21","Jul21","Aug21","Sep21","Oct21","Nov21","Dec21"
),ordered = TRUE)

season <- factor(c("Breeding","Breeding","Breeding","Breeding","Breeding",
                   "Breeding","Non","Non","Non","Non","Non","Non","Breeding","Breeding",
                   "Breeding","Breeding","Breeding","Breeding","Non","Non","Non","Non",
                   "Non","Non","Breeding","Breeding"
),ordered = TRUE)


#chick
ch <- c(0,0,0,0,0,0,36,65,147,152,112,41,0,0,0,0,0,0,50,50,81,113,24,16,0,0
)
#juvenile
ju <- c(22,35,61,62,18,21,36,2,1,0,0,0,26,18,23,63,25,44,19,6,5,0,0,0,60,0	
)
#breeder
br <- c(9,23,31,25,12,10,35,34,83,78,56,25,14,9,11,20,13,27,
        29,28,47,57,14,8,0,0
)


#Survival rates & fecundity
ch.s <- 0.3        #Chick survival
ju.s <- 0.5        #Juvenile survival
br.s <- 0.7        #Breeder survival
br.f <- 1.9        #Average Number of chicks/Breeder


#Population Growth rate (shall be rho??)
#Use totall number of Pool population in Nov 2019(308) - Oct 2020(501) = 809
#minor with number of adult females in Nov 2020(292) - Oct 2021(377) = 669
#overall.growth.rate <- N[t+1] / N[t]
overall.growth.rate <- 669 / 809

rho <- overall.growth.rate
rho


jags.data=list(ch=ch, ju=ju, br=br, ch.s=ch.s, ju.s=ju.s, br.s=br.s, br.f=br.f, rho=rho)


#jag.m <- jags.model( file = "modelgp.txt", data, n.chains=3, n.adapt=2000 )


#jags.model
cat (file="modelgp.txt", "
model {
#T = 20 years following GP age

#Priors
#Survival & Fecundity
for (ch.s in 1:1) {
 for (ju.s in 1:1) {
  for (br.s in 1:1) {
    for (br.f in 1:1){
     for (rho in 1:1) {
     s[ch.s, ju.s, br.s, br.f, rho] ~ dnorm(1, 1)
    }
   }
  }
 }
}


#Process model ####### First year
for (t in 1: 20) {
# Process model: in non breeding 
# Total number of chicks
F[t] ~ dpois(br[t] * rho[t])
# Allocate chicks to sex (pick only female)
ch[t] ~ dbin(gamma, F[t]) # or just (F[t]/2)
ch[t] <- (br.s * br[t]) * br.f
}

# Process model: in breeding
for (t in 1: 20) {
ju[t+1] <- ch.s * ch[t] 
br[t+1] <- br.s * br[t]
}


##### Protocal for second year and so on...........

###### t+1 ################### 
### Population got more breeder from previous year (ju,t+1)
### Due to (ju,t+1) is fully mature they have same survival with breeder

for (t in 1: 20) {
# Process model: in non breeding 
# Total number of chicks
F[t+1] ~ dpois(br[t+1] * rho[t+1])
# Allocate chicks to sex (pick only female)
ch[t+1] ~ dbin(gamma, F[t+1]) # or just (F[t]/2)
ch[t+1] <- ((br.s * br[t+1]) + (ju.s * ju[t+1])) * br.f
}

# Process model: in breeding
for (t in 1: 20) {
ju[t+2] <- ch.s * ch[t+1]
br[t+2] <- br.s * (br[t+1] + ju[t+1]) 
}


##### t+2 ###################
for (t in 1: 20) {
# Process model: in non breeding 
# Total number of chicks
F[t+2] ~ dpois(br[t+2] * rho[t+2])
# Allocate chicks to sex (pick only female)
ch[t+2] ~ dbin(gamma, F[t+2]) # or just (F[t]/2)
ch[t+2] <- (br.s * br[t+2]) * br.f
}

# Process model: in breeding
for (t in 1: 20) {
ju[t+3] <- ch.s * ch[t+2]
br[t+3] <- br.s * br[t+2]
}


##### t+3 ###################
for (t in 1: 20) {
# Process model: in non breeding 
# Total number of chicks
F[t+3] ~ dpois(br[t+3] * rho[t+3])
# Allocate chicks to sex (pick only female)
ch[t+3] ~ dbin(gamma, F[t+3]) # or just (F[t]/2)
ch[t+3] <- ((br.s * br[t+3]) + (ju.s * ju[t+3])) * br.f
}

for (t in 1: 20) {
# Process model: in breeding
ju[t+4] <- ch.s * ch[t+3]
# number of br+ju might needed to /2 due to ratio for breeder:un breeder = 50:50
# mean un breeder might refer to failed breeder + immature females
br[t+4] <- br.s * (br[t+3] + ju[t+3]) 
}


# Observation models
for (t in 1: 20) {
ch[t] ~ dpois(Nch[t]) 
ju[t] ~ dpois(Nju[t]) 
br[t] ~ dpois(Nbr[t])
}
")


#Model for initial state 
#Use the number from population Nov-19 - Apr 20 and 
#generate random value based on the population number.
inits <- function (){
  Nall <- array(NA, dim=c(3, 3, 20))
  Nall[1,,]<- round(runif(4, 100, 250))
  s <- array(runif(3, 0.3, 0.7), dim=c(3, 3, 6))
  list(Nall=Nall, s=s)
}


# Parameters monitored
parameters <- c("s", "mean.rho", "sigma.rho", "rho", "ch", "ju", "br")


# MCMC settings
ni <- 100000; nb <- 20000; nc <- 3; nt <- 20; na <- 10000


# Call JAGS (ART 8 min), check convergence and summarize posteriors
set.seed(101) # Reproducible results without crashing
out1 <- jagsUI::jags(jags.data, inits, parameters, "modelgp.txt", n.iter=ni, n.burnin=nb
             , n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE)

traceplot(out1)
