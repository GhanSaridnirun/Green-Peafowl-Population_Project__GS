library(bayesplot)
library(rjags)
library(bayestestR)
library(tidyverse)
#GPF data
y<-c(3,7,4,4,4,11,3,5,5,10)
data=list(y=y,
         N=length(y))
jags.m <- jags.model( file = "JAGS_mod.txt", data, n.chains=1, n.adapt=500 )

params=c("mu", "sigma")
samps <- coda.samples( jags.m, params, n.iter=10000 )
posterior<-as.mcmc(samps)
describe_posterior(posterior)
hdi(posterior)
dens=estimate_density(posterior)

dens %>% filter(Parameter=="mu") %>% 
  ggplot(mapping = aes(x = x, y = y)) +
  theme_classic() +
  geom_area(fill = "lightblue")+
  geom_vline(xintercept = mean(posterior[,1]), lty=2)+
  labs(x="mu", y="")+
  annotate("text", x = mean(posterior[,1]+2), y = 0.4, label = paste0("mean= ",round(mean(posterior[,1]),2)))

dens %>% filter(Parameter=="sigma") %>% 
  ggplot(mapping = aes(x = x, y = y)) +
  theme_classic() +
  geom_area(fill = "lightblue")+
  geom_vline(xintercept = mean(posterior[,2]), lty=2)+
  labs(x="sigma", y="")+
  annotate("text", x = mean(posterior[,2]+1), y = 0.4, label = paste0("mean= ",round(mean(posterior[,2]),2)))
