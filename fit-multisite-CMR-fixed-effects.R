setwd("/Users/Brian/Dropbox/Gopher Tortoise Demography/Tortoise-PVA-Conecuh")
## for running on macbook

setwd("C:/Users/bpf0006/Desktop/Tortoise-PVA")
## for running on remote desktop PC

#setwd("~/Dropbox/Gopher Tortoise Demography/Tortoise-IPM-GitHub")

library(imputeTS) # Data interpolation
library(statmod)
library(jagsUI) # to run mcmc in parallel
library(pushoverr)

source("functions.R")

# cmr data ----
## Load mark-recapture data and wrangle into formats for multi-state mark-recapture analysis
captures = read.csv("captures.csv", header=TRUE)

## Subset to six study sites, combine into a DF, convert site to numeric
one = droplevels(subset(captures, captures$Site == "1"))
one = droplevels(subset(one, one$Year != "1995")) # remove 1995, when site was not sampled thoroughly
two = droplevels(subset(captures, captures$Site == "2"))
three = droplevels(subset(captures, captures$Site == "3"))
four = droplevels(subset(captures, captures$Site == "4"))
four = droplevels(subset(four, four$Year != "2020"))  # remove one dead-recovery from 2020
five = droplevels(subset(captures, captures$Site == "5"))
six = droplevels(subset(captures, captures$Site == "6"))

captures = rbind(one,two,three,four,five,six)
captures$Site = as.numeric(as.character(captures$Site))

## Create mark-recapture histories for individuals in each site
for (i in c(1:6)){	# For each i in sites 1-6
  
  # Set up a matrix to save capture histories
  site = droplevels(subset(captures, captures$Site == i))
  inds = sort(as.numeric(as.character(levels(factor(site$Number)))))
  years = sort(as.numeric(as.character(levels(factor(site$Year)))))
  #interval = c(min(years):(min(years)+(max(years)-min(years))))
  interval = c(1991:2020)	
  ### Hard-code the mark-recapture interval from 1991-2020
  ch = matrix(NA, length(inds), length(interval))
  rownames(ch) = inds
  colnames(ch) = interval
  
  # Mark-recapture histories should be state-specific,
  # where individuals observed each year are specified as:
  # juveniles = 1, females = 2, males = 3
  for (j in 1:length(inds)){
    ind = subset(site, site$Number == inds[j])
    for (k in 1:length(interval)){
      yr = subset(ind, ind$Year == interval[k])
      # Specify J = 1, F = 2, M = 3, or NA for unsampled years
      if(interval[k] %in% levels(as.factor(site$Year)) == FALSE) {
        ch[j,k] = NA} else {
          if(length(yr$Year) > 0 && yr$ReproSex[1] == "J") {
            ch[j,k] = 1} else {
              if(length(yr$Year) > 0 && yr$ReproSex[1] == "F") {
                ch[j,k] = 2} else {
                  if(length(yr$Year) > 0 && yr$ReproSex[1] == "M") {
                    ch[j,k] = 3} else {ch[j,k] = 0}}}}
    }
  }
  
  # Assign capture histories a unique object by site
  assign(paste0("Site", i), ch) 
}

# bind all capture hx together
CH = rbind(Site1, Site2, Site3, Site4, Site5, Site6)
sites = list(Site1, Site2, Site3, Site4, Site5, Site6)


# site ID for each individual
nind.site = sapply(sites, nrow)
site = rep(c(1:6), nind.site)

nyears = dim(CH)[2]
nsites = length(unique(sites))

CH <- CH  # Recoded CH
CH[CH==0] <- 4

get.first <- function(x) min(which(x!=4))  # 4 is "not observed"
f <- apply(CH, 1, get.first) 

CH.NA <- CH
CH.NA[CH.NA==4] <- NA
z.inits <- ms.init.z(CH.NA,NA) 

# change all NAs to 4s (not seen)
CH.test <- CH
CH.test[is.na(CH.test)] <- 4

# survey matrix == 1 if site j was surveyed in year t and 0 otherwise
surveyCMR = matrix(1, nrow = length(sites), ncol = dim(CH)[2])
for(i in 1:length(sites)){
  ch = sites[[i]]
  missingCMR = which(apply(ch, 2, function(x) is.na(sum(x))))
  surveyCMR[i,missingCMR] <- 0
}

# fit model ----

inits <- function(){list(mean.phiJ=runif(6,0,1), 
                         mean.phiF=runif(6,0,1), 
                         mean.phiM=runif(6,0,1),
                         mean.tauJA=runif(6,0,1), 
                         mean.pFemale=runif(6,0,1),
                         overall.mean.pJ = runif(1, 0, 1),
                         overall.mean.pF = runif(1, 0, 1),
                         overall.mean.pM = runif(1, 0, 1),
                         z = z.inits)} 

# R = prior for variance-covariance matrices (large values on diag = large variances)
# K&S recommend a sensitivity analysis for choice of values in R matrix
R = matrix(0, nrow = nsites, ncol = nsites)
diag(R) <- 1

jags.data <- list(y = CH.test, 
                  f = f, 
                  nyears = nyears, 
                  nsites = 6,
                  nind = dim(CH.test)[1],
                  surveyCMR = surveyCMR,
                  site = site,
                  R = R,
                  zeros = rep(0, nsites)
)

parameters <- c("mean.phiJ", "mean.phiF", "mean.phiM",
                "mean.tauJA", "mean.pFemale",
                "mean.pJ", "mean.pF", "mean.pM", 
                "sigma2.y", "sigma.J", "sigma.F", "sigma.M",
                "sigma.tauJA", "sigma.pFemale",
                "beta.site.J", "beta.site.F", "beta.site.M",
                "beta.site.tauF", "beta.site.tauM")

ni <- 150000
nt <- 10
nb <- 100000
nc <- 3
na <- 50000 


msCMR <- jags(jags.data, inits, parameters, "multisite-CMR_fixed-site-effects.jags", 
              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt = na, parallel = TRUE)


#library(beepr)
#beep(sound=8)

pushover("Multisite CMR done - Carribean", 
         user="ufxwvtmc4jby3zgsetfsya6n8wzmkj", 
         app="asbrmzwsprh7r4isr6pra5itgjeeng")

summary(msCMR)

# Save the RDS file for use in population projection (PVA) analysis

saveRDS(msCMR, "multisite-CMR_fixed-site-effects.rds")



#### Looking at the results
msCMR = readRDS("multisite-CMR_fixed-site-effects.rds")

# Look at results
options(scipen=999)
msCMR$summary

# Traceplots
par(mfrow = c(3,3))
plot(msCMR)

# Check Rhat
par(mfrow=c(1,1))
hist(msCMR$summary[,8])
length(msCMR$summary[,8])
length(which(msCMR$summary[,8]< 1.1))

# Explore mean, LCL, UCL for parameter estimates at each site
msCMR$summary[1:6,c(1,3,7)]   #phiJ 
msCMR$summary[7:12,c(1,3,7)]  #phiF 
msCMR$summary[13:18,c(1,3,7)] #phiM

msCMR$summary[19:24,c(1,3,7)]    #tauJA
msCMR$summary[25:30,c(1,3,7)]    #pFemale
1-msCMR$summary[25:30,c(1,3,7)]  #pMale

# Create tauJF and tauJM by multiplying tauJA*pFemale and tauJA*(1-pFemale)
str(msCMR)
tauJF = msCMR$sims.list$mean.tauJA*msCMR$sims.list$mean.pFemale
tauJM = msCMR$sims.list$mean.tauJA*(1-msCMR$sims.list$mean.pFemale)

TauJF = matrix(NA,6,3,dimnames=list(c(1:6),c("Mean","LCL","UCL")))
for (i in 1:6){
  site = tauJF[,i]
  TauJF[i,1] = mean(site)
  TauJF[i,2] = apply(data.frame(site),2,quantile,probs=c(0.025))
  TauJF[i,3] = apply(data.frame(site),2,quantile,probs=c(0.975))
}

TauJM = matrix(NA,6,3,dimnames=list(c(1:6),c("Mean","LCL","UCL")))
for (i in 1:6){
  site = tauJM[,i]
  TauJM[i,1] = mean(site)
  TauJM[i,2] = apply(data.frame(site),2,quantile,probs=c(0.025))
  TauJM[i,3] = apply(data.frame(site),2,quantile,probs=c(0.975))
}












