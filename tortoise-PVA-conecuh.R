############################################################################
############################################################################
##########    Modeling population viability of gopher tortoises   ##########
##########      in Conecuh National Forest, Alabama using         ##########
##########            population projection models                ##########
##########               B Folt, October 2020                     ##########
############################################################################
############################################################################

# Clear workspace 
rm(list=ls())

# Set working directory
setwd("/Users/Brian/Dropbox/Gopher Tortoise Demography/Tortoise-PVA-Conecuh")

# Objective: use matrix population projection models to estimate
# population growth and extinction risk for six tortoise populations
# at Conecuh National Forest


########## Section 1)
########## Summarize mark-recapture data and 
########## describe changes in population structure through time

### Summarize mark-recapture results at each site
captures = read.csv("captures.csv", header=TRUE)

res = matrix(NA,6,5, dimnames=
        list(c("Site1","Site2","Site3","Site4","Site5","Site6"),
        c("Site","Captures","Individuals","Recaps","MeanRecaps")))
for (i in 1:6){
  site = droplevels(subset(captures, captures$Site == i))
  caps = length(site[,1])
  inds = length(table(as.numeric(site$Number)))
  recaps = length(subset(site, site$Status == "R")[,1])
  mrecaps = mean(table(as.numeric(site$Number)))
  res[i,] = c(i,caps,inds,recaps,mrecaps)
}
res
sum(res[,2]) # total captures
sum(res[,3]) # total individuals


#install.packages("clipr")
library(clipr)
write_clip(res)
# Paste 'res' into manuscript prep Excel file to make Table 2 


### Summarize population structure at sites through time
res = matrix(NA,1,7)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

res1 = matrix(NA,1,7)
colnames(res1) = c("Era","Site","Number","Sex","AgeSex","ReproSex","CLmax")

for (h in c(1,2)){
  if(h == 1){caps = droplevels(subset(captures, captures$Year < 2004))}
      else(caps = droplevels(subset(captures, captures$Year > 2004)))
  for (i in 1:6){
      site = droplevels(subset(caps, caps$Site == i))
  for (j in levels(as.factor(as.numeric(site$Number)))){
      ind = droplevels(subset(site, site$Number == j))
      if(h==1){res[1,1] = c("1991-2003")}else(res[1,1] = c("2013-2020"))
      res[1,2:7] = c(i,j,as.character(Mode(ind$Sex)),as.character(Mode(ind$AgeSex)),
          as.character(Mode(ind$ReproSex)),mean(na.omit(ind$Clmax)))
      res1 = rbind(res1,res)
    }
  }
}

result=data.frame(res1[-c(1),])
head(result)
str(result) 
# all vars are categorical, need to reformat continuous as such

result$Site <- as.numeric(as.character(result$Site))
result$Number <- as.numeric(as.character(result$Number))
result$CLmax <- as.numeric(as.character(result$CLmax))
str(result)

# Plot histograms of population structure
library(ggplot2)
library(easyGgplot2)
# had to first update R and Rstudio, and then manually installed using: 
# install.packages("devtools")
# devtools::install_github("kassambara/easyGgplot2")


plot = ggplot2.histogram(data=result, xName='CLmax',
  groupName='ReproSex', legendPosition="top",
  alpha=0.5, position="stack", faceting=TRUE,
  facetingVarNames=c("Site","Era"), 
  facetingRect=list(lineType="solid",lineColor="black"),
  backgroundColor="white")
plot = ggplot2.customize(plot, xtitle="        Carapace length (cm)",
  ytitle="Number observed", legendTitle="Age class")
(plot = plot + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()))


# Proportion of males, females, and juveniles
early = subset(result, result$Era == "1991-2003")
popstr1 = matrix(NA,6,8)
colnames(popstr1) = c("Site","Juveniles","Females","Males","Adults","Total","M-F-ChSq","Ad-J-ChiSq")
for (i in 1:6){
  site = subset(early, early$Site == i)
  popstr1[i,] = c(i,as.numeric(length(subset(site, ReproSex=="J")[,6])),
                  as.numeric(length(subset(site, ReproSex=="F")[,6])),
                  as.numeric(length(subset(site, ReproSex=="M")[,6])),
                  as.numeric(length(subset(site, ReproSex=="F")[,6])+length(subset(site, ReproSex=="M")[,6])),
                  as.numeric(length(site[,6])), 
                  as.numeric(chisq.test(c(length(subset(site, ReproSex=="F")[,6]),length(subset(site, ReproSex=="M")[,6])))$p.value),
                  as.numeric(chisq.test(c(length(subset(site, ReproSex=="F")[,6])+length(subset(site, ReproSex=="M")[,6]),
                                          length(subset(site, ReproSex=="J")[,6])))$p.value))
}

recent = subset(result, result$Era == "2013-2020")
popstr2 = matrix(NA,6,8)
colnames(popstr2) = c("Site","Juveniles","Females","Males","Adults","Total","M-F-ChSq","Ad-J-ChiSq")
for (i in 1:6){
  site = subset(recent, recent$Site == i)
  popstr2[i,] = c(i,as.numeric(length(subset(site, ReproSex=="J")[,6])),
                  as.numeric(length(subset(site, ReproSex=="F")[,6])),
                  as.numeric(length(subset(site, ReproSex=="M")[,6])),
                  as.numeric(length(subset(site, ReproSex=="F")[,6])+length(subset(site, ReproSex=="M")[,6])),
                  as.numeric(length(site[,6])), 
                  as.numeric(chisq.test(c(length(subset(site, ReproSex=="F")[,6]),length(subset(site, ReproSex=="M")[,6])))$p.value),
                  as.numeric(chisq.test(c(length(subset(site, ReproSex=="F")[,6])+length(subset(site, ReproSex=="M")[,6]),
                               length(subset(site, ReproSex=="J")[,6])))$p.value))
}

popstr = rbind(popstr1,popstr2)
Era = c("Old","Old","Old","Old","Old","Old",
        "New","New","New","New","New","New")
POPSTR = cbind(Era,popstr)
# Dataframe of population structure among sites through time 

pop1 = data.frame(t(popstr[c(1,7),c(2:4)]))
colnames(pop1) = c("Old","New")

pop2 = data.frame(t(popstr[c(2,8),c(2:4)]))
colnames(pop2) = c("Old","New")

pop3 = data.frame(t(popstr[c(3,9),c(2:4)]))
colnames(pop3) = c("Old","New")

pop4 = data.frame(t(popstr[c(4,10),c(2:4)]))
colnames(pop4) = c("Old","New")

pop5 = data.frame(t(popstr[c(5,11),c(2:4)]))
colnames(pop5) = c("Old","New")

pop6 = data.frame(t(popstr[c(6,12),c(2:4)]))
colnames(pop6) = c("Old","New")

(test = fisher.test(pop1))
(test = fisher.test(pop2))
(test = fisher.test(pop3))
(test = fisher.test(pop4))
(test = fisher.test(pop5))
(test = fisher.test(pop6))

### Manually annotate facets to indicate signif pop. structure
### add * to Sites 5+6 in left column for biased sex ratio
### add ** to Site 5 in right column to indicate fewer juvenile:adult than 1:1
### add cross to Site 5 to indicate shift in structure through time


## Subset to six study sites, combine into a DF, convert site to numeric
one = droplevels(subset(captures, captures$Site == "1"))
two = droplevels(subset(captures, captures$Site == "2"))
three = droplevels(subset(captures, captures$Site == "3"))
four = droplevels(subset(captures, captures$Site == "4"))
five = droplevels(subset(captures, captures$Site == "5"))
six = droplevels(subset(captures, captures$Site == "6"))

# Remove one individual from Site 1 in 1995, when the 
# site was not sampled thoroughly
one = droplevels(subset(one, one$Year != "1995"))

# Remove a dead-recovery from Site 4 in 2020
four = droplevels(subset(four, four$Year != "2020"))

# Bind capture histories together
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

# Bind all capture histories together
CH = rbind(Site1, Site2, Site3, Site4, Site5, Site6)

# Indexed list of each site
sites = list(Site1, Site2, Site3, Site4, Site5, Site6)

# Survey matrix, where
# values == 1 if site j was surveyed in year t,
# and 0 otherwise
surveys = matrix(1, nrow = 6, ncol = 30, dimnames=list(c(1:6),c(1991:2020)))
for(i in 1:6){
  ch = sites[[i]]
  missingCMR = which(apply(ch, 2, function(x) is.na(sum(x))))
  surveys[i,missingCMR] <- 0
}

#library(clipr)
write_clip(surveys)
  # copy and paste table into manuscript prep Excel to 
  # tabulate years with CMR data in
  # Supplementary Table 1



########## Section 2)
########## Graphically analyze results from survival analysis

# Load the results of survival analysis
conecuh = readRDS("multisite-CMR_fixed-site-effects.rds")
options(scipen=999)
conecuh$summary

## Copy survival estimates and CI into a table for supplementary document
apparentsurvival = conecuh$summary[1:18,c(1,3,7)]

#library(clipr)
write_clip(apparentsurvival)
# copy and paste table into manuscript prep Excel file to 
# summarize apparent survival results by site in Supplementary Table 2


### Plot population vital rates
par(mfrow=c(1,1)) 
graphics.off()

for (i in 1:6){
mean.phiJ <- mean(conecuh$sims.list$mean.phiJ[,i])
lclPhiJ <- quantile(conecuh$sims.list$mean.phiJ[,i], probs = 0.025)
uclPhiJ <- quantile(conecuh$sims.list$mean.phiJ[,i], probs = 0.975)

mean.phiF <- mean(conecuh$sims.list$mean.phiF[,i])
lclPhiF <- quantile(conecuh$sims.list$mean.phiF[,i], probs = 0.025)
uclPhiF <- quantile(conecuh$sims.list$mean.phiF[,i], probs = 0.975)

mean.phiM <- mean(conecuh$sims.list$mean.phiM[,i])
lclPhiM <- quantile(conecuh$sims.list$mean.phiM[,i], probs = 0.025)
uclPhiM <- quantile(conecuh$sims.list$mean.phiM[,i], probs = 0.975)

phi = c(conecuh$sims.list$mean.phiJ[,i],conecuh$sims.list$mean.phiF[,i], conecuh$sims.list$mean.phiM[,i])
stages = as.factor(c(c(rep("Juveniles",length(phi)/3), c(rep("Females", length(phi)/3)), c(rep("Males", length(phi)/3)))))
site = rep(i,length(stages))
survival = data.frame(site,stages,phi)
colnames(survival) = c("Site","Stage","Phi")
assign(paste0("survival",i), survival)
}
PHI = rbind(survival1,survival2,survival3,survival4,survival5,survival6)
PHI$Stage = factor(PHI$Stage, levels=c('Juveniles','Females','Males'))

library(ggplot2)
library(magrittr)
library(ggpubr)
# had to update R and Rstudio in order to install this using: 
# install.packages("devtools")
# devtools::install_github("kassambara/ggpubr")


# Violin plot of survival by stage for one site;
# e.g., Site1 = survival1
ggviolin(survival1, x = "Stage", y = "Phi", fill = "Stage", 
         palette = c("pink","lightgreen","lightblue"), xlab=FALSE,
         ylim=c(0.4,1), lty=2, ylab="Survival", add="mean_sd")

# Multi-panel violin plot that describes survival by stage among sites
ggviolin(PHI, x = "Site", y = "Phi", combine=TRUE, facet.by="Stage",
         fill = "Site", font.y=16, xlab="Site", font.x=16,
         #palette = c("pink","lightgreen","lightblue"),
         ylim=c(0.3,1), lty=2, ylab=expression(paste("Apparent survival probability (", phi, ") ")),
         add="mean_sd")


# Create tauJF and tauJM by multiplying tauJA*pFemale and tauJA*(1-pFemale)
str(conecuh)
tauJF = conecuh$sims.list$mean.tauJA*conecuh$sims.list$mean.pFemale
tauJM = conecuh$sims.list$mean.tauJA*(1-conecuh$sims.list$mean.pFemale)

TauJF = matrix(NA,6,4,dimnames=list(c(1:6),c("Mean","SD","LCL","UCL")))
for (i in 1:6){
  site = tauJF[,i]
  TauJF[i,1] = mean(site)
  TauJF[i,2] = sd(site)
  TauJF[i,3] = apply(data.frame(site),2,quantile,probs=c(0.025))
  TauJF[i,4] = apply(data.frame(site),2,quantile,probs=c(0.975))
}

TauJM = matrix(NA,6,4,dimnames=list(c(1:6),c("Mean","SD","LCL","UCL")))
for (i in 1:6){
  site = tauJM[,i]
  TauJM[i,1] = mean(site)
  TauJM[i,2] = sd(site)
  TauJM[i,3] = apply(data.frame(site),2,quantile,probs=c(0.025))
  TauJM[i,4] = apply(data.frame(site),2,quantile,probs=c(0.975))
}

#library(clipr)
write_clip(TauJF)
write_clip(TauJM)
  # Copy and paste transition probability estimates into manuscript prep Excel file
  # to make Supplementary Table 3


########## Section 3)
########## Construct a population projection model

## Define number of replicate simulations and time-period for population project
r = 1000      # Number of simulation replicates
t = 100       # Number of years to project simulation replicates

## Population growth & quasi-extinction estimates
lam = matrix(0,r,t)   # Population growth
Pext = matrix(0,r,t)  # Probability of quasi-extinction

## Summarize demographic estimates for study populations
s=6 # Number of research sites

params = matrix(NA, nrow=s, ncol=10)
dimnames(params) = list(c("Site1","Site2","Site3","Site4","Site5","Site6"), 
                        c("phiH","phiJ","phiF","phiM","gJF","gJM","f","nJ","nF","nM"))
params[,1] <- rep(0.15,6)                 # Hatchling survival
params[,2] <- conecuh$summary[1:6,1]      # Juvenile survival
params[,3] <- conecuh$summary[7:12,1]     # Female survival
params[,4] <- conecuh$summary[13:18,1]    # Male survival
params[,5] <- TauJF[,1]    # Transition probability from j - ad female
params[,6] <- TauJM[,1]    # Transition probability from j - ad male
params[,7] <- rep(5.4,6)                  # Mean fecundity
params[,8] <- c(6,4,4,10,2,6)        # Initial abundance of juveniles
params[,9] <- c(10,7,9,16,5,8)     # Initial abundance of adult females
params[,10] <- c(10,7,9,16,5,8)     # Initial abundance of adult males

paramsSE = matrix(NA, nrow=s, ncol=10)
dimnames(paramsSE) = list(c("Site1","Site2","Site3","Site4","Site5","Site6"), 
                          c("phiHse","phiJse","phiFse","phiMse","gJFse","gJMse","fse","nJse","nFse","nMse"))
paramsSE[,1] <- rep(0.03,6)                 # Hatchling survival SE
paramsSE[,2] <- conecuh$summary[1:6,2]      # Juvenile survival SE
paramsSE[,3] <- conecuh$summary[7:12,2]     # Female survival SE
paramsSE[,4] <- conecuh$summary[13:18,2]    # Male survival SE
paramsSE[,5] <- TauJF[,2]    # Transition probability from j - ad female SD
paramsSE[,6] <- TauJM[,2]    # Transition probability from j - ad male SD
paramsSE[,7] <- rep(0.5,6)                    # Mean fecundity SE
paramsSE[,8] <- c(2,2,2,4,2,2)      # Initial abundance of juveniles SE
paramsSE[,9] <- c(3,2,3,5,2,3)     # Initial abundance of adult females SE
paramsSE[,10] <- c(3,2,3,5,2,3)     # Initial abundance of adult males SE

#### Use for-loops to simulate and projection population demographics across sites
library(lognorm)

## Start a timer!
ptm <- proc.time()

## Population-level loop
for (h in c(1:(s))){    # Subset to research site, i

  ## Sh - survival of hatchlings
  mSh = params[h,1]                 # Mean hatchling survival		
  varSh = paramsSE[h,1]							# Variance of mean
  aSh = mSh*((mSh*(1-mSh)/(varSh^2))-1) 	 	
  bSh = (1-mSh)*((mSh*(1-mSh)/(varSh^2))-1) 
  Sh = matrix(rbeta(r,aSh,bSh),r,1)	# Parametric uncertainty 
  #hist(Shi)
  #SDmShi = matrix(rinvgauss(r,varSh^2,1),r,1)
  AShi = matrix(0,r,1)				# beta distribution shape parameters
  BShi = matrix(0,r,1)				# beta distribution shape parameters
  Sht = matrix(0,r,t)					# annual variation survival

  ## Sj - survival of juveniles 
  #mSj = params[h,2] 		
  #varSj = paramsSE[h,2]		
  #aSj = mSj*((mSj*(1-mSj)/(varSj^2))-1) 	 	
  #bSj = (1-mSj)*((mSj*(1-mSj)/(varSj^2))-1) 
  #Sj = matrix(rbeta(r,aSj,bSj),r,1)			 
  samples = matrix(sample(1:length(conecuh$sims.list$mean.phiJ[,1]),r,replace=FALSE),r)
  Sj = matrix(0,r,1)  
  varSj = matrix(0,r,1)
  # Draw survival and error from empirical results at each site 
  # for this and other parameters
  for (q in 1:length(samples)){   
    Sj[q,1] = conecuh$sims.list$mean.phiJ[samples[q],h]
    varSj[q,1] = conecuh$sims.list$sigma.J[samples[q],h]
    }
  #SDmSji = matrix(rinvgauss(r,varSj^2,1),r,1)
  ASji = matrix(0,r,1)				
  BSji = matrix(0,r,1)				
  Sjt = matrix(0,r,t)					

  ## Sf - survival of females 
  #mSf = params[h,3] 		
  #varSf = paramsSE[h,3]		#
  #aSf = mSf*((mSf*(1-mSf)/(varSf^2))-1) 	 	
  #bSf = (1-mSf)*((mSf*(1-mSf)/(varSf^2))-1) 
  #Sf = matrix(rbeta(r,aSf,bSf),r,1)	
  Sf = matrix(0,r,1)
  varSf = matrix(0,r,1)
  for (q in 1:length(samples)){
    Sf[q,1] = conecuh$sims.list$mean.phiF[samples[q],h]
    varSf[q,1] = conecuh$sims.list$sigma.F[samples[q],h]
  }  
  #SDmSfi = matrix(rinvgauss(r,varSf^2,1),r,1)
  ASfi = matrix(0,r,1)
  BSfi = matrix(0,r,1)
  Sft = matrix(0,r,t)

  ## Sm - survival of males 
  #mSm = params[h,4] 		
  #varSm = paramsSE[h,4]		#
  #aSm = mSm*((mSm*(1-mSm)/(varSm^2))-1) 	 	
  #bSm = (1-mSm)*((mSm*(1-mSm)/(varSm^2))-1) 
  #Sm = matrix(rbeta(r,aSm,bSm),r,1)	
  Sm = matrix(0,r,1)
  varSm = matrix(0,r,1)
  for (q in 1:length(samples)){
    Sm[q,1] = conecuh$sims.list$mean.phiM[samples[q],h]
    varSm[q,1] = conecuh$sims.list$sigma.M[samples[q],h]
  }  
  #SDmSmi = matrix(rinvgauss(r,varSf^2,1),r,1)
  ASmi = matrix(0,r,1)
  BSmi = matrix(0,r,1)
  Smt = matrix(0,r,t)
  
  ## gJF - transition from juvenile female to adult female
  #mTjf = params[h,5]  			
  #varTjf = paramsSE[h,5]	
  #aTjf = mTjf*((mTjf*(1-mTjf)/varTjf^2)-1)
  #bTjf = (1-mTjf)*((mTjf*(1-mTjf)/varTjf^2)-1)
  #Tjf = matrix(rbeta(r,aTjf,bTjf),r,1)
  Tjf = matrix(0,r,1)
  for (q in 1:length(samples)){
    Tjf[q,1] = tauJF[samples[q],h]
  }    
  #hist(Tjf)
  ATjfi = matrix(0,r,1)
  BTjfi = matrix(0,r,1)
  Tjft = matrix(0,r,t)
  
  ## gJM - transition from juvenile male to adult male
  #mTjm = params[h,6]  			
  #varTjm = paramsSE[h,6]	
  #aTjm = mTjm*((mTjm*(1-mTjm)/varTjm^2)-1)
  #bTjm = (1-mTjm)*((mTjm*(1-mTjm)/varTjm^2)-1)
  #Tjm = matrix(rbeta(r*t,aTjm,bTjm),r,1)
  Tjm = matrix(0,r,1)
  for (q in 1:length(samples)){
    Tjm[q,1] = tauJM[samples[q],h]
  }        
  #hist(Tjm)
  ATjmi = matrix(0,r,1)
  BTjmi = matrix(0,r,1)
  Tjmt = matrix(0,r,t)
  
  ## Pfb - proportion of females breeding
  mPfb = 0.95
  varPfb = 0.04
  aPfb = mPfb*((mPfb*(1-mPfb)/varPfb^2)-1)
  bPfb = (1-mPfb)*((mPfb*(1-mPfb)/varPfb^2)-1)
  Pfb = matrix(rbeta(r*t,aPfb,bPfb),r,1)
  #hist(Pfb)
  APfbi = matrix(0,r,1)
  BPfbi = matrix(0,r,1)
  Pfbt = matrix(0,r,t)
  
  ## Fa - fecundity of adult females (clutch size)
  muFa = params[h,7]
  sdFa = paramsSE[h,7]
  #mean(matrix(round(rlnorm(5000,getParmsLognormForMoments(muFa,sdFa^2)[1],getParmsLognormForMoments(muFa,sdFa^2)[2])),r,1))
  #hist(round(rlnorm(5000,getParmsLognormForMoments(muFa,sdFa^2)[1],getParmsLognormForMoments(muFa,sdFa^2)[2])),r,1)
  Fi = matrix(0,r,1)
  Fit = matrix(0,r,t)
  eggs = matrix(0,r,t)		
  
  ## Pns - probability of nest surviving from predators
  mPns = 0.35
  varPns = 0.04
  aPns = mPns*((mPns*(1-mPns)/varPns^2)-1)
  bPns = (1-mPns)*((mPns*(1-mPns)/varPns^2)-1)
  Pns = matrix(rbeta(r*t,aPns,bPns),r,1)
  #hist(Pns)
  APnsi = matrix(0,r,1)
  BPnsi = matrix(0,r,1)
  Pnst = matrix(0,r,t) 

  ## Ph - probability of eggs hatching (hatching success)
  mPh = 0.85
  varPh = 0.05
  aPh = mPh*((mPh*(1-mPh)/varPh^2)-1)
  bPh = (1-mPh)*((mPh*(1-mPh)/varPh^2)-1)
  Ph = matrix(rbeta(r*t,aPh,bPh),r,1)
  #hist(Ph)
  APhi = matrix(0,r,1)
  BPhi = matrix(0,r,1)
  Pht = matrix(0,r,t)
  
  ## Nj - initial abundance of juveniles
  muNj = params[h,8]
  sdNj = paramsSE[h,8]
  abundShape2Nj = log((sdNj^2)/(muNj^2)+1)
  abundShape1Nj = log(muNj)-1/2*abundShape2Nj
  #round(rlnorm(1,abundShape1Nj,abundShape2Nj))  
  #hist(round(rlnorm(100,abundShape1Nj,abundShape2Nj)))
  Nj = matrix(0,r,t)				
  Nj[,1] = round(rlnorm(r,abundShape1Nj,abundShape2Nj))  # Initial juv abundance
  
  ## Nf - initial abundance of adult females
  muNf = params[h,9]
  sdNf = paramsSE[h,9]
  abundShape2Nf = log((sdNf^2)/(muNf^2)+1)
  abundShape1Nf = log(muNf)-1/2*abundShape2Nf
  #round(rlnorm(1,abundShape1Nf,abundShape2Nf))  
  #hist(round(rlnorm(100,abundShape1Nf,abundShape2Nf)))
  Nf = matrix(0,r,t)				
  Nf[,1] = round(rlnorm(r,abundShape1Nf,abundShape2Nf))  # Initial female abundance
  
  ## Nm - initial abundance of adult males
  muNm = params[h,10]
  sdNm = paramsSE[h,10]
  abundShape2Nm = log((sdNm^2)/(muNm^2)+1)
  abundShape1Nm = log(muNm)-1/2*abundShape2Nm
  #round(rlnorm(1,abundShape1Nm,abundShape2Nm))  
  #hist(round(rlnorm(100,abundShape1Nm,abundShape2Nm)))
  Nm = matrix(0,r,t)				
  Nm[,1] = round(rlnorm(r,abundShape1Nm,abundShape2Nm))  # Initial male abundance
  
## Replication loop; draws replicate-level means for demographics at each stage
  for (i in 1:r){

    AShi[i] = 100*Sh[i]		 	# Hatchling survival	
    BShi[i] = 100*(1-Sh[i])
    
    ASji[i] = 100*Sj[i]			# Juvenile survival 		
    BSji[i] = 100*(1-Sj[i])
    
    ASfi[i] = 100*Sf[i]			# Female survival		
    BSfi[i] = 100*(1-Sf[i])
    
    ASmi[i] = 100*Sm[i]			# Male survival		
    BSmi[i] = 100*(1-Sm[i])

    ATjfi[i] = 100*Tjf[i]			# Transition from juvenile to female 		
    BTjfi[i] = 100*(1-Tjf[i])		
    
    ATjmi[i] = 100*Tjm[i]			# Transition from juvenile to male 		
    BTjmi[i] = 100*(1-Tjm[i])		

    APfbi[i] = 100*Pfb[i]			# Probability of females breeding		
    BPfbi[i] = 100*(1-Pfb[i])
    
    Fi[i] = rlnorm(1,
      getParmsLognormForMoments(muFa,sdFa^2)[1],
      getParmsLognormForMoments(muFa,sdFa^2)[2])

    APnsi[i] = 100*Pns[i]			# Probability of nests survival from predation		
    BPnsi[i] = 100*(1-Pns[i])		
    
    APhi[i] = 100*Ph[i]			# Probability of eggs hatching (hatching success)		
    BPhi[i] = 100*(1-Ph[i])		
    
    ## Projection loop; drawing annual temporal variation in demographic rates
    for(j in 1:t){			

      Sht[i,j] = rbeta(1,AShi[i],BShi[i])			# Hatchling survival
      Sjt[i,j] = rbeta(1,ASji[i],BSji[i])			# Juvenile survival
      Sft[i,j] = rbeta(1,ASfi[i],BSfi[i])			# Female survival
      Smt[i,j] = rbeta(1,ASmi[i],BSmi[i])			# Male survival
      Tjft[i,j] = rbeta(1,ATjfi[i],BTjfi[i])	# Transition juv to female
      Tjmt[i,j] = rbeta(1,ATjmi[i],BTjmi[i])	# Transition juv to male
      Pfbt[i,j] = rbeta(1,APfbi[i],BPfbi[i])  # Prob of females breeding
      Fit[i,j] = rlnorm(1,getParmsLognormForMoments(Fi[i],sdFa^2)[1],
                         getParmsLognormForMoments(Fi[i],sdFa^2)[2])
      Pnst[i,j] = rbeta(1,APnsi[i],BPnsi[i])  # Prob of nests surviving predation
      Pht[i,j] = rbeta(1,APhi[i],BPhi[i])     # Prob. of hatching

      # Projection equation for adult females
      if (j>1) Nf[i,j] = round((Nf[i,j-1]*Sft[i,j-1]+Nj[i,j-1]*Sjt[i,j-1]*Tjft[i,j-1])) 

      # Projection equation for adult males
      if (j>1) Nm[i,j] = round((Nm[i,j-1]*Smt[i,j-1]+Nj[i,j-1]*Sjt[i,j-1]*Tjmt[i,j-1])) 
      
      # Calculate the no. of eggs produced per year, 
      # including females and males; account for prob. of females breeding (Pfb), 
      # whether nest survival from predation (Pns), & imperfect hatching success (Ph)
      eggs[i,j] = round(round(round(round(Nf[i,j]*Pfbt[i,j])*Pnst[i,j])*Fit[i,j])*Pht[i,j])

      # Projection equation for juveniles
      if (j>1) Nj[i,j] = round((Nj[i,j-1]*Sjt[i,j-1]*(1-Tjft[i,j-1]-Tjmt[i,j-1]) + eggs[i,j-1]*Sht[i,j-1])) 
      
      # Calculate population growth rate
      if(j>1) lam[i,j] = (Nf[i,j]+Nm[i,j]+Nj[i,j])/(Nf[i,j-1]+Nm[i,j-1]+Nj[i,j-1])
      
      # Calculate extinction risk, the proportion of replicates that end with 
      # <3 females OR <1 males
      if (Nf[i,j]<3 | Nm[i,j]<1) Pext[i,j]=1 else Pext[i,j]=0 
      
    } # Close projection loop
  } # Close replication loop
  
  N = Nf + Nm + Nj  # Total number 
  
  # Calculate median lambda 
  medlam = apply(lam, 2, median, na.rm=TRUE)
  
  # Calculate extinction risk
  PE = apply(Pext,2,sum)/r
  PEt = PE[t]
  
  # Save all important objects for each site
  assign(paste0("Nf", h), Nf)
  assign(paste0("Nm", h), Nm)
  assign(paste0("Nj", h), Nj)
  assign(paste0("N", h), N)
  assign(paste0("lam", h), lam)
  assign(paste0("medlam", h), medlam)
  assign(paste0("PE", h), PE)
  assign(paste0("PEt", h), PEt)

} # Close site-scenario loop

# Stop the clock
proc.time() - ptm

# Examine abundance estimates; e.g., at Site 1
median(N1)
head(Nf1,20)
head(Nm1,20)
head(N1,20)

# Examine extinction risks at each site
PEt1; PEt2; PEt3; PEt4; PEt5; PEt6
PEt1 = 0.001

### Plot projected abundance at Sites 1-6 
res = rbind(N1,N2,N3,N4,N5,N6)
#res = rbind(Na1,Na2,Na3,Na4,Na5,Na6)
colnames(res)=c(1:100)
res = cbind(c(rep(1,1000),rep(2,1000),rep(3,1000),rep(4,1000),rep(5,1000),rep(6,1000)),res)
max = max(res)

par(mfrow=c(2,3), oma=c(5,5,0,0)+0.5, mar=c(0,0,1,1)+0.5, cex.lab=1.7) 
sites = c("1","2","3","4","5","6")
pe = list(PEt1,PEt2,PEt3,PEt4,PEt5,PEt6)
for (i in 1:s){
  plot(1:t, apply(subset(res, res[,1] == i)[,-c(1)],2,median), type="l", lty=1, lwd=4, ylim=c(0,300), 
       axes=FALSE, cex.lab=1.4, cex.axis=2, xlab="Time (years)", ylab="Abundance")
  # plots median abundance
  if (i > 3){	# restricts axis labels to bottom row
  axis(side=1, lwd=4, cex.axis=1.7)
    } else {axis(side=1, labels=FALSE, lwd=4, cex.axis=1.3)}
  if (i %in%  c(1,4)){	# restricts y-axis labels to panels in c(1,4)
  axis(side=2, lwd=4, cex.axis=1.55)
    } else {axis(side=2, labels=FALSE,lwd=4, cex.axis=1.3)}	
  lines(1:t, apply(subset(res, res[,1] == i)[,-c(1)],2,quantile,probs=c(0.025)), lty=3, lwd=4, col="darkgrey")	# adds lower CI
  lines(1:t, apply(subset(res, res[,1] == i)[,-c(1)],2,quantile,probs=c(0.975)), lty=3, lwd=4, col="darkgrey")	# adds upper CI
  if(i == 2) {legend(2,0.98*max, c("Median abundance","95% CI"), inset=0.05, cex=1,
                               lty=c(1,3), col=c("black","darkgrey"),
                               lwd=c(3,3), box.lwd=2)}
  text(80,0.93*max, paste0("Site ", i), cex=2)
  text(80,0.86*lim[i], parse(text = paste0('P[e] == ', pe[i])), cex=1.5)
  title(xlab = "Time (years)", ylab = "Abundance", outer=TRUE, line=3, cex.sub=3, cex.lab=2.7)
}


### Plot projected abundance at Sites 1-6 
### But re-scale y-axes to have limits that are reasonable
lim = c(40,25,50,70,20,40) 
  # visually examine the results to find reasonable ylims
par(mfrow=c(2,3), oma=c(5,5,0,0)+0.5, mar=c(0,0,1,1)+0.5, cex.lab=1.7) 
sites = c("1","2","3","4","5","6")
pe = list(PEt1,PEt2,PEt3,PEt4,PEt5,PEt6)
for (i in 1:s){
  plot(1:t, apply(subset(res, res[,1] == i)[,-c(1)],2,median), type="l", lty=1, lwd=4,
       ylim=c(0,lim[i]), 
       axes=FALSE, cex.lab=1.4, cex.axis=2, xlab="Time (years)", ylab="Population size (N)")
  # plots median abundance
  if (i > 3){	# restricts axis labels to bottom row
    axis(side=1, lwd=4, cex.axis=1.8)
  } else {axis(side=1, labels=FALSE, lwd=4, cex.axis=2)}
  axis(side=2, lwd=4, cex.axis=1.8)
  lines(1:t, apply(subset(res, res[,1] == i)[,-c(1)],2,quantile,probs=c(0.075)), lty=3, lwd=4, col="darkgrey")	# adds lower CI
  lines(1:t, apply(subset(res, res[,1] == i)[,-c(1)],2,quantile,probs=c(0.925)), lty=3, lwd=4, col="darkgrey")	# adds upper CI
  if(i == 2) {legend(5,0.98*lim[2], c("Median abundance","85% CI"), inset=0.05, cex=1,
                     lty=c(1,3), col=c("black","darkgrey"),
                     lwd=c(3,3), box.lwd=2)}
  text(80,0.95*lim[i], paste0("Site ", i), cex=2)
  #text(20,0.84*lim[i], paste0("Pe = ", pe[i]), cex=1.5)
  text(80,0.86*lim[i], parse(text = paste0('P[e] == ', pe[i])), cex=1.5)
  title(xlab = "Time (years)", ylab=expression(Population~size~italic((N))), outer=TRUE, line=3, cex.sub=3, cex.lab=2.7)
}


# Examine extinction risks and abundances at each site
PEt1; PEt2; PEt3; PEt4; PEt5; PEt6

median(N1[,100])
apply(N1,2,quantile,probs=c(0.025))[100] #lower 95 CI
apply(N1,2,quantile,probs=c(0.975))[100] #upper 95 CI

median(N2[,100])
apply(N2,2,quantile,probs=c(0.025))[100] #lower 95 CI
apply(N2,2,quantile,probs=c(0.975))[100] #upper 95 CI

median(N3[,100])
apply(N3,2,quantile,probs=c(0.025))[100] #lower 95 CI
apply(N3,2,quantile,probs=c(0.975))[100] #upper 95 CI

median(N4[,100])
apply(N4,2,quantile,probs=c(0.025))[100] #lower 95 CI
apply(N4,2,quantile,probs=c(0.975))[100] #upper 95 CI

median(N5[,100])
apply(N5,2,quantile,probs=c(0.025))[100] #lower 95 CI
apply(N5,2,quantile,probs=c(0.975))[100] #upper 95 CI

median(N6[,100])
apply(N6,2,quantile,probs=c(0.025))[100] #lower 95 CI
apply(N6,2,quantile,probs=c(0.975))[100] #upper 95 CI



########## Section 4)
########## Use R packages to perform sensitivity analysis
library(popbio)

# Use the 'params' and 'paramsSE' objects from previous Section  
params
paramsSE

# Run a series of for-loops to analyze demographic sensitivities for each population 1-6

#### Use for-loops to simulate and projection population demographics across sites

r = 1000 # replicates

## Start a timer!
ptm <- proc.time()

## Population-level loop
for (h in c(1:(s))){    # Subset to research site, i
  
  ## Sh - survival of hatchlings
  mSh = params[h,1]                 # Mean hatchling survival		
  varSh = paramsSE[h,1]							# Variance of mean
  aSh = mSh*((mSh*(1-mSh)/(varSh^2))-1) 	 	
  bSh = (1-mSh)*((mSh*(1-mSh)/(varSh^2))-1) 
  Shi = matrix(rbeta(r,aSh,bSh),r,1)	# Parametric uncertainty 
  #hist(Shi)
  #SDmShi = matrix(rinvgauss(r,varSh^2,1),r,1)
  AShi = matrix(0,r,1)				# beta distribution shape parameters
  BShi = matrix(0,r,1)				# beta distribution shape parameters
  Sht = matrix(0,r,t)					# annual variation survival
  
  ## Sj - survival of juveniles 
  # sampled from the posterior distribution of our estimation analysis above
  #sample(conecuh$sims.list$mean.phiJ[,h],1)
  
  ## Sf - survival of females 
  #sample(conecuh$sims.list$mean.phiF[,h],1)
  
  ## Sm - survival of males 
  #sample(conecuh$sims.list$mean.phiM[,h],1)
  
  ## gJF - transition from juvenile to female
  #sample(conecuh$sims.list$mean.tauJF[,h],1)
  
  ## gJM - transition from juvenile to male
  #sample(conecuh$sims.list$mean.tauJM[,h],1)
  
  ## Pfb - proportion of females breeding
  mPfb = 0.95
  varPfb = 0.04
  aPfb = mPfb*((mPfb*(1-mPfb)/varPfb^2)-1)
  bPfb = (1-mPfb)*((mPfb*(1-mPfb)/varPfb^2)-1)
  Pfb = matrix(rbeta(r*t,aPfb,bPfb),r,1)
  #hist(Pfb)
  APfbi = matrix(0,r,1)
  BPfbi = matrix(0,r,1)
  Pfbt = matrix(0,r,t)
  
  ## Pns - probability of nest surviving from predators
  mPns = 0.35
  varPns = 0.04
  aPns = mPns*((mPns*(1-mPns)/varPns^2)-1)
  bPns = (1-mPns)*((mPns*(1-mPns)/varPns^2)-1)
  Pns = matrix(rbeta(r*t,aPns,bPns),r,1)
  #hist(Pns)
  APnsi = matrix(0,r,1)
  BPnsi = matrix(0,r,1)
  Pnst = matrix(0,r,t) 

  ## Fa - fecundity of adult females (clutch size)
  muFa = params[h,7]
  sdFa = paramsSE[h,7]
  #mean(matrix(round(rlnorm(5000,getParmsLognormForMoments(muFa,sdFa^2)[1],getParmsLognormForMoments(muFa,sdFa^2)[2])),r,1))
  #hist(round(rlnorm(5000,getParmsLognormForMoments(muFa,sdFa^2)[1],getParmsLognormForMoments(muFa,sdFa^2)[2])),r,1)
  Fi = matrix(0,r,1)
  Fit = matrix(0,r,t)
  eggs = matrix(0,r,t)		
  
  ## Ph - probability of eggs hatching (hatching success)
  mPh = 0.85
  varPh = 0.05
  aPh = mPh*((mPh*(1-mPh)/varPh^2)-1)
  bPh = (1-mPh)*((mPh*(1-mPh)/varPh^2)-1)
  Ph = matrix(rbeta(r*t,aPh,bPh),r,1)
  #hist(Ph)
  APhi = matrix(0,r,1)
  BPhi = matrix(0,r,1)
  Pht = matrix(0,r,t)
  
  ## Nj - initial abundance of juvenile females
  muNj = params[h,8]
  sdNj = paramsSE[h,8]
  abundShape2Nj = log((sdNj^2)/(muNj^2)+1)
  abundShape1Nj = log(muNj)-1/2*abundShape2Nj
  #round(rlnorm(1,abundShape1Nj,abundShape2Nj))  
  #hist(round(rlnorm(100,abundShape1Nj,abundShape2Nj)))
  Nj = matrix(0,r,t)				
  Nj[,1] = round(rlnorm(r,abundShape1Nj,abundShape2Nj))  # Initial juv abundance
  
  ## Nf - initial abundance of adult females
  muNf = params[h,9]
  sdNf = paramsSE[h,9]
  abundShape2Nf = log((sdNf^2)/(muNf^2)+1)
  abundShape1Nf = log(muNf)-1/2*abundShape2Nf
  #round(rlnorm(1,abundShape1Nf,abundShape2Nf))  
  #hist(round(rlnorm(100,abundShape1Nf,abundShape2Nf)))
  Nf = matrix(0,r,t)				
  Nf[,1] = round(rlnorm(r,abundShape1Nf,abundShape2Nf))  # Initial female abundance
  
  ## Nm - initial abundance of adult males
  muNm = params[h,10]
  sdNm = paramsSE[h,10]
  abundShape2Nm = log((sdNm^2)/(muNm^2)+1)
  abundShape1Nm = log(muNm)-1/2*abundShape2Nm
  #round(rlnorm(1,abundShape1Nm,abundShape2Nm))  
  #hist(round(rlnorm(100,abundShape1Nm,abundShape2Nm)))
  Nm = matrix(0,r,t)				
  Nm[,1] = round(rlnorm(r,abundShape1Nm,abundShape2Nm))  # Initial male abundance
  
  ## Matrices to save results
  lam = matrix(0,r)  # Population growth; lambda
  SSD = matrix(0,r,3) # Stable-stage distribution
  colnames(SSD)=c("J","F","M")
  GT = matrix(0,r)    # Generation time
  RV = matrix(0,r,3)  # Reproductive values
  elas = matrix(0,r,9)  # Elasticity values
  colnames(elas) = c("JJ","JF","JM",
                      "FJ","FF","FM",
                      "MJ","MF","MM")
  
  #### Build a matrix model for each population, h
  #### Run an elasticity analysis, replicated r time, and save the mean results
  for (i in 1:r){
  stages = c("J","F","M")

  lam[i,] = eigen.analysis(matrix(c(
      sample(conecuh$sims.list$mean.phiJ[,h],1)*(1-sample(tauJF[,h],1)-sample(tauJM[,h],1)),rbeta(1,aPfb,bPfb)*rlnorm(1,getParmsLognormForMoments(muFa,sdFa^2)[1],getParmsLognormForMoments(muFa,sdFa^2)[2])*rbeta(1,aPns,bPns)*rbeta(1,aPh,bPh)*rbeta(1,aSh,bSh),0,
      sample(conecuh$sims.list$mean.phiJ[,h],1)*sample(tauJF[,h],1),sample(conecuh$sims.list$mean.phiF[,h],1),0,
      sample(conecuh$sims.list$mean.phiJ[,h],1)*sample(tauJM[,h],1),0,sample(conecuh$sims.list$mean.phiM[,h],1)),
      nrow=3,byrow=TRUE,dimnames=list(stages,stages)))$lambda1
  
  SSD[i,] = eigen.analysis(matrix(c(
    sample(conecuh$sims.list$mean.phiJ[,h],1)*(1-sample(tauJF[,h],1)-sample(tauJM[,h],1)),rbeta(1,aPfb,bPfb)*rlnorm(1,getParmsLognormForMoments(muFa,sdFa^2)[1],getParmsLognormForMoments(muFa,sdFa^2)[2])*rbeta(1,aPns,bPns)*rbeta(1,aPh,bPh)*rbeta(1,aSh,bSh),0,
    sample(conecuh$sims.list$mean.phiJ[,h],1)*sample(tauJF[,h],1),sample(conecuh$sims.list$mean.phiF[,h],1),0,
    sample(conecuh$sims.list$mean.phiJ[,h],1)*sample(tauJM[,h],1),0,sample(conecuh$sims.list$mean.phiM[,h],1)),
    nrow=3,byrow=TRUE,dimnames=list(stages,stages)))$stable.stage

  GT[i,] = generation.time(matrix(c(
    sample(conecuh$sims.list$mean.phiJ[,h],1)*(1-sample(tauJF[,h],1)-sample(tauJM[,h],1)),rbeta(1,aPfb,bPfb)*rlnorm(1,getParmsLognormForMoments(muFa,sdFa^2)[1],getParmsLognormForMoments(muFa,sdFa^2)[2])*rbeta(1,aPns,bPns)*rbeta(1,aPh,bPh)*rbeta(1,aSh,bSh),0,
    sample(conecuh$sims.list$mean.phiJ[,h],1)*sample(tauJF[,h],1),sample(conecuh$sims.list$mean.phiF[,h],1),0,
    sample(conecuh$sims.list$mean.phiJ[,h],1)*sample(tauJM[,h],1),0,sample(conecuh$sims.list$mean.phiM[,h],1)),
    nrow=3,byrow=TRUE,dimnames=list(stages,stages)))
  
  RV[i,] = eigen.analysis(matrix(c(
    sample(conecuh$sims.list$mean.phiJ[,h],1)*(1-sample(tauJF[,h],1)-sample(tauJM[,h],1)),rbeta(1,aPfb,bPfb)*rlnorm(1,getParmsLognormForMoments(muFa,sdFa^2)[1],getParmsLognormForMoments(muFa,sdFa^2)[2])*rbeta(1,aPns,bPns)*rbeta(1,aPh,bPh)*rbeta(1,aSh,bSh),0,
    sample(conecuh$sims.list$mean.phiJ[,h],1)*sample(tauJF[,h],1),sample(conecuh$sims.list$mean.phiF[,h],1),0,
    sample(conecuh$sims.list$mean.phiJ[,h],1)*sample(tauJM[,h],1),0,sample(conecuh$sims.list$mean.phiM[,h],1)),
    nrow=3,byrow=TRUE,dimnames=list(stages,stages)))$repro.value
  
  elas[i,] = eigen.analysis(matrix(c(
    sample(conecuh$sims.list$mean.phiJ[,h],1)*(1-sample(tauJF[,h],1)-sample(tauJM[,h],1)),rbeta(1,aPfb,bPfb)*rlnorm(1,getParmsLognormForMoments(muFa,sdFa^2)[1],getParmsLognormForMoments(muFa,sdFa^2)[2])*rbeta(1,aPns,bPns)*rbeta(1,aPh,bPh)*rbeta(1,aSh,bSh),0,
    sample(conecuh$sims.list$mean.phiJ[,h],1)*sample(tauJF[,h],1),sample(conecuh$sims.list$mean.phiF[,h],1),0,
    sample(conecuh$sims.list$mean.phiJ[,h],1)*sample(tauJM[,h],1),0,sample(conecuh$sims.list$mean.phiM[,h],1)),
    nrow=3,byrow=TRUE,dimnames=list(stages,stages)))$elasticities

  }
  # Save all important objects for each site
  site = list(lam,SSD,GT,RV,elas)
  names(site) = c("lam","SSD","GT","RV","elas")
  assign(paste0("site", h), site)
  
} # Close site loop

# Stop the clock
proc.time() - ptm

# E.g., 
apply(site1$lam,2,median)
apply(site1$lam,2,mean)
apply(site1$SSD,2,median)
apply(site1$RV,2,median)
apply(site1$elas,2,median)

apply(site2$lam,2,median)
apply(site2$SSD,2,median)
apply(site2$RV,2,median)
apply(site2$elas,2,median)

apply(site3$lam,2,median)
apply(site3$SSD,2,median)
apply(site3$RV,2,median)
apply(site3$elas,2,median)

apply(site4$lam,2,median)
apply(site4$SSD,2,median)
apply(site4$RV,2,median)
apply(site4$elas,2,median)

apply(site5$lam,2,median)
apply(site5$SSD,2,median)
apply(site5$RV,2,median)
apply(site5$elas,2,median)

apply(site6$lam,2,median)
apply(site6$SSD,2,median)
apply(site6$RV,2,median)
apply(site6$elas,2,median)



#### Graph summary results for demographic parameters at each site
#library(ggplot2)
library(ggpubr)

## Lamda -- Plot and compare lambda among sites
x = 1000
sites = c(rep(1,x),rep(2,x),rep(3,x),rep(4,x),rep(5,x),rep(6,x))
df = rbind(site1$lam,site2$lam,site3$lam,site4$lam,site5$lam,site6$lam)
lam = data.frame(cbind(sites,df))
colnames(lam) = c("Site","Lambda")

ggviolin(lam, x = "Site", y = "Lambda", fill = "Site", 
         palette = c("lightblue","red","lightblue","blue","red","pink"),
         ylim=c(0.85,1.15), lty=2, add="median", font.y=16, font.x=16,
         ylab=expression(paste("Population growth (", lambda, ")")))

# Save median (95% CI) for each site in a df
mLAM = matrix(NA,6,1) # MEDIAN
for (i in 1:6){
  site = subset(lam, lam$Site == i)
  mLAM[i,1] = signif(median(site[,2]),3)
}

mLAMci = matrix(NA,6,2) # 95% CI
for (i in 1:6){
  site = subset(lam, lam$Site == i)
  mLAMci[i,1] = signif(apply(site,2,quantile,probs=c(0.025))[2],3)
  mLAMci[i,2] = signif(apply(site,2,quantile,probs=c(0.975))[2],3)
}
colnames(mLAMci) = c("Low 95%", "High 95%")

## Stable-stage distribution (SSD) -- Plot and compare SSD among sites
library(reshape2) # note - this graph doesn't work well
x=1000
Site = c(rep(1,x),rep(2,x),rep(3,x),rep(4,x),rep(5,x),rep(6,x))
df = rbind(site1$SSD,site2$SSD,site3$SSD,site4$SSD,site5$SSD,site6$SSD) 
ssd = data.frame(cbind(Site,df))
ssd = melt(ssd, id.vars=c("Site"))
colnames(ssd) = c("Site","Stage","Value")
    
ggviolin(ssd, x = "Site", y = "Value", fill = "Stage", font.y=16, font.x=16,
         palette = c("lightgreen","red","blue"),
         lty=2, ylab="Stable-stage distribution", add="mean_sd")

# save results in a df
mSSD = matrix(NA,6,3) #MEDIAN
for (i in 1:6){
  site = subset(ssd, ssd$Site == i)
  j = subset(site, site$Stage == "J")
  fem = subset(site, site$Stage == "F")
  male = subset(site, site$Stage == "M")
  mSSD[i,] = signif(c(median(j[,3]),median(fem[,3]),median(male[,3])),3)
}

mSSDciLOW = matrix(NA,6,3) # LOW CI
for (i in 1:6){
  site = subset(ssd, lam$Site == i)
  j = subset(site, site$Stage == "J")
  fem = subset(site, site$Stage == "F")
  male = subset(site, site$Stage == "M") 
  mSSDciLOW[i,1:3] = c(apply(matrix(j[,3]),2,quantile,probs=c(0.025)),
                     apply(matrix(fem[,3]),2,quantile,probs=c(0.025)),
                     apply(matrix(male[,3]),2,quantile,probs=c(0.025)))
}
colnames(mSSDciLOW) = c("J", "F", "M")

mSSDciHIGH = matrix(NA,6,3) # LOW CI
for (i in 1:6){
  site = subset(ssd, lam$Site == i)
  j = subset(site, site$Stage == "J")
  fem = subset(site, site$Stage == "F")
  male = subset(site, site$Stage == "M") 
  mSSDciHIGH[i,1:3] = c(apply(matrix(j[,3]),2,quantile,probs=c(0.975)),
                       apply(matrix(fem[,3]),2,quantile,probs=c(0.975)),
                       apply(matrix(male[,3]),2,quantile,probs=c(0.975)))
}
colnames(mSSDciHIGH) = c("J", "F", "M")


## Reproductive value (RV) of females -- Plot and compare RV among sites
x=1000
Site = c(rep(1,x),rep(2,x),rep(3,x),rep(4,x),rep(5,x),rep(6,x))
df = rbind(site1$RV,site2$RV,site3$RV,site4$RV,site5$RV,site6$RV) 
rv = data.frame(cbind(Site,df))
colnames(rv) = c("Site","J","F","M")
rv = melt(rv, id.vars=c("Site"))
colnames(rv) = c("Site","Stage","Value")

# save results in a df
mRV = matrix(NA,6,1)
for (i in 1:6){
  site = subset(rv, rv$Site == i & rv$Stage == "F")
  mRV[i,] = signif(median(site[,3]),3)
}
colnames(mRV) = "ReproValue"

mRVci = matrix(NA,6,2) # 95% CI
for (i in 1:6){
  site = subset(rv, rv$Site == i & rv$Stage == "F")
  site = site[,-c(2)]
  mRVci[i,1] = signif(apply(site,2,quantile,probs=c(0.025))[2],3)
  mRVci[i,2] = signif(apply(site,2,quantile,probs=c(0.975))[2],3)
}
colnames(mRVci) = c("Low 95%", "High 95%")

## Generation time (G) -- Plot and compare G  among sites
x = 1000
sites = c(rep(1,x),rep(2,x),rep(3,x),rep(4,x),rep(5,x),rep(6,x))
df = rbind(site1$GT,site2$GT,site3$GT,site4$GT,site5$GT,site6$GT)
gt = data.frame(cbind(sites,df))
colnames(gt) = c("Site","GT")

ggviolin(gt, x = "Site", y = "GT", fill = "Site", font.y=16, font.x=16,
         palette = c("blue","red","blue","lightblue","red","pink"),
         ylim=c(0,600), lty=2, ylab="Generation time (yr)", add="median")

# save results in a df
mG = matrix(NA,6,1)
for (i in 1:6){
  site = subset(gt, gt$Site == i)
  mG[i,] = signif(c(median(site[,2])),3)
}

mGci = matrix(NA,6,2) # 95% CI
for (i in 1:6){
  site = subset(gt, gt$Site == i)
  mGci[i,1] = signif(apply(site,2,quantile,probs=c(0.025))[2],3)
  mGci[i,2] = signif(apply(site,2,quantile,probs=c(0.975))[2],3)
}
colnames(mGci) = c("Low 95%", "High 95%")


## Elasticity -- Plot and compare elasticity values
## at 1) a single site, and 2) among all sites

# 1) single site
site = site4 
x=1000
Parameter = c(rep("S_J",x),rep("T_JF",x),rep("Recruitment",x),rep("S_F",x))
df = c(site$elas[,c(1)],site$elas[,c(2)],site$elas[,c(4)],site$elas[,c(5)])
elas = data.frame(Parameter,df)
colnames(elas) = c("Parameter","Value")
str(elas)

ggviolin(elas, x = "Parameter", y = "Value", fill = "Parameter", 
         palette = c("pink","lightblue","orange","red"), 
         font.y = 18, xlab=FALSE, order=c("Recruitment","S_JF","T_JF","S_F"),
         ylim=c(0,1.2), lty=2, ylab="Elasticity", add="median")

# 2) All sites; using multi-facet plot
sites = list(site1,site2,site3,site4,site5,site6)
for (i in 1:6){
  site = sites[[i]]
  x=1000
  Parameter = c(rep("S_J",x),rep("T_JF",x),rep("Recruitment",x),rep("S_F",x))
  df = c(site$elas[,c(1)],site$elas[,c(2)],site$elas[,c(4)],site$elas[,c(5)])
  Site = rep(i,x*4)
  elas = data.frame(Site,Parameter,df)
  colnames(elas) = c("Site","Parameter","Value")
  assign(paste0("elas", i), elas) 
}
elas = rbind(elas1,elas2,elas3,elas4,elas5,elas6)
head(elas)

# Faceted by site
(p0 = ggviolin(elas, x = "Parameter", y = "Value", fill = "Parameter", 
         palette = c("pink","lightblue","orange","red"), 
         facet.by = "Site", legend=NULL,
         xlab="Parameter", font.y = 18, font.x=18, font.tickslab=16, 
         order=c("Recruitment","S_J","T_JF","S_F"),
         ylim=c(0,1), lty=2, ylab="Elasticity", add="median"))

(p1 = p0 + scale_x_discrete(labels = 
        c('Recruitment' = expression(italic('R')),
        'S_J' = bquote(phi^j),
        'T_JF' = bquote(tau^jf), 
        'S_F' = bquote(phi^f))))


# Faceted by parameter - not my favorite
(p0 = ggviolin(elas, x = "Site", y = "Value", fill = "Site", 
        #palette = c("pink","lightblue","orange","red"), 
        facet.by = "Parameter", legend=NULL,
        xlab="Site", font.y = 18, font.x=18, font.tickslab=16, 
        #order=c("Recruitment","PhiJ","TauJF","PhiF"),
        panel.labs=list(Parameter = c(expression(italic('R')),
            bquote(phi[j]),bquote(tau[jf]),bquote(phi[f]))),
        ylim=c(0,1), lty=2, ylab="Elasticity", add="median"))



##### TABLE summarizing matrix results
##### i.e., lamda, SSD, reproductive values, and generation times

# Median estimates
res = data.frame(cbind(c(1:6),mLAM,mSSD,mRV,mG))
colnames(res) = c("Site","Lamda","SSD_J","SSD_F","SSD_M",
                    "Repro_F","GenTime")
res

#install.packages("clipr")
write_clip(res)

# 95% CI
ci = data.frame(cbind(c(1:6),mLAMci,mSSDciLOW[,1],mSSDciHIGH[,1],
        mSSDciLOW[,2],mSSDciHIGH[,2],mSSDciLOW[,3],mSSDciHIGH[,3],
        mRVci,mGci))
colnames(ci) = c("Site","LamLow","LamHigh","SSD_J_Low","SSD_J_High",
                 "SSD_F_Low","SSD_F_High","SSD_M_Low","SSD_M_High",
                 "ReproV_Low","ReproV_High","GenTime_Low","GenTime_High")
ci
write_clip(ci)







