#model simplu #

##################################################################################
############################ CODING BLOCK ########################################
##################################################################################
# Run-Time Environment:               R version 3.5.1
# Author:			      M.J. van Esdonk
# Adjustments:			S. Cristea
# Short title:			      Pharmacokinetic simulations in R for amikacin
# Version:			      V.1.0
##################################################################################
##################################################################################
###
rm(list = ls(all.names = TRUE))
###

version <- 1.0

###############################################
############ Load libraries
library(mrgsolve)
library(dplyr)
library(ggplot2)

##############################################
## Covariate estimates 
bw  <- 3520 # Birth weight (g)
pna  <- 5 # Postnatal age (days)
cw  <- 3400 # Current weight (g)
nsaid <- 0 # Ibuprofen administration (1 - yes, 0 - no)


###############################################
### Dosing

inf <- T
dose <- 15 * (cw/1000) # mg/kg
ii1 <- 24 # dosing interval in hours

## Infusion duration

d2 <- 20 # min 
dur <- d2 / 60 # duration in hours 


###############################################
## Parameter estimates
# these are provided in the model file




###############################################
## Insert etas and sigmas

etacl  <- 0.0899 # there is only IIV on CL
etavd  <- 0
etavd2 <- 0
etaq1  <- 0

#################################################
## Insert residual error

sigmaprop <- 0.0614 # Proportional error
sigmaadd  <- 0.267 # Additive error

###############################3
## Simulation info

nsamples <- 1000 ### Number of simulated individuals
sim_time <- 120 ## Time of simulation

# set probabilities of ribbon in figure
minprobs=0.1
maxprobs=0.9

###################################################
############# Set Dosing objects

if(inf){
  ## IV infusion 
  Administration <-  as.data.frame(ev(ID=1:nsamples,ii=ii1, cmt=1, addl=9999, amt=dose, rate = dose/dur,time=0)) 
}else{
  ## IV BOLUS
  Administration <-  as.data.frame(ev(ID=1:nsamples,ii=ii1, cmt=1, addl=9999, amt=dose, rate = 0,time=0)) 
}

## Sort by ID and time
data <- Administration[order(Administration$ID,Administration$time),]


######################################################################
### Load in model for mrgsolve
mod <- mread_cache("amik_popPK")

## Specify the omegas and sigmas matrices
omega <- cmat(etacl,
              0,etavd,
              0,0,etavd2,
              0,0,0,etaq1)


sigma <- cmat(sigmaprop,
              0,sigmaadd)

## Set parameters in dataset
# these will taken from the model file

## Set covariates in dataset
data$BW <- bw
data$PNA <- pna
data$CW <- cw
data$NSAID <- nsaid

#################################################################
###### Perform simulation
out <- mod %>%
  data_set(data) %>%
  omat(omega) %>%
  smat(sigma) %>%
  mrgsim(end=sim_time,delta=sim_time/100, obsonly=TRUE) # Simulate 100 observations, independent of the total simulated time


### Output of simulation to dataframe
df <- as.data.frame(out)

#################################################################  
## Calculate summary statistics

sum_stat <- df %>%
  group_by(time) %>%
  summarise(Median_C=median(DV),
            Low_percentile=quantile(DV,probs=minprobs),
            High_percentile=quantile(DV,probs=maxprobs)
  )  


################################### Graphical output

plot_pk <- ggplot(sum_stat, aes(x=time,y=Median_C)) +
  
  ## Add ribbon for variability
  geom_ribbon(aes(ymin=Low_percentile, ymax=High_percentile, x=time), alpha = 0.15, linetype=0)+
  
  ## Add median line
  geom_line(size=2) +
  
  # scale_y_log10()+
  ## Add therapeutic target:
  geom_abline(slope = 0, intercept = 35, col = "red", linetype = "dashed")+
  geom_abline(slope = 0, intercept = 24, col = "green", linetype = "dashed")+
  geom_abline(slope = 0, intercept = 3, col = "red", linetype = "dashed")+
  geom_abline(slope = 0, intercept = 1.5, col = "green", linetype = "dashed")+
  
  # Set axis and theme
  ylab(paste("Concentration",sep=""))+
  xlab("Time after dose (h)")+
  theme_bw()+
  
  # Remove legend
  theme(legend.position="none")

