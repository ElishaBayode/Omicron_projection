# Base Least squares code taken from BC_projection.R script - differs at line 83

########## Libraries and source files
require(deSolve)
require(CanCovidData)
library(plyr)
library(tidyverse)
library(tidyr)
require(ggplot2)
require(reshape2)
library(lubridate)
library(dplyr)
library(data.table)

set.seed(3242)

source("analysis-new/functions.R")
source("analysis-new/likelihood_func.R")
source("analysis-new/BC_data.R")

forecasts_days <- 30 
intro_date <-  ymd("2021-11-20")
stop_date <- ymd("2022-01-24") # make it uniform 
#import data 
#run BC_data.R (preferably line by line to check if there are 0 cases and NA's)
dat = readRDS("data/BC-dat.rds")
#include Omicron wave only
dat <- dat %>% filter(date >= intro_date &  date <= stop_date)
dat_omic <- dat
dat_omic <- filter(dat_omic, date >= intro_date) %>% select(c("day", "value"))
dat_omic$day <- 1:nrow(dat_omic)
test_prop_BC <- filter(mytest_BC, date >= intro_date)$test_prop


#fit (by eyeballing) test_prop to a sigmoid function 
#test_prop_BC1 <- c(test_prop_BC[1:length(dat_omic$value)], rep(last(test_prop_BC),forecasts_days)) #ensuring the length is consistent
fake_test_prop_BC <- (1 - (1-0.2)/(1 + exp(-0.25*(1:(nrow(dat_omic)+forecasts_days)-45))))
plot(fake_test_prop_BC)
lines(test_prop_BC)
abline(v=length(dat_omic$day)) # where it will be cut off at the next line:

#subset for fitting alone (length of data)
test_prop_BC <- fake_test_prop_BC[1:length(dat_omic$day)]
test_prop <- test_prop_BC 


#set values to generate initial conditions with make_init()
N=5.07e6
N_pop=N
vaxlevel_in = 0.82 # portion of the pop vaccinated at start time 
port_wane_in = 0.04 # portion boosted at start tie 
past_infection_in = 0.1  #increased this from 0.1 to 0.18 # total in R at start time
incres_in = 470 # resident strain (delta) incidence at start 
incmut_in = 3 # new (omicron) inc at stat 
simu_size = 1e5 # number of times to resample the negative binom (for ribbons)
forecasts_days =30 # how long to forecast for 
times = 1:nrow(dat_omic)




#declaring fixed parameters 
eff_date <-   ymd("2021-12-29")  # intervention date 
parameters <-         c(sigma=1/3, # incubation period (3 days) (to fixed)
                        gamma=1/(5), #recovery rate (fixed)
                        nu =0.007, #vax rate: 0.7% per day (fixed)
                        mu=1/(82*365), # 1/life expectancy (fixed)
                        w1= 1/(0.5*365),# waning rate from R to S (fixed)
                        w2= 1/(0.5*365), # waning rate from Rv to V (fixed)
                        w3= 1/(0.5*365),# waning rate Rw to W (fixed)
                        ve=1, # I think this should be 1. it is not really efficacy  ( fixed)
                        beta_r=0.555, #transmission rate (to estimate) (0.35)
                        beta_m=1.1, #transmission rate (to estimate)(*1.9)
                        epsilon_r = (1-0.8), # % this should be 1-ve 
                        epsilon_m = 1-0.3, #(1-0.25)?(1-0.6), # % escape capacity #(fixed)
                        b= 0.006, # booster rate  (fixed)
                        beff = 0.7, # booster efficacy
                        wf=0.05, # protection for newly recovered #0.2
                        N=5e6,
                        stngcy= 0.4,#0.78, #(*%(reduction)) strength of intervention (reduction in beta's)
                        eff_t = as.numeric(eff_date - intro_date),
                        p = 0.5, #negative binomial mean
                        theta = 0.1 #negative binomial dispersion
                        
)
init <- make_init()   #generate initial states
parameters_BC <- parameters


source("analysis-new/mod_fitting_setup.R")




# ---- Least squares fitting



### Super simple grid search:
beta_m_vals <- seq(0, 5, length.out=100)
result<-0
for (i in 1:length(beta_m_vals)){
  parameters["beta_m"] <- beta_m_vals[i]
  # Get model simulation with that beta_n
  model_from_beta <- as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                                       parms=parameters)) 
  # Extract the incident cases per day
  getinc = get_total_incidence(model_from_beta, parameters = parameters)
  getinc$date = dat$date[1]+getinc$time
  tprop_vec = c(test_prop, rep(test_prop[length(test_prop)], nrow(model_from_beta)-length(test_prop))) 
  getinc$inc_reported = getinc$inc_tot*tprop_vec
  result[i] <- (sum((getinc$inc_reported - dat_omic$value)^2))
}

beta_m_vals[which.min(result)]
#this catches estimated parameter values from MLE , and adds them to 'parameters' structure
parameters["beta_m"] <- beta_m_vals[which.min(result)]




# OR
### Slightly-better optim!
# fitting any subset of parameters
guess <- c(beta_m = 1.8, p = 0.1)

ls_fit <- function(par,test_prop,dat_omic,parameters){
  
  if(!all(names(par) %in% names(parameters))) stop('Names of pars to fit do not match names in parameters object')
  parameters[names(par)] <- par
  times <- 1:nrow(dat_omic)


  # Get model simulation with that set of params
  model_from_beta <- as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                                                parms=parameters)) 
  # Extract the incident cases per day
  getinc = get_total_incidence(model_from_beta, parameters = parameters)
  tprop_vec = c(test_prop, rep(test_prop[length(test_prop)], nrow(model_from_beta)-length(test_prop))) 
  getinc$inc_reported = getinc$inc_tot*tprop_vec
  result[i] <- (sum((getinc$inc_reported - dat_omic$value)^2))
}

fit_BC <- optim(fn=ls_fit,  par=guess, lower=c(0), 
                upper = c(Inf), method = "L-BFGS-B", parameters = parameters_BC,
                test_prop=test_prop, dat_omic=dat_omic)

fit_BC
#this catches estimated parameter values from MLE , and adds them to 'parameters' structure
parameters[names(guess)] <- fit_BC$par








#make prediction and projection with estimated parameters 
times <- 1:(nrow(dat_omic) + forecasts_days)
out_BC <- as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                                     parms=parameters)) 



# Decision here about whether to project test_prop forward (use fake_test_prop_BC) or keep same test_prop
# as last day (use thisvec)

#with test_prop 
incidence_BC =  parameters[["sigma"]]*(out_BC$Er + out_BC$Erv + out_BC$Erw +
                                         out_BC$Em + out_BC$Emv +
                                         out_BC$Emw)*fake_test_prop_BC



uncert_bound_BC = raply(simu_size,rnbinom(n=length(incidence_BC),
                                          mu=parameters[["p"]]*incidence_BC,
                                          size=1/parameters[["theta"]]))

project_dat_BC =  as.data.frame(aaply(uncert_bound_BC 
                                      ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
  mutate(date=seq.Date(ymd(intro_date),
                       ymd(intro_date)-1+length(times), 1))

#####################check fit ######################

#add dat to data for plotting 
dat_reported <- dat_omic  %>% mutate(date=seq.Date(ymd(intro_date),
                                                   ymd(intro_date)-1+length(dat_omic$day), 1))

ggplot() + geom_line(data=project_dat_BC,aes(x=date,y=`50%`), col="green",size=1.5,alpha=0.4) +
  geom_ribbon(data=project_dat_BC,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='darkgreen',alpha=0.1, size = 1.5)+
  geom_point(data=dat_reported,aes(x=date, y=value),color='grey48', alpha=0.8, size = 1.5)



