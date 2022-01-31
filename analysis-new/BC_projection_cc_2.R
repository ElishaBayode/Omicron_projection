require(deSolve)
require(ggplot2)
require(reshape2)
library(lubridate)
library(dplyr)
library(data.table)


source("analysis-new/mod_fitting_setup.R")
source("analysis-new/likelihood_func.R")


# cc to do list
# setup stuff 
# - set up default parameters better (efficacy, beff etc) 
# - if we are not fitting beta_r, ensure that the initial delta decline is reasonable
# - alternatively we may want to start the model earlier
# so that we have more information about beta_r, because omicron has not happened yet 

# fit approach stuff 
# - figure out what is identifiable: what can we really fit? 
# - decide what parameters are being fitted
# - show fit compared to data - even just the fit without resampling
# - set up a section that does some  further sanity checking for the fit: 
# --- is it just getting the bounds we put in? 
# --- does the parameter seem reasonable based on what we know about omicron vs delta? 
# --- are the trajectories sensible? what about portion infections in vax vs unvax? 

# - check against positivity, wastewater to ensure that where we see a decline it 
# is mirrored in these other data (waiting to hear on wastewater. have source for positivity) 
# - make longer-term plots exploring the duration of immunity and its impact
# - when doing that, explore the impact of vaccination rate


forecasts_days <- 30 
intro_date <-  ymd("2021-11-20")
stop_date <- ymd("2022-01-24")  
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
fake_test_prop_BC <- (1 - (1-0.05)/(1 + exp(-0.25*(1:(length(test_prop_BC)+forecasts_days)-35))))
plot(fake_test_prop_BC)
lines(test_prop_BC)

#subset for fitting alone (length of data)
test_prop_BC <- fake_test_prop_BC[1:length(dat_omic$day)]
test_prop <- test_prop_BC 


# ---- initialization ---- 
# set values to generate initial conditions with make_init()
# http://www.bccdc.ca/Health-Info-Site/Documents/VoC/VoC_Weekly_20220121.pdf
# note: maybe 10% at dec 1 and about 50% by dec 12 
N=5.07e6
N_pop=N
#ascFrac <- 0.5
vaxlevel_in = 0.82
port_wane_in = 0.04 # this is the portion *boosted* at the start time 
past_infection_in = 0.1  #increased this from 0.1 to 0.18 
incres_in = 470
incmut_in = 3
simu_size = 1e5
forecasts_days =30
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
                        #beta_m=0.8*2.2, #transmission rate (to estimate)(*1.9)
                        epsilon_r = (1-0.8), # % this should be 1-ve 
                        epsilon_m = 1-0.3, #(1-0.25)?(1-0.6), # % escape capacity #(fixed)
                        b= 0.006, # booster rate  (fixed)
                        beff = 0.7, # booster efficacy
                        wf=0.05, # protection for newly recoverd #0.2
                        N=5e6,
                        stngcy= 0.4,#0.78, #(*%(reduction)) strength of intervention (reduction in beta's)
                        eff_t = as.numeric(eff_date - intro_date)

)

parameters_BC <- parameters
init <- make_init()   #generate initial states.

# ---- run the model and compare the model to the data ----
times <- 1:(nrow(dat_omic) + forecasts_days)
allpars = c(parameters, beta_m=1 , p=0.5) 
outtest <- as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                                     parms=allpars)) 
inctest = get_total_incidence(outtest, parameters = allpars)
inctest$date = dat$date[1]+inctest$time
thisvec = c(test_prop, rep(test_prop[length(test_prop)], nrow(outtest)-length(test_prop))) 
inctest$inc_reported = inctest$inc_tot*thisvec

# sanity check - should have delta in a decline of about 2% /day (-0.02) and around 
# a 0.2 to 0.25 difference between the two, so omicron at about 0.2 
get_growth_rate(outtest, startoffset = 2, duration = 10)

ggplot(data =inctest, aes(x=date, y = inc_reported))+geom_line() +
  geom_line(aes(x=date, y= inc_res), color = "blue") +
    geom_line(aes(x=date, y= inc_mut), color = "red") +
  geom_point(data = dat, aes(x=date, y=cases), alpha=0.5) +
  ylim(c(0,50000)) + xlim(c(ymd("2021-11-20"), ymd("2022-02-28")))


# ---- fit the model  ---- CC -have not changed from here. JS: have updated from here (no param transforms), 
# but it doesn't work because of '1: In (x$Er + x$Erv + x$Erw + x$Em + x$Emv + x$Emw) * test_prop :
# longer object length is not a multiple of shorter object length' issue

#fitting beta_r, beta_m, p and dispersion parameter 
#guess <- c(log(0.7), logit(0.8),log(2.1),log(0.01))
# JS: Switched log and logit to match likelihood_func.R
guess <- c(1.5, 0.1)

#the parameters are constrained  accordingly (lower and upper)
rm(test_prop) #to check that it's being passed to LK

fit_BC <- optim(fn=func_loglik,  par=guess, lower=c(0, 0), 
                upper = c(Inf, 1), method = "L-BFGS-B", 
                test_prop=fake_test_prop_BC[1:nrow(dat_omic)], dat_omic=dat_omic)

# JS: testing out other ways to optimize:
# No bounds
#fit_BC <- optim(fn=func_loglik,  par=guess, method = "L-BFGS-B")
#function 'nlm' instead of optim
#fit_BC <- nlm(f=func_loglik,  p=guess, typsize=guess)

fit_BC 

#JS: Quick plotting likelihood surfaces
#z <- matrix(NA, 50,50)
#x <- c(seq(0, 5, length.out=50))
#y <- c(seq(0, 5, length.out=50))
#for (i in 1:50){
#  for (j in 1:50){
#   z[i,j] <- func_loglik(par=c(log(x[i]), log(y[j]), log(0.5), log(0.1)),test_prop,dat_omic)
#  }
#}
#image(x,y,z)


# beta_m <- seq(from=2.2,to=3.5,length=50)
# beta_r <- seq(0.9,3,length=50)
# grid <- expand.grid(x=beta_r,y=beta_m)
# 
# 
# grid <- ddply(grid,~x+y,mutate,loglik=func_loglik(par=c(x, y, logit(0.5),
#                                                 log(0.1)),test_prop,dat_omic))
# 
# 
# grid <- subset(grid,is.finite(loglik))
# ggplot(grid,aes(x=x,y=y,z=loglik,fill=loglik))+
#   geom_tile()+geom_contour(binwidth=1)
# 


#this catches estimated parameter values from MLE 
mle_est_BC <- c(beta_m=fit_BC$par[1], p=fit_BC$par[2],theta=0.1)

#check parameter estimates 
mle_est_BC

parameters <- c(parameters_BC,mle_est_BC)

#make prediction and projection with estimated parameters 

out_BC <- as.data.frame(deSolve::ode(y=init_BC,time=times,func= sveirs,
                                     parms=parameters)) 


#test_prop is shorter than projection, so we'll use the last value of test_prop for the rest of the simulation 
#test_prop_BC <- c(test_prop[1:length(dat_omic$value)], rep(last(test_prop),forecasts_days))


#with test_prop 
incidence_BC =  parameters[[1]]*(out_BC$Er + out_BC$Erv + out_BC$Erw +
                                   out_BC$Em + out_BC$Emv +
                                   out_BC$Emw)*fake_test_prop_BC



uncert_bound_BC = raply(simu_size,rnbinom(n=length(incidence_BC),
                                          mu=parameters[["p"]]*incidence_BC,
                                          size=1/parameters[["theta"]]))

project_dat_BC =  as.data.frame(aaply(uncert_bound_BC,2,quantile,
                                      na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
  mutate(date=seq.Date(ymd(intro_date),
                       ymd(intro_date)-1+length(times), 1))

#####################check fit ######################

ggplot() + geom_line(data=project_dat_BC,aes(x=date,y=`50%`), col="green",size=1.5,alpha=0.4) +
  geom_ribbon(data=project_dat_BC,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='darkgreen',alpha=0.1, size = 1.5)+
  geom_point(data=dat_reported,aes(x=date, y=value),color='grey48', alpha=0.8, size = 1.5)



##########################

#re-estimate parameters without test_prop 

guess_2  <- c(log(0.81), logit(0.8),log(2.1),log(0.01))

#the parameters are constrained  accordingly (lower and upper)

fit_BC_2 <- optim(fn=func_loglik_2,  par=guess_2, lower=c(log(0.6), log(2), 0.0001,  log(0.1)), 
                upper = c(log(0.8),log(2.5), 0.001,  log(0)), method = "L-BFGS-B")


#this catches estimated parameter values from MLE 
mle_est_BC_2 <- c(beta_r=exp(fit_BC$par[1]),beta_m=exp(fit_BC$par[2]), p=expit(fit_BC$par[3])
                ,theta=exp(fit_BC$par[4]))

#check parameter estimates 
mle_est_BC_2

parameters_2 <- c(parameters_BC,mle_est_BC_2)

#make prediction and projection with estimated parameters 

out_BC_2 <- as.data.frame(deSolve::ode(y=init_BC,time=times,func= sveirs,
                                     parms=parameters_2)) 




#without test_prop 
incidence_BC_rel =  parameters[[1]]*(out_BC_2$Er + out_BC_2$Erv + out_BC_2$Erw +
                                       out_BC_2$Em + out_BC_2$Emv +
                                       out_BC_2$Emw)

uncert_bound_BC_rel = raply(simu_size,rnbinom(n=length(incidence_BC_rel),
                                              mu=parameters[["p"]]*incidence_BC_rel,
                                              size=1/parameters[["theta"]]))

project_dat_BC_rel <- as.data.frame(aaply(uncert_bound_BC_rel,2
                                          ,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
  mutate(date=seq.Date(ymd(intro_date),
                       ymd(intro_date)-1+length(times), 1))



#add dat to data for plotting 
dat_reported <- dat_omic  %>% mutate(date=seq.Date(ymd(intro_date),
                              ymd(intro_date)-1+length(dat_omic$day), 1))





saveRDS(project_dat_BC, file.path("data/BC_test_constraints.rds"))
project_dat_BC =readRDS(file.path("data/BC_test_constraints.rds"))

saveRDS(project_dat_BC_rel, file.path("data/BC_no_constraints.rds"))
project_dat_BC_rel =readRDS(file.path("data/BC_no_constraints.rds"))


get_true_incidence_plot(times, start_date=intro_date, 
                        parameters_base=c(parameters_BC,mle_est_BC), init=init_BC)

get_true_incidence_prop_plot(times, start_date=intro_date, 
                             parameters_base=c(parameters_BC,mle_est_BC), init=init_BC)




cols <- c("Current, with testing constraints" = "darkgreen", "Without testing constraints"="orange")

gg_BC <- ggplot() + geom_line(data=project_dat_BC,aes(x=date,y=`50%`, colour = "Current, with testing constraints"),size=1.2,alpha=0.4) +
  geom_ribbon(data=project_dat_BC,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='darkgreen',alpha=0.1)+
  geom_line(data=project_dat_BC_rel,aes(x=date,y=`50%`, color="Without testing constraints"),size=1.2,alpha=0.4) +
  geom_ribbon(data=project_dat_BC_rel,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='orange',alpha=0.1)+
  geom_point(data=dat_reported,aes(x=date, y=value),color='grey48', alpha=0.8) + 
  #geom_line(aes(y=typical),color='blue') +
  labs(y="Reported cases",x="Date") + ylim(c(0,35000)) + 
  scale_x_date(date_breaks = "8 days", date_labels = "%b-%d-%y") +theme_light() +
  scale_color_manual(values = cols) +  theme(axis.text=element_text(size=15),
                                             plot.title = element_text(size=15, face="bold"),
                                             axis.text.x = element_text(angle = 45, hjust = 1),
                                             legend.position = "bottom", legend.title = element_text(size=15),
                                             legend.text = element_text(size=15),
                                             axis.title=element_text(size=15,face="bold")) +
  
  labs(color = " ",title="BC") +  geom_vline(xintercept=eff_date, linetype="dashed", 
                                             color = "grey", size=1)


gg_BC

ggsave(file="figs/BC_proj.png", gg_BC, width = 10, height = 8)
saveRDS(gg_BC, file.path("figs/BC-fig.rds"))

#things to try  
# fix beta_r and probably p 
# estimate beta_m alone 
#reduce population immunity and increase duration 




