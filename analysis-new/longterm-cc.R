require(deSolve)
require(ggplot2)
require(reshape2)
library(lubridate)
library(dplyr)
library(data.table)


source("analysis-new/mod_fitting_setup.R")
source("analysis-new/likelihood_func.R")

# --- -set up data 
forecasts_days <- 30 
intro_date <-  ymd("2021-11-20")
stop_date <- ymd("2022-01-24") # 
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
test_prop_BC1 <- c(test_prop_BC[1:length(dat_omic$value)], rep(last(test_prop_BC),forecasts_days)) #ensuring the length is consistent
fake_test_prop_BC <- (1 - (1-0.05)/(1 + exp(-0.25*(1:length(test_prop_BC1)-35))))

test_prop_BC <- fake_test_prop_BC[1:length(dat_omic$day)]
test_prop <- test_prop_BC 
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
init <- make_init()   #generate initial states. this function now uses the above 
# variables for its default input. 



#declaring fixed parameters 
eff_date <-   ymd("2021-12-29")  # intervention date 
mypars <-         c(sigma=1/3, # incubation period (3 days) (to fixed)
                        gamma=1/(5), #recovery rate (fixed)
                        nu =0.007, #vax rate: 0.7% per day (fixed)
                        mu=1/(82*365), # 1/life expectancy (fixed)
                        w1= 1/(0.5*365),# waning rate from R to S (fixed)
                        w2= 1/(0.5*365), # waning rate from Rv to V and Rw to W (fixed)
                        w3= 1/(0.5*365),# waning rate from boosted back to unboosted vax
                        ve=1, # I think this should be 1. it is not really efficacy  ( fixed)
                        beta_r=0.555, #transmission rate (to estimate) (0.35)
                        #beta_m=0.8*2.2, #transmission rate (to estimate)(*1.9)
                        epsilon_r = (1-0.8), # % this should be 1-ve 
                        epsilon_m = 1-0.3, #(1-0.25)?(1-0.6), # % escape capacity #(fixed)
                        b= 0.006, # booster rate  (fixed)
                        beff = 0.1, # booster efficacy
                        wf=0.15, # protection for newly recoverd #0.2
                        N=5e6,
                        stngcy= -0.1,#0.78, #(*%(reduction)) strength of intervention (reduction in beta's)
                        eff_t = as.numeric(eff_date - intro_date)
                        
)


pars1=c(mypars ,beta_m=1 , p=0.5) 

pars2=pars1
pars2["wf"]=0.3

pars2["w1"] = 1/365
pars2["w2"] = 1/365
pars2["w3"] = 1/365

gg = compare_two_parsets(pars1, pars2,
                         name1="Recovery 85% protection", 
                         name2 = "Recovery 70% projection", 
                         numdays =600,dispar=0.05) 
gg 

c <- 1 - mypars["stngcy"]/(1+ exp(-1.25*(1:300-mypars["eff_t"]))) 
plot(1:300, c)

compare_two_parsets = function(pars1, pars2,name1 = "first", name2="second", 
                               returnplot = T,numsamples = 1e3, numdays = 300, dispar=0.1) { 
    # run the first 
    out1 <- as.data.frame(deSolve::ode(y=init_BC,time=1:numdays,func= sveirs,
                                       parms=pars1)) 
    inc1 =  pars1["sigma"]*(out1$Er + out1$Erv + out1$Erw +
                                         out1$Em + out1$Emv +
                                         out1$Emw)
    unc1 = raply(numsamples,rnbinom(n=length(inc1),
                                              mu=pars1[["p"]]*inc1,
                                              size=1/dispar))
    
    proj1 =  as.data.frame(aaply(unc1,2,quantile,
                                          na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
        mutate(date=seq.Date(ymd(intro_date),
                             ymd(intro_date)-1+numdays, 1))
    proj1$name = name1
    # run the second 
    out2 <- as.data.frame(deSolve::ode(y=init_BC,time=1:numdays,func= sveirs,
                                       parms=pars2)) 
    inc2 =  pars2["sigma"]*(out2$Er + out2$Erv + out2$Erw +
                                 out2$Em + out2$Emv +
                                 out2$Emw)
    unc2 = raply(numsamples,rnbinom(n=length(inc2),
                                   mu=pars2[["p"]]*inc2,
                                   size=1/dispar))
    
    proj2 =  as.data.frame(aaply(unc2,2,quantile,
                                 na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
        mutate(date=seq.Date(ymd(intro_date),
                             ymd(intro_date)-1+numdays, 1))
    proj2$name = name2
    proj = rbind(proj1, proj2)
    if (returnplot) { 
        gg =   ggplot(data = proj) + geom_line(aes(x=date,y=`50%`, color=name),
                                               size=1.5,alpha=0.6) +
                        geom_ribbon(aes(x=date,ymin=`2.5%`,ymax=`97.5%`, fill=name),
                                    alpha=0.2, size = 1.5)+
            ylab("Incident infections") + theme_light() + theme(legend.position = "bottom")

    } else {
        return(proj) 
    }
}
