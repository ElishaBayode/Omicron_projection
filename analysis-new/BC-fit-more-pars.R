########## Source required files
source("analysis-new/functions.R")
source("analysis-new/likelihood_func.R")
#run BC_data.R line by line (if possible) :sometimes case data are incomplete, with 0's  and NA's 
source("analysis-new/BC_data.R") 




########## Data set up
forecasts_days <- 30 
intro_date <-  ymd("2021-11-20")
stop_date <- ymd("2022-01-30")  
#import data 
dat = readRDS("data/BC-dat.rds")
#include Omicron wave only
dat <- dat %>% filter(date >= intro_date &  date <= stop_date)
dat_omic <- dat
dat_omic <- filter(dat_omic, date >= intro_date) %>% select(c("day", "value", "date"))
dat_omic$day <- 1:nrow(dat_omic)



########## Test_prop set up
#-------EB: test_prop now has dates and this ensures it starts at the right date 

#tp_approx <- tp_approx_fit(mytest_BC=mytest_BC, dat= dat_omic, forecasts_days=forecasts_days,
#                          howlow = 0.2, 

tp_approx <- tp_approx_fit(mytest=mytest_BC, dat= dat_omic, 
                           forecasts_days=forecasts_days, howlow = 0.1, 
                           slope = 0.2,  midpoint=45, intro_date= intro_date, stop_date=stop_date)
plot_fit <- tp_approx[1]
plot_fit 

# test_prop for fit and projection 
test_prop_proj <- data.frame(tp_approx[2])

#subset for fitting alone (length of data)
test_prop <- test_prop_proj$tp[1:length(dat_omic$day)] #or set forecasts_days to 0 in tp_approx_fit() above 


# ---- cc's function to modify test_prop ----
# start lowering it even more when it hits the startlower value
# lower it by a linear function , ending at endfraction of its final value
# ie if endfraction is 0.5, the final value after adjustment is 1/2 the value before
# adjustment 
modifytp = function(startlower=0.6, endfraction=0.55, test_prop) {
  tt = min(which(test_prop < startlower))
  myline = seq(1, endfraction, length.out = length(test_prop)-tt)
  myfrac = c(rep(1, tt), myline)
  return(test_prop*myfrac)
}

mytp = modifytp(startlower = 0.5, endfraction = 0.35, test_prop)
ggplot(data = data.frame(date = intro_date+1:length(test_prop), 
                         test_prop=test_prop, 
                         mytp = mytp), 
       aes(x=date, y=test_prop))+geom_line() + 
  geom_line(aes(x=date, y=mytp), color="blue") + ylim(c(0,1))

test_prop = mytp 
# NOTE leaving as is for now 


# --- cc's function to append values to test_prop to stop getting length errors ----
# see changes to how test_prop is used in the likelihood functions
# basically if test_prop does not start at the same date as the simulation, you are hooped 
extendtp <- function(n=100, test_prop=test_prop){
  return(c(test_prop, rep(test_prop[length(test_prop)], n-length(test_prop))))
}

simu_size = 1e5 # number of times to resample the negative binom (for ribbons)
forecasts_days =30 # how long to forecast for 
times = 1:nrow(dat_omic)

#declaring  parameters 
eff_date <-   ymd("2021-12-31")  # intervention date 
parameters <-         c(sigma=1/3, # incubation period (days) 
                        gamma=1/5, #recovery rate 
                        nu =0.007, #vax rate: 0.7% per day 
                        mu=1/(82*365), # 1/life expectancy 
                        w1= 1/(0.5*365),# waning rate from R to S 
                        w2= 1/(0.5*365), # waning rate from Rv to V 
                        w3= 1/(0.5*365),# waning rate Rw to W 
                        ve=1, # I think this should be 1. it is not really efficacy  
                        beta_r=0.6, #transmission rate 
                        beta_m=1, #transmission rate 
                        epsilon_r = (1-0.8), # % this should be 1-ve 
                        epsilon_m = 1-0.3, # % 1-ve omicron 
                        b = 0.006, # booster rate
                        beff = 0.7, # booster efficacy
                        wf=0.05, # protection for newly recovered
                        N=5e6,
                        stngcy= 0.4, #(*%(reduction)) strength of intervention (reduction in beta's)
                        eff_t = as.numeric(eff_date - intro_date),
                        p = 0.5, #negative binomial mean
                        theta = 0.1 #negative binomial dispersion
                        
)

init <- make_init(N=5.07e6, vaxlevel = 0.88,
                  port_wane = 0.04 , 
                  past_infection = 0.18, incres = 470, incmut =  6, 
                  pars=as.list(parameters))



########## Fit  more parameters (epsilon_m, epsilon_r, ve sigma are used in make_init() to generate initial states)
# REMEMBER: params must be named and need to ensure there's an upper and lower bound for each in optim

guess <- c( beta_m=1.2, beta_r=0.6, epsilon_m=0.7, epsilon_r =0.2,  sigma=1/3, stngcy =0.4) 
#guess <- c( beta_m=1.2, beta_r=0.6, epsilon_m=0.7,
#            epsilon_r =0.2,  sigma=0.3) ,sigma estimated to be 0.197, add p=0.3 (p estimated as 1)

#beta_m    beta_r epsilon_m epsilon_r     sigma     gamma    stngcy 
#1.0101158 0.5433433 0.5221261 0.1608351 0.3107221 0.1295124 0.3426588  

#the parameters are constrained  accordingly (lower and upper)
fit_BC <- optim(fn=func_loglik,  par=guess, lower=c(0,0,0,0,0,0), #if theta's range includes 0, it generates NA's
                upper = c(Inf, Inf,1,1, 1,1), method = "L-BFGS-B", parameters = parameters,
                test_prop=test_prop, dat_omic=dat_omic, init = init)
# check the values:
fit_BC
func_loglik(fit_BC$par, test_prop, dat_omic,parameters) 






#this catches estimated parameter values from MLE , and adds them to 'parameters' structure
parameters[names(guess)] <- fit_BC$par

plot.loglik.info(parameters, 1:nrow(dat_omic), test_prop,init=init) # cc's sanity check plot 


#use estimated params to generate make_init()

init_estim <- make_init(N=5.07e6, vaxlevel = 0.88,
                                port_wane = 0.04 , 
                                past_infection = 0.18, incres = 470, incmut =  6, 
                                pars=as.list(parameters))



guess2 <- c(beta_m=1.2, beta_r=0.6, stngcy = 0.4, theta=0.1, gamma=1/5 ) 

#the parameters are constrained  accordingly (lower and upper)
fit_BC2 <- optim(fn=func_loglik,  par=guess2, lower=c(0,0,0,0.001,-Inf,0), #if theta range includes 0, it generates NA's
                upper = c(Inf,Inf,1, Inf,Inf,1), method = "L-BFGS-B", parameters = parameters,
                test_prop=test_prop, dat_omic=dat_omic, init = init_estim)
# check the values:
fit_BC2
func_loglik(fit_BC2$par, test_prop, dat_omic,parameters) 



parameters[names(guess2)] <- fit_BC2$par

plot.loglik.info(params=parameters, 1:nrow(dat_omic), test_prop,init = init_estim)

gg = simple_prev_plot(parameters, numdays = 190, mode = "both"); gg  # cc's simple prevalence plot 



########## Check fit: make prediction and projection with estimated parameters
times <- 1:(nrow(dat_omic) + forecasts_days)
out_BC <- as.data.frame(deSolve::ode(y=init_estim,time=times,func= sveirs,
                                     parms=parameters)) 
#with test_prop 
incidence_BC =  parameters[["p"]]*parameters[["sigma"]]*(out_BC$Er + out_BC$Erv + out_BC$Erw +
                                                           out_BC$Em + out_BC$Emv +
                                                           out_BC$Emw)*extendtp(nrow(out_BC), test_prop)

plot(1:nrow(dat_omic), incidence_BC[1:nrow(dat_omic)])
uncert_bound_BC = raply(simu_size,rnbinom(n=length(incidence_BC),
                                          mu=incidence_BC,
                                          size=1/parameters[["theta"]]))

project_dat_BC =  as.data.frame(aaply(uncert_bound_BC 
                                      ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
  mutate(date=seq.Date(ymd(intro_date),
                       ymd(intro_date)-1+length(times), 1))
#add dat to data for plotting 

dat_reported <- dat_omic  %>% mutate(date=seq.Date(ymd(intro_date),
                                                   ymd(intro_date)-1+length(dat_omic$day), 1))

ggplot() + geom_line(data=project_dat_BC,aes(x=date,y=`50%`), col="green",size=1.5,alpha=0.4) +
  geom_ribbon(data=project_dat_BC,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='darkgreen',alpha=0.1, size = 1.5)+
  geom_point(data=dat_reported,aes(x=date, y=value),color='grey48', alpha=0.8, size = 1.5)


