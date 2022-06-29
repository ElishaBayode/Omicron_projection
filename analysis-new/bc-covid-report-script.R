
########## cc plot for bc covid group, taken from bc projectoin 


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



########## Parameters set up 
#set values to generate initial conditions with make_init()
N=5.07e6
N_pop=N
vaxlevel_in = 0.88 # portion of the pop vaccinated at start time 
port_wane_in = 0.04 # portion boosted at start tie 
past_infection_in = 0.18  #increased this from 0.1 to 0.18 # total in R at start time
incres_in = 650 # resident strain (delta) incidence at start 
incmut_in = 15# new (omicron) inc at stat 
simu_size = 1e5 # number of times to resample the negative binom (for ribbons)
forecasts_days =30 # how long to forecast for 
times = 1:nrow(dat_omic)

#declaring  parameters 
eff_date <-   ymd("2021-12-31")  # intervention date 
parameters <-         c(sigma=1/3, # incubation period (days) 
                        gamma=1/5, #recovery rate 
                        nu =0.007, #vax rate: 0.7% per day 
                        mu=1/(82*365), # 1/life expectancy 
                        w1= 1/(0.25*365),# waning rate from R to S 
                        w2= 1/(0.5*365), # waning rate from Rv to V 
                        w3= 1/(0.5*365),# waning rate Rw to W 
                        ve=1, # I think this should be 1. it is not really efficacy  
                        beta_r=0.6, #transmission rate 
                        beta_m=1, #transmission rate 
                        epsilon_r = (1-0.8), # % this should be 1-ve 
                        epsilon_m = 1-0.3, # % 1-ve omicron 
                        b= 0.006, # booster rate
                        beff = 0.7, # booster efficacy
                        wf=0.15, # protection for newly recovered
                        N=5e6,
                        stngcy= 0.4, #(*%(reduction)) strength of intervention (reduction in beta's)
                        eff_t = as.numeric(eff_date - intro_date),
                        p = 0.25, #negative binomial mean
                        theta = 0.1 #negative binomial dispersion
                        
)
init <- make_init()   #generate initial states



########## Fit the model
# fitting any subset of parameters. 
# REMEMBER: params must be named and need to ensure there's an upper and lower bound for each in optim
guess <- c( beta_m=1, stngcy=0.4,beta_r=0.6) 

#the parameters are constrained  accordingly (lower and upper)
fit_BC <- optim(fn=func_loglik,  par=guess, lower=c(0,0,0), 
                upper = c(Inf,1,Inf), method = "L-BFGS-B", parameters = parameters,
                test_prop=test_prop, dat_omic=dat_omic)
# check the values:
fit_BC
func_loglik(fit_BC$par, test_prop, dat_omic,parameters) 

#this catches estimated parameter values from MLE , and adds them to 'parameters' structure
parameters[names(guess)] <- fit_BC$par
plot.loglik.info(parameters, 1:nrow(dat_omic), test_prop) # cc's sanity check plot 
gg = simple_prev_plot(parameters, numdays = 190, mode = "both"); gg  # cc's simple prevalence plot 

# --- make plot of reported cases 
times <- 1:(nrow(dat_omic) + forecasts_days)
out_BC <- as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                                     parms=parameters)) 
#with test_prop 
incidence_BC =  parameters[["p"]]*parameters[["sigma"]]*(out_BC$Er + out_BC$Erv + out_BC$Erw +
                                                           out_BC$Em + out_BC$Emv +
                                                           out_BC$Emw)*test_prop_proj$tp
uncert_bound_BC = raply(simu_size,rnbinom(n=length(incidence_BC),
                                          mu=incidence_BC,
                                          size=1/parameters[["theta"]]))

project_dat_BC =  as.data.frame(aaply(uncert_bound_BC 
                                      ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975)))
project_dat_BC$date =  intro_date-1 + 1:nrow(project_dat_BC)
#add dat to data for plotting 
dat_reported <- dat_omic  %>% mutate(date=intro_date-1 + 1:nrow(dat_omic))

ggplot() + geom_line(data=filter(project_dat_BC, date < ymd("2022-02-02")),
                     aes(x=date,y=`50%`), col="green",size=1.5,alpha=0.4) +
  geom_ribbon(data=filter(project_dat_BC, date < ymd("2022-02-02")),aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='darkgreen',alpha=0.1, size = 1.5)+
  geom_point(data=filter(dat_reported, date < ymd("2022-02-02")),aes(x=date, y=value),color='grey48', alpha=0.8, size = 1.5)+
   theme_light() +theme(axis.title.x = element_blank())+ylab("Reported cases")


# make a plot of incidence without test prop and p 
# --- make plot of reported cases 
#with test_prop 
true_incidence_BC =  parameters[["sigma"]]*(out_BC$Er + out_BC$Erv + out_BC$Erw +
                                                           out_BC$Em + out_BC$Emv +
                                                           out_BC$Emw)
uncert_bound_BC = raply(simu_size,rnbinom(n=length(true_incidence_BC),
                                          mu=true_incidence_BC,
                                          size=1/parameters[["theta"]]))

project_dat_BC =  as.data.frame(aaply(uncert_bound_BC, 2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
  mutate(date=seq.Date(ymd(intro_date),
                       ymd(intro_date)-1+length(times), 1))
#add dat to data for plotting 
dat_reported <- dat_omic  %>% mutate(date=seq.Date(ymd(intro_date),
                                                   ymd(intro_date)-1+length(dat_omic$day), 1))

ggplot() + geom_line(data=project_dat_BC,aes(x=date-5,y=`50%`), col="blue",size=1.5,alpha=0.4) +
  geom_ribbon(data=project_dat_BC,aes(x=date-5,ymin=`2.5%`,ymax=`97.5%`),fill='blue',alpha=0.1, size = 1.5)+
  theme_light() +theme(axis.title.x = element_blank())+ylab("Infections")

