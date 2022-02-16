########## Source required files
source("analysis-new/functions.R")
source("analysis-new/likelihood_func.R")
#run BC_data.R line by line (if possible) :sometimes case data are incomplete, with 0's  and NA's 
source("analysis-new/BC_data.R") 




########## Data set up
forecasts_days <- 30 
intro_date <-  ymd("2021-11-30") #NOTE: now starting on Nov 30, for PHAC forecasts (it was 2021-11-20 )
stop_date <- ymd("2022-02-12")  
#import data 
dat = readRDS("data/BC-dat.rds")
#include Omicron wave only
dat <- dat %>% filter(date >= intro_date &  date <= stop_date)
dat_omic <- dat
dat_omic <- filter(dat_omic, date >= intro_date) %>% dplyr::select(c("day", "value", "date"))
dat_omic$day <- 1:nrow(dat_omic)



########## Test_prop set up
#-------EB: test_prop now has dates and this ensures it starts at the right date 

#tp_approx <- tp_approx_fit(mytest_BC=mytest_BC, dat= dat_omic, forecasts_days=forecasts_days,
#                          howlow = 0.2, 

tp_approx <- tp_approx_fit(mytest=mytest_BC, dat= dat_omic, 
                           forecasts_days=forecasts_days, howlow = 0.1, 
                           slope = 0.2,  midpoint=35, intro_date= intro_date, stop_date=stop_date)
plot_fit <- tp_approx[1]
plot_fit 

# test_prop for fit and projection 
test_prop_proj_BC <- data.frame(tp_approx[2])

#subset for fitting alone (length of data)
test_prop <- test_prop_proj_BC$tp[1:length(dat_omic$day)] #or set forecasts_days to 0 in tp_approx_fit() above 


# ---- cc's function to modify test_prop ----
# start lowering it even more when it hits the startlower value
# lower it by a linear function , ending at endfraction of its final value
# ie if endfraction is 0.5, the final value after adjustment is 1/2 the value before
# adjustment 
modifytp = function(startlower=1, endfraction=1, test_prop) {
  tt = min(which(test_prop < startlower))
  myline = seq(1, endfraction, length.out = length(test_prop)-tt)
  myfrac = c(rep(1, tt), myline)
  return(test_prop*myfrac)
}

mytp = modifytp(startlower = 1, endfraction = 1, test_prop)
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


#sources
#http://www.bccdc.ca/health-info/diseases-conditions/covid-19/data#COVID-19Dashboard 

dat %>% filter(date == intro_date )
########## Parameters set up 
#set values to generate initial conditions with make_init()
N=5.07e6
N_pop=N
vaxlevel_in = 0.88 # portion of the pop vaccinated at start time 
port_wane_in = 0.04 # portion boosted at start tie 
past_infection_in = 0.2 #  #increased this from 0.1 to 0.18 # total in R at start time
incres_in = 433*2 #( 389, reported cases on Nov 30) # resident strain (delta) incidence at start 
incmut_in = 6# new (omicron) inc at stat (#first cases of Omicron reported on 30 Nov)
simu_size = 1e5 # number of times to resample the negative binom (for ribbons)
forecasts_days =30 # how long to forecast for 
times = 1:nrow(dat_omic)

#declaring  parameters 
eff_date <-   ymd("2021-12-31")  # intervention date 
intv_date <-  ymd("2022-02-10")
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
                        epsilon_m = 1-0.15, # % 1-ve omicron 
                        b= 0.006, # booster rate
                        beff = 0.7, # booster efficacy
                        wf=0.1, # protection for newly recovered
                        N=5.07e6,
                        stngcy= 0.4, #(*%(reduction)) strength of intervention (reduction in beta's)
                        eff_t = as.numeric(eff_date - intro_date),
                        relx_level = 0,
                        rlx_t = as.numeric(intv_date - intro_date),
                        p = 0.25, #negative binomial mean
                        theta = 0.1 #negative binomial dispersion
                        
)
init <- make_init()   #generate initial states




########## Run the model and compare the model to the data
#          (i.e. without fitting anything)
times <- 1:(nrow(dat_omic) + forecasts_days)
outtest <- as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                                      parms=parameters)) 
inctest = get_total_incidence(outtest, parameters = parameters) # has p but not test_prop! 
inctest$date = intro_date-1+1:nrow(inctest)
thisvec=extendtp(nrow(inctest), test_prop)
inctest$inc_reported = inctest$inc_tot*thisvec

# sanity check - should have delta in a decline of about 2% /day (-0.02) and around 
# a 0.2 to 0.25 difference between the two, so omicron at about 0.2 
get_growth_rate(outtest, startoffset = 2, duration = 10)

ggplot(data =inctest, aes(x=date, y = inc_reported))+geom_line() +
  geom_line(aes(x=date, y= inc_res), color = "blue") +
  geom_line(aes(x=date, y= inc_mut), color = "red") +
  geom_point(data = dat, aes(x=date, y=cases), alpha=0.5) +
  ylim(c(0,5000)) + xlim(c(ymd("2021-11-20"), ymd("2022-02-28")))


#################### Alternative: penalized log likelihood
# Penalize loglh by (a) distance of %residentstrain from 50% on Dec 12th and
#                   (b) distance of mutant relative growth advantage from 0.2 during December 5th-15th
# according to penalty weight, pen.weight (a new input)

# Establish what the penalties are
# 1. 50/50 resident vs mutant on Dec 12th
known_prop <- 0.5
date_known_prop <- "2021-12-12"
# 2. growth advantage of mutant strain
known_growth <- 0.2
period_known_growth <- c("2021-12-05", "2021-12-15")
penalties <- list(known_prop = known_prop, date_known_prop = date_known_prop, 
                  known_growth = known_growth, period_known_growth = period_known_growth)

# Determine weight of penalty. 
# Qu: how strong should penalty be on scale of 0-1? 0 = no penalty. 1 = relatively as impactful as the likelihood
pen.size <- 0.6

guess <- c( beta_m=1, stngcy=0.4,beta_r=0.6, theta=0.1, p=0.2) 
pen.fit_BC <- optim(fn=func_penloglik,  par=guess, lower=c(0,0,0,0.001,0), upper = c(Inf,1,Inf,Inf,1), 
                    method = "L-BFGS-B", 
                    parameters = parameters, test_prop=test_prop, 
                    pen.size=pen.size, penalties = penalties, dat_omic=dat_omic, hessian=T)
pen.fit_BC
####################
func_loglik(pen.fit_BC$par, test_prop, dat_omic,parameters) 
parameters[names(guess)] <- pen.fit_BC$par


gg = simple_prev_plot(parameters, numdays = 190, mode = "both"); gg  # cc's simple prevalence plot 





########## Check fit: make prediction and projection with estimated parameters

#NOTE
#EB----------- use fitted test_prop instead of extendtp() for PHAC forecats 

mytp = mytp = modifytp(startlower = 1, endfraction = 1, 
                       test_prop = test_prop_proj_BC$tp[1:(length(dat_omic$day)+forecasts_days)])
test_prop <- mytp


################################# re-sampling CI

library(MASS)
times <- 1:(nrow(dat_omic) + forecasts_days)
# Sample from asymptotic distribution of MLEs, 100 times
resampled <- mvrnorm(n = 100, mu = pen.fit_BC$par, Sigma = (1/nrow(dat_omic))*solve(pen.fit_BC$hessian)) # don't need -ve because we MINimised NEGloglike

resample_incidence <- function(x, parameters){
  # For each set of resampled MLEs, do a model simulation and output the daily incident reported cases
  parameters[names(guess)] <- x
  out_samp <- as.data.frame(deSolve::ode(y=init,time=times,func=sveirs,
                                         parms=parameters)) 
  incidence_samp =  parameters[["p"]]*parameters[["sigma"]]*(out_samp$Er + out_samp$Erv + out_samp$Erw +
                                                               out_samp$Em + out_samp$Emv +
                                                               out_samp$Emw)*test_prop#extendtp(nrow(out_samp), test_prop)
}
incidence_resampled <- apply(resampled, 1, resample_incidence, parameters=parameters)


bound_mlesample_BC <-apply(incidence_resampled,2,  function(x){raply(simu_size/100,rnbinom(n=length(x),
                                                           mu=x, size=1/parameters[["theta"]]))}, simplify=FALSE)
bound_mlesample_BC <- do.call(rbind, bound_mlesample_BC)

project_dat_BC  =  as.data.frame(aaply(bound_mlesample_BC 
                                         ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
  mutate(date=seq.Date(ymd(intro_date),
                       ymd(intro_date)-1+length(times), 1))


#add dat to data for plotting 

dat_reported <- dat_omic  %>% mutate(date=seq.Date(ymd(intro_date),
                                                   ymd(intro_date)-1+length(dat_omic$day), 1))

ggplot() + geom_line(data=project_dat_BC,aes(x=date,y=`50%`), col="green",size=1.5,alpha=0.4) +
  geom_ribbon(data=project_dat_BC,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='darkgreen',alpha=0.1, size = 1.5)+
  geom_point(data=dat_reported,aes(x=date, y=value),color='grey48', alpha=0.8, size = 1.5)

# ---- end check fit ---- 

 
mu1 <- mean(dat_omic$value)
mu2 <- mu1/last(test_prop)

theta_2 = 1 / ( 1/mu1 - 1/mu2 + 1/parameters[["theta"]])  #m1 / end of test_prop

parameters[["theta"]] <- theta_2 

times <- 1:(nrow(dat_omic) + forecasts_days)
out_BC_rel <- as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                                         parms=parameters)) 

incidence_BC_rel =  parameters[["p"]]*parameters[["sigma"]]*(out_BC_rel$Er + out_BC_rel$Erv + out_BC_rel$Erw +
                                                               out_BC_rel$Em + out_BC_rel$Emv +
                                                               out_BC_rel$Emw)

uncert_bound_BC_rel = raply(simu_size,rnbinom(n=length(incidence_BC_rel),
                                              mu=incidence_BC_rel,
                                              size=1/parameters[["theta"]]))
project_dat_BC_rel =  as.data.frame(aaply(uncert_bound_BC_rel 
                                          ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
mutate(date=seq.Date(ymd(intro_date),
                       ymd(intro_date)-1+length(times), 1))

plot(project_dat_BC_rel$`50%`)


#########################################vary parameters  ##########
#issues--- can't re-estimate without tests prop... parameters become unreasonable 

relx_level = data.frame(relx_level=c(0.4,0.8))

vary_parameter <- function(x, parameters,init1=init){
  parameters["relx_level"] <- x
  out_var <- as.data.frame(deSolve::ode(y=init1,time=times,func=sveirs,
                                         parms=parameters)) 
  incidence_var =  parameters[["p"]]*parameters[["sigma"]]*(out_var$Er + out_var$Erv + out_var$Erw +
                                                               out_var$Em + out_var$Emv +
                                                               out_var$Emw)
}
incidence_var <- data.frame(apply(relx_level, 1, vary_parameter, parameters=parameters))
plot(NA,NA, xlim=c(0,100), ylim=c(0,50000))
lines(incidence_BC_rel)
lines(incidence_var$X1)
lines(incidence_var$X2)



uncert_bound_BC_rel_1 = raply(simu_size,rnbinom(n=length(incidence_var$X1),
                                              mu=incidence_var$X1,
                                              size=1/parameters[["theta"]]))

project_dat_BC_rel_1 =  as.data.frame(aaply(uncert_bound_BC_rel_1 
                                          ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
                                           mutate(date=seq.Date(ymd(intro_date),
                                           ymd(intro_date)-1+length(times), 1))



uncert_bound_BC_rel_2 = raply(simu_size,rnbinom(n=length(incidence_var$X2),
                                                mu=incidence_var$X2,
                                                size=1/parameters[["theta"]]))

project_dat_BC_rel_2 =  as.data.frame(aaply(uncert_bound_BC_rel_2 
                                            ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
  mutate(date=seq.Date(ymd(intro_date),
                       ymd(intro_date)-1+length(times), 1))





##################


saveRDS(project_dat_BC, file.path("data/BC_test_constraints.rds"))
project_dat_BC =readRDS(file.path("data/BC_test_constraints.rds"))

saveRDS(project_dat_BC_rel, file.path("data/BC_no_constraints.rds"))
project_dat_BC_rel =readRDS(file.path("data/BC_no_constraints.rds"))

saveRDS(project_dat_BC_rel_1, file.path("data/BC_40percent_inc.rds"))
project_dat_BC_rel_1 = readRDS(file.path("data/BC_40percent_inc.rds"))

saveRDS(project_dat_BC_rel_2, file.path("data/BC_80percent_inc.rds"))
project_dat_BC_rel_2 = readRDS(file.path("data/BC_80percent_inc.rds"))
# "#", 
#make figures 

cols <- c("Current, TC" = "darkgreen",
          "NTC"="orange", "40% increase"="#0072B2", "80% increase"= "#CC79A7")

gg_BC <- ggplot() + geom_line(data=project_dat_BC,aes(x=date,y=`50%`, colour = "Current, TC"),size=1.2,alpha=0.4) +
  geom_ribbon(data=project_dat_BC,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='darkgreen',alpha=0.1)+
  geom_line(data=project_dat_BC_rel,aes(x=date,y=`50%`, color="NTC"),size=1.2,alpha=0.4) +
  geom_ribbon(data=project_dat_BC_rel,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='orange',alpha=0.1)+
  geom_line(data=project_dat_BC_rel_1,aes(x=date,y=`50%`, color="Without testing constraints"),size=1.2,alpha=0.4) +
  geom_ribbon(data=project_dat_BC_rel_1,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='yellow',alpha=0.1)+
  geom_point(data=dat_reported,aes(x=date, y=value),color='grey48', alpha=0.8) + 
  
  geom_line(data=project_dat_BC_rel_1,aes(x=date,y=`50%`, color="40% increase"),size=1.2,alpha=0.4) +
  geom_ribbon(data=project_dat_BC_rel_1,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill="purple",alpha=0.1)+
  
  geom_line(data=project_dat_BC_rel_2,aes(x=date,y=`50%`, color="80% increase"),size=1.2,alpha=0.4) +
  geom_ribbon(data=project_dat_BC_rel_2,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill="#D55E00",alpha=0.1)+
  
  #geom_line(aes(y=typical),color='blue') +
  labs(y="Reported cases",x="Date") + ylim(c(0,30000)) + 
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

ggsave(file="figs/BC_proj.png", gg_BC, width = 14, height = 8)
saveRDS(gg_BC, file.path("figs/BC-fig.rds"))






