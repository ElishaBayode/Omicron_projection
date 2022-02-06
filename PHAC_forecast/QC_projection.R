

#set values to generate initial conditions with make_init()

source("analysis-new/functions.R")
source("analysis-new/likelihood_func.R")
#run QC_data.R line by line (if possible) :sometimes case data are incomplete, with 0's  and NA's 
source("analysis-new/QC_data.R") 

#sources
#https://angusreid.org/omicron-incidence-restrictions/ (interesting stats)
#https://www.quebec.ca/en/health/health-issues/a-z/2019-coronavirus/situation-coronavirus-in-quebec/covid-19-vaccination-data

#dat %>% filter(`Date reported` == ymd("2021-11-30"))
########## Data set up
forecasts_days <- 30 
intro_date <-  ymd("2021-11-30")
stop_date <- ymd("2022-01-30")  
#import data 
dat = readRDS("data/QC-dat.rds")
#include Omicron wave only
#dat <- dat %>% filter(date >= intro_date &  date <= stop_date)
dat_omic <- dat
dat_omic <- filter(dat_omic, date >= intro_date) %>% select(c("day", "value", "date"))
dat_omic$day <- 1:nrow(dat_omic)



########## Test_prop set up
#-------EB: test_prop now has dates and this ensures it starts at the right date 
tp_approx <- tp_approx_fit(mytest=mytest_QC, dat= dat_omic, forecasts_days=forecasts_days, howlow = 0.1, 
                           slope = 0.12,  midpoint=34, intro_date= intro_date, stop_date=stop_date)
plot_fit <- tp_approx[1]
plot_fit 

# test_prop for fit and projection 
test_prop_proj <- data.frame(tp_approx[2])

#subset for fitting alone (length of data)
#or set forecasts_days to 0 in tp_approx_fit() above 

modifytp = function(startlower=0.6, endfraction=0.5, test_prop=test_prop_proj$tp[1:length(dat_omic$day)]) {
  tt = min(which(test_prop < startlower))
  myline = seq(1, endfraction, length.out = length(test_prop)-tt)
  myfrac = c(rep(1, tt), myline)
  return(test_prop*myfrac)
}

### EB____ we should modify test_prop within modifytp()

mytp = modifytp(startlower = 0.5, endfraction = 0.4, test_prop)
ggplot(data = data.frame(date = intro_date+1:length(test_prop), 
                         test_prop=test_prop, 
                         mytp = mytp), 
       aes(x=date, y=test_prop))+geom_line() + 
  geom_line(aes(x=date, y=mytp), color="blue") + ylim(c(0,1))

test_prop = mytp

# --- cc's function to append values to test_prop to stop getting length errors ----

extendtp <- function(n=100, test_prop=test_prop){
  return(c(test_prop, rep(test_prop[length(test_prop)], n-length(test_prop))))
}

filter(dat, date == ymd("2021-11-30") )

########## Parameters set up 
#set values to generate initial conditions with make_init()

N=8.485e6
N_pop=N
vaxlevel_in =  0.76 # 
port_wane_in = 0.05 # portion boosted at start time .... on Nov 30
past_infection_in = 0.2  #increased this from 0.1 to 0.18 # total in R at start time (620K so far, with 50% asc rate(guess it's lower in QC))
incres_in = 784  ## # resident strain (delta) incidence at start (784 , reported cases on Nov 30) 
incmut_in = 25 # new (omicron) inc at stat (presumably very low)---28% in QC as of Dec 14, frist detected on 28th Nov 
simu_size = 1e5 # number of times to resample the negative binom (for ribbons)
forecasts_days =30 # how long to forecast for 
times = 1:nrow(dat_omic)

#declaring  parameters 
eff_date <-   ymd("2021-12-31")  # intervention date # restrictions were updated on 24 December 
parameters <-         c(sigma=1/3, # incubation period (days) 
                        gamma=1/(5), #recovery rate 
                        nu =0.007, #vax rate: 0.7% per day 
                        mu=1/(82*365), # 1/life expectancy 
                        w1= 1/(0.5*365),# waning rate from R to S 
                        w2= 1/(0.5*365), # waning rate from Rv to V 
                        w3= 1/(0.5*365),# waning rate Rw to W 
                        ve=1, # I think this should be 1. it is not really efficacy  
                        beta_r=0.3, #transmission rate 
                        beta_m=1.2, #transmission rate 
                        epsilon_r = (1-0.8), # % this should be 1-ve 
                        epsilon_m = 1-0.3, # % escape capacity 
                        b= 0.006, # booster rate
                        beff = 0.7, # booster efficacy
                        wf=0.1, # protection for newly recovered
                        N=8.485e6,
                        stngcy= 0.4, #(*%(reduction)) strength of intervention (reduction in beta's)
                        eff_t = as.numeric(eff_date - intro_date),
                        p = 0.5, #negative binomial mean
                        theta = 0.1 #negative binomial dispersion
                        
)
init <- make_init()   #generate initial states




########## Run the model and compare the model to the data
#          (i.e. without fitting anything)
times <- 1:(nrow(dat_omic) + forecasts_days)
outtest <- as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                                      parms=parameters)) 
inctest = get_total_incidence(outtest, parameters = parameters)
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
  ylim(c(0,30000)) + xlim(c(ymd("2021-11-20"), ymd("2022-02-28")))



########## Fit the model
# fitting any subset of parameters. 
# REMEMBER: params must be named and need to ensure there's an upper and lower bound for each in optim
guess <- c(beta_m = 1.1, stngcy=0.5, beta_r=0.6) 


#the parameters are constrained  accordingly (lower and upper)
fit_QC <- optim(fn=func_loglik,  par=guess, lower=c(0,0,0,-Inf), 
                upper = c(Inf,1,Inf, Inf), method = "L-BFGS-B", parameters = parameters,
                test_prop=test_prop, dat_omic=dat_omic)
# check the values:
fit_QC
func_loglik(fit_QC$par, test_prop, dat_omic,parameters) 

#this catches estimated parameter values from MLE , and adds them to 'parameters' structure
parameters[names(guess)] <- fit_QC$par
plot.loglik.info(parameters, 1:nrow(dat_omic), test_prop) # cc's sanity check plot 
gg_QC = simple_prev_plot(parameters, numdays = 190) + ylab("Prevalence")  + 
  labs(title="QC"); gg_QC  # cc's simple prevalence plot 
ggsave(file="figs/QC_proj.png", gg_QC, width = 8, height = 8)

#----------- use fitted test_prop instead of extendtp()

mytp = mytp = modifytp(startlower = 0.5, endfraction = 0.4, 
                       test_prop = test_prop_proj$tp[1:(length(dat_omic$day)+forecasts_days)])
test_prop <- mytp




########## Check fit: make prediction and projection with estimated parameters
times <- 1:(nrow(dat_omic) + forecasts_days)
out_QC <- as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                                     parms=parameters)) 
#with test_prop 
incidence_QC =  parameters[["p"]]*parameters[["sigma"]]*(out_QC$Er + out_QC$Erv + out_QC$Erw +
                                                           out_QC$Em + out_QC$Emv +
                                                           out_QC$Emw)*test_prop #extendtp(nrow(out_QC), test_prop)

plot(1:nrow(dat_omic), incidence_QC[1:nrow(dat_omic)])
uncert_bound_QC = raply(simu_size,rnbinom(n=length(incidence_QC),
                                          mu=incidence_QC,
                                          size=1/parameters[["theta"]]))
project_dat_QC =  as.data.frame(aaply(uncert_bound_QC 
                                      ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
  mutate(date=seq.Date(ymd(intro_date),
                       ymd(intro_date)-1+length(times), 1))
#add dat to data for plotting 
dat_reported <- dat_omic  %>% mutate(date=seq.Date(ymd(intro_date),
                                                   ymd(intro_date)-1+length(dat_omic$day), 1))



ggplot() + geom_line(data=project_dat_QC,aes(x=date,y=`50%`), col="green",size=1.5,alpha=0.4) +
  geom_ribbon(data=project_dat_QC,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='darkgreen',alpha=0.1, size = 1.5)+
  geom_point(data=dat_reported,aes(x=date, y=value),color='grey48', alpha=0.8, size = 1.5)




#####____ fit and make projections without test_prop 

guess <- c(beta_m = 1.1, stngcy=0.5, beta_r=0.6) 


fit_QC_rel <- optim(fn=func_loglik,  par=guess, lower=c(0,0,0), 
                    upper = c(Inf,1,Inf), method = "L-BFGS-B", parameters = parameters,
                    test_prop=rep(1,nrow(dat_omic)), dat_omic=dat_omic) #set test_prop=1
# check the values:

parameters[names(guess)] <- fit_QC_rel$par

times <- 1:(nrow(dat_omic) + forecasts_days)
out_QC_rel <- as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                                         parms=parameters)) 

incidence_QC_rel =  parameters[["p"]]*parameters[["sigma"]]*(out_QC$Er + out_QC$Erv + out_QC$Erw +
                                                               out_QC$Em + out_QC$Emv +
                                                               out_QC$Emw)

uncert_bound_QC_rel = raply(simu_size,rnbinom(n=length(incidence_QC_rel),
                                              mu=incidence_QC_rel,
                                              size=1/parameters[["theta"]]))
project_dat_QC_rel =  as.data.frame(aaply(uncert_bound_QC_rel 
                                          ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
  mutate(date=seq.Date(ymd(intro_date),
                       ymd(intro_date)-1+length(times), 1))



saveRDS(project_dat_QC, file.path("data/QC_test_constraints.rds"))
project_dat_QC =readRDS(file.path("data/QC_test_constraints.rds"))

saveRDS(project_dat_QC_rel, file.path("data/QC_no_constraints.rds"))
project_dat_QC_rel =readRDS(file.path("data/QC_no_constraints.rds"))



#make figures 

cols <- c("Current, with testing constraints" = "darkgreen", "Without testing constraints"="orange")

gg_QC <- ggplot() + geom_line(data=project_dat_QC,aes(x=date,y=`50%`, colour = "Current, with testing constraints"),size=1.2,alpha=0.4) +
  geom_ribbon(data=project_dat_QC,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='darkgreen',alpha=0.1)+
  geom_line(data=project_dat_QC_rel,aes(x=date,y=`50%`, color="Without testing constraints"),size=1.2,alpha=0.4) +
  geom_ribbon(data=project_dat_QC_rel,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='orange',alpha=0.1)+
  geom_point(data=dat_reported,aes(x=date, y=value),color='grey48', alpha=0.8) + 
  #geom_line(aes(y=typical),color='blue') +
  labs(y="Reported cases",x="Date") + ylim(c(0,70000)) + 
  scale_x_date(date_breaks = "8 days", date_labels = "%b-%d-%y") +theme_light() +
  scale_color_manual(values = cols) +  theme(axis.text=element_text(size=15),
                                             plot.title = element_text(size=15, face="bold"),
                                             axis.text.x = element_text(angle = 45, hjust = 1),
                                             legend.position = "bottom", legend.title = element_text(size=15),
                                             legend.text = element_text(size=15),
                                             axis.title=element_text(size=15,face="bold")) +
  
  labs(color = " ",title="QC") +  geom_vline(xintercept=eff_date, linetype="dashed", 
                                             color = "grey", size=1)


gg_QC

ggsave(file="figs/QC_proj.png", gg_QC, width = 10, height = 8)
saveRDS(gg_QC, file.path("figs/QC-fig.rds"))

