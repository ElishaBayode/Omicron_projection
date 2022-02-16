
source("analysis-new/functions.R")
source("analysis-new/likelihood_func.R")
#run SK_data.R line by line (if possible) :sometimes case data are incomplete, with 0's  and NA's 
source("PHAC_forecast/SK_data.R") 

#sources
#https://www.luminultra.com/covid-19-testing/wastewater-testing/?creative=550110891097&keyword=wastewater%20monitoring&matchtype=b&network=g&device=c&gclid=Cj0KCQiA9OiPBhCOARIsAI0y71DRpb51-NjOilN22prf0GoKZLgldkdRlJGmOnKMRz2NYca0-YQC2NIaAjL7EALw_wcB
#https://water.usask.ca/covid-19/#PrinceAlbertWastewaterData
#https://water.usask.ca/covid-19/#NorthBattlefordWastewaterData 
#https://nccid.ca/wastewater-surveillance-for-covid-19/ 

#dat %>% filter(`Date reported` == ymd("2021-11-30"))
########## Data set up
forecasts_days <- 30 
intro_date <-  ymd("2021-11-30")
stop_date <- ymd("2022-02-07")  
#import data 
dat = readRDS("data/SK-dat.rds")
#include Omicron wave only
dat <- dat %>% filter(date >= intro_date &  date <= stop_date)
dat_omic <- dat
dat_omic <- filter(dat_omic, date >= intro_date) %>% dplyr::select(c("day", "value", "date"))
dat_omic$day <- 1:nrow(dat_omic)

#filter(dat, date==intro_date)

########## Test_prop set up
#-------EB: test_prop now has dates and this ensures it starts at the right date 
###########SK age-structured data (????), we'll use adjusted AB test_prop this week 

tp_approx <- tp_approx_fit(mytest=mytest_AB, dat= dat_omic, forecasts_days=forecasts_days, howlow = 0.1, 
                           slope = 0.18,  midpoint=35, intro_date= intro_date, stop_date=stop_date)
#should decline later compared to AB 

plot_fit <- tp_approx[1]
plot_fit 

# test_prop for fit and projection 
test_prop_proj_SK <- data.frame(tp_approx[2])


test_prop <- test_prop_proj_SK$tp[1:length(dat_omic$day)]
#subset for fitting alone (length of data)
#or set forecasts_days to 0 in tp_approx_fit() above 

modifytp = function(startlower=1, endfraction=1, test_prop=test_prop_proj$tp[1:length(dat_omic$day)]) {
  tt = min(which(test_prop < startlower))
  myline = seq(1, endfraction, length.out = length(test_prop)-tt)
  myfrac = c(rep(1, tt), myline)
  return(test_prop*myfrac)
}

### EB____ we should modify test_prop within modifytp()

mytp = modifytp(startlower = 1, endfraction = 1, test_prop)
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


########## Parameters set up 
#set values to generate initial conditions with make_init()

N=1.174e6 
N_pop=N
vaxlevel_in = 0.72 # portion of the pop vaccinated at start time () on Nov 30, 72% double vaxed 76% single dose
port_wane_in = 0.086 # portion boosted at start time 8.6% on Nov 30
past_infection_in = 0.15  #increased this from 0.1 to 0.18 # total in R at start time (336K so far, with 50% asc rate)
incres_in = 42*2#866#(made up factor) # resident strain (delta) incidence at start (433, reported cases) 
incmut_in = 2 # new (omicron) inc at stat (presumably very low)
simu_size = 1e5 # number of times to resample the negative binom (for ribbons)
forecasts_days =30 # how long to forecast for 
times = 1:nrow(dat_omic)




#declaring  parameters 
eff_date <-   ymd("2022-01-15")  # intervention date # restrictions were updated on 24 December 
intv_date <-  ymd("2022-02-10")
parameters <-         c(sigma=1/3, # incubation period (days) 
                        gamma=1/(5), #recovery rate 
                        nu =0.007, #vax rate: 0.7% per day 
                        mu=1/(82*365), # 1/life expectancy 
                        w1= 1/(0.5*365),# waning rate from R to S 
                        w2= 1/(0.5*365), # waning rate from Rv to V 
                        w3= 1/(0.5*365),# waning rate Rw to W 
                        ve=1, # I think this should be 1. it is not really efficacy  
                        beta_r=0.5, #transmission rate 
                        beta_m=1.05, #transmission rate 
                        epsilon_r = (1-0.8), # % this should be 1-ve 
                        epsilon_m = 1-0.3, # % escape capacity 
                        b= 0.006, # booster rate
                        beff = 0.7, # booster efficacy
                        wf=0.1, # protection for newly recovered
                        N=1.174e6,
                        stngcy= 0.4, #(*%(reduction)) strength of intervention (reduction in beta's)
                        relx_level = 0,
                        rlx_t = as.numeric(intv_date - intro_date),
                        eff_t = as.numeric(eff_date - intro_date),
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
pen.fit_SK <- optim(fn=func_penloglik,  par=guess, lower=c(0,0,0,0.001,0), upper = c(Inf,1,Inf,Inf,1), 
                    method = "L-BFGS-B", 
                    parameters = parameters, test_prop=test_prop, 
                    pen.size=pen.size, penalties = penalties, dat_omic=dat_omic, hessian=T)
pen.fit_SK
####################
func_loglik(pen.fit_SK$par, test_prop, dat_omic,parameters) 
parameters[names(guess)] <- pen.fit_SK$par


gg = simple_prev_plot(parameters, numdays = 190, mode = "both"); gg  # cc's simple prevalence plot 





########## Check fit: make prediction and projection with estimated parameters

#NOTE
#EB----------- use fitted test_prop instead of extendtp() for PHAC forecats 

mytp = mytp = modifytp(startlower = 1, endfraction = 1, 
                       test_prop = test_prop_proj_SK$tp[1:(length(dat_omic$day)+forecasts_days)])
test_prop <- mytp


################################# re-sampling CI

library(MASS)
times <- 1:(nrow(dat_omic) + forecasts_days)
# Sample from asymptotic distribution of MLEs, 100 times
resampled <- mvrnorm(n = 100, mu = pen.fit_SK$par, Sigma = (1/nrow(dat_omic))*solve(pen.fit_SK$hessian)) # don't need -ve because we MINimised NEGloglike

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


bound_mlesample_SK <-apply(incidence_resampled,2,  function(x){raply(simu_size/100,rnbinom(n=length(x),
                                                                                           mu=x, size=1/parameters[["theta"]]))}, simplify=FALSE)
bound_mlesample_SK <- do.call(rbind, bound_mlesample_SK)

project_dat_SK  =  as.data.frame(aaply(bound_mlesample_SK 
                                       ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
  mutate(date=seq.Date(ymd(intro_date),
                       ymd(intro_date)-1+length(times), 1))


#add dat to data for plotting 

dat_reported <- dat_omic  %>% mutate(date=seq.Date(ymd(intro_date),
                                                   ymd(intro_date)-1+length(dat_omic$day), 1))

ggplot() + geom_line(data=project_dat_SK,aes(x=date,y=`50%`), col="green",size=1.5,alpha=0.4) +
  geom_ribbon(data=project_dat_SK,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='darkgreen',alpha=0.1, size = 1.5)+
  geom_point(data=dat_reported,aes(x=date, y=value),color='grey48', alpha=0.8, size = 1.5)

# ---- end check fit ---- 


mu1 <- mean(dat_omic$value)
mu2 <- mu1/last(test_prop)

theta_2 = 1 / ( 1/mu1 - 1/mu2 + 1/parameters[["theta"]])  #m1 / end of test_prop

parameters[["theta"]] <- theta_2 

times <- 1:(nrow(dat_omic) + forecasts_days)
out_SK_rel <- as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                                         parms=parameters)) 

incidence_SK_rel =  parameters[["p"]]*parameters[["sigma"]]*(out_SK_rel$Er + out_SK_rel$Erv + out_SK_rel$Erw +
                                                               out_SK_rel$Em + out_SK_rel$Emv +
                                                               out_SK_rel$Emw)

uncert_bound_SK_rel = raply(simu_size,rnbinom(n=length(incidence_SK_rel),
                                              mu=incidence_SK_rel,
                                              size=1/parameters[["theta"]]))
project_dat_SK_rel =  as.data.frame(aaply(uncert_bound_SK_rel 
                                          ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
  mutate(date=seq.Date(ymd(intro_date),
                       ymd(intro_date)-1+length(times), 1))

plot(project_dat_SK_rel$`50%`)


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
lines(incidence_SK_rel)
lines(incidence_var$X1)
lines(incidence_var$X2)



uncert_bound_SK_rel_1 = raply(simu_size,rnbinom(n=length(incidence_var$X1),
                                                mu=incidence_var$X1,
                                                size=1/parameters[["theta"]]))

project_dat_SK_rel_1 =  as.data.frame(aaply(uncert_bound_SK_rel_1 
                                            ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
  mutate(date=seq.Date(ymd(intro_date),
                       ymd(intro_date)-1+length(times), 1))



uncert_bound_SK_rel_2 = raply(simu_size,rnbinom(n=length(incidence_var$X2),
                                                mu=incidence_var$X2,
                                                size=1/parameters[["theta"]]))

project_dat_SK_rel_2 =  as.data.frame(aaply(uncert_bound_SK_rel_2 
                                            ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
  mutate(date=seq.Date(ymd(intro_date),
                       ymd(intro_date)-1+length(times), 1))





##################


saveRDS(project_dat_SK, file.path("data/SK_test_constraints.rds"))
project_dat_SK =readRDS(file.path("data/SK_test_constraints.rds"))

saveRDS(project_dat_SK_rel, file.path("data/SK_no_constraints.rds"))
project_dat_SK_rel =readRDS(file.path("data/SK_no_constraints.rds"))

saveRDS(project_dat_SK_rel_1, file.path("data/SK_40percent_inc.rds"))
project_dat_SK_rel_1 = readRDS(file.path("data/SK_40percent_inc.rds"))

saveRDS(project_dat_SK_rel_2, file.path("data/SK_80percent_inc.rds"))
project_dat_SK_rel_2 = readRDS(file.path("data/SK_80percent_inc.rds"))
# "#", 
#make figures 

cols <- c("Current, TC" = "darkgreen",
          "NTC"="orange", "40% increase"="#0072B2", "80% increase"= "#CC79A7")

gg_SK <- ggplot() + geom_line(data=project_dat_SK,aes(x=date,y=`50%`, colour = "Current, TC"),size=1.2,alpha=0.4) +
  geom_ribbon(data=project_dat_SK,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='darkgreen',alpha=0.1)+
  geom_line(data=project_dat_SK_rel,aes(x=date,y=`50%`, color="NTC"),size=1.2,alpha=0.4) +
  geom_ribbon(data=project_dat_SK_rel,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='orange',alpha=0.1)+
  geom_line(data=project_dat_SK_rel_1,aes(x=date,y=`50%`, color="Without testing constraints"),size=1.2,alpha=0.4) +
  geom_ribbon(data=project_dat_SK_rel_1,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='yellow',alpha=0.1)+
  geom_point(data=dat_reported,aes(x=date, y=value),color='grey48', alpha=0.8) + 
  
  geom_line(data=project_dat_SK_rel_1,aes(x=date,y=`50%`, color="40% increase"),size=1.2,alpha=0.4) +
  geom_ribbon(data=project_dat_SK_rel_1,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill="purple",alpha=0.1)+
  
  geom_line(data=project_dat_SK_rel_2,aes(x=date,y=`50%`, color="80% increase"),size=1.2,alpha=0.4) +
  geom_ribbon(data=project_dat_SK_rel_2,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill="#D55E00",alpha=0.1)+
  
  #geom_line(aes(y=typical),color='blue') +
  labs(y="Reported cases",x="Date") + ylim(c(0,40000)) + 
  scale_x_date(date_breaks = "8 days", date_labels = "%b-%d-%y") +theme_light() +
  scale_color_manual(values = cols) +  theme(axis.text=element_text(size=15),
                                             plot.title = element_text(size=15, face="bold"),
                                             axis.text.x = element_text(angle = 45, hjust = 1),
                                             legend.position = "bottom", legend.title = element_text(size=15),
                                             legend.text = element_text(size=15),
                                             axis.title=element_text(size=15,face="bold")) +
  
  labs(color = " ",title="SK") +  geom_vline(xintercept=eff_date, linetype="dashed", 
                                             color = "grey", size=1)


gg_SK

ggsave(file="figs/SK_proj.png", gg_SK, width = 14, height = 8)
saveRDS(gg_SK, file.path("figs/SK-fig.rds"))


