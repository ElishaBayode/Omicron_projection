#source relevant scripts
source("analysis-new/functions.R") #functions
source("analysis-new/BC_data.R") # data
source("analysis-new/likelihood_func.R") # likelihood function
source("analysis-new/function-new.R")
########## Data setup  ######
#import data 
dat = readRDS("data/BC-dat.rds")

#resize data to include Omicron wave only---stopping in March to
#enable sensible comparison between model output and seroprevalence data

forecasts_days <- 1 
intro_date <-  ymd("2021-11-30")
stop_date <- ymd("2022-06-25")#last(dat$date)# #last date of data   

dat <- dat %>% filter(date >= intro_date &  date <= stop_date)
dat_omic <- dat
dat_omic <- filter(dat_omic, date >= intro_date) %>% dplyr::select(c("day", "value", "date"))
dat_omic$day <- 1:nrow(dat_omic)



#check spline 
# ggplot(pivot_longer(test_BCspline$upred,
#        cols = c(1,3),values_to = "number",names_to = "type"), 
#        aes(x=Reported_Date, y=number,color=type))+geom_point()

#test prop goes down to 0.33 during the Omicron wave 
x=intro_date+1:nrow(dat_omic)
tptest = 1 -0.67/(1+exp(-0.2*as.numeric((x - ymd("2021-12-20")))))
plot(x,tptest)
test_prop= tptest


########## Parameters setup 
#set values to generate initial conditions with make_init()
N=5.07e6
N_pop=N
vaxlevel_in = 0.88 # portion of the pop vaccinated at start time 
port_wane_in = 0.04 # portion boosted at start tie 
past_infection_in = 0.14  #increased this from 0.1 to 0.18 # total in R at start time 
#New-----changed to past_infection_in to 14% (by end of Nov 2021), seroprev was 9% in Sept/Oct 

incres_in = 389*1.78 #( 389, reported cases on Nov 30) # resident strain (delta) incidence at start 
incmut_in = 30# new (omicron) inc at stat (#first cases of Omicron reported on 30 Nov)
simu_size = 1e5 # number of times to resample the negative binom (for ribbons)
#forecasts_days =30 # how long to forecast for 
times = 1:nrow(dat_omic)

#declaring  parameters 
eff_date <-   ymd("2021-12-31")  # intervention date 
intv_date <-  ymd("2022-03-05") #  increase due to BA.2 (not needed here)
fur_intv_date <- ymd("2022-07-04") #increase due to reopening (not needed here)


parameters <-         c(sigma=1/3, # incubation period (days) 
                        gamma=1/5, #recovery rate 
                        nu =0.007, #vax rate: 0.7% per day 
                        mu=1/(82*365), # 1/life expectancy 
                        w_m= 1/(0.33*365),# waning rate from R to S 
                        w_r= 1/(0.33*365),
                        w_b= 1/(0.33*365), # waning rate from Rv to V 
                        ve=1, # I think this should be 1. it is not really efficacy  
                        beta_r=0.6, #transmission rate 
                        beta_m=1.12, #transmission rate 
                        c_m = 0.005,  #(1-%protection) protection from mutant variants when individuals  just recovered  from it (made-up (for now)) 
                        c_r = 0.005, #protection from resident variants when individuals  just recovered  from it (made-up (for now)) 
                        c_mr = 0.2, #cross immunity of resident  from mutant 
                        c_rm = 0.001, #cross immunity of mutant  from resident 
                        epsilon_r = (1-0.8), # % this should be 1-ve 
                        epsilon_m = 1-0.3, # % 1-ve omicron 
                        b= 0.018,#0.018, # booster rate
                        beff = 0.88, # booster efficacy
                        wf=0.01, # protection for newly recovered
                        N=5.07e6,
                        stngcy= 0.4, #(*%(reduction)) strength of intervention (reduction in beta's)
                        eff_t = as.numeric(eff_date - intro_date),
                        relx_level = 0.45,
                        fur_relx_level = 0,
                        rlx_t = as.numeric(intv_date - intro_date),
                        fur_rlx_t = as.numeric(fur_intv_date - intro_date),
                        p = 0.2, #negative binomial mean (from pre-Omicron seroprevalence estimates)
                        theta = 0.1 #negative binomial dispersion
                        
)
init <- make_init()   #generate initial states





########## Run the model and compare the model to the data at least for IC checking
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

# I think this just shows the initial conditions and parameters and what they do
ggplot(data =inctest, aes(x=date, y = inc_reported))+geom_line() +
  geom_line(aes(x=date, y= inc_res), color = "blue") +
  geom_line(aes(x=date, y= inc_mut), color = "red") +
  geom_point(data = dat, aes(x=date, y=cases), alpha=0.5) +
  ylim(c(0,55000))# + xlim(c(ymd("2021-11-30"), ymd("2022-02-28")))


# Fitting 

known_prop <- 0.5 # known the proportion of resident and new strain
date_known_prop <- "2021-12-12" # what time was the proportion known? 
# 2. growth advantage of mutant strain
known_growth <- 0.2 # daily rise of omicron on dec 12 
period_known_growth <- c("2021-12-05", "2021-12-15")
penalties <- list(known_prop = known_prop, date_known_prop = date_known_prop, 
                  known_growth = known_growth, period_known_growth = period_known_growth)

# Determine weight of penalty. 
# Qu: how strong should penalty be on scale of 0-1? 0 = no penalty. 1 = relatively as impactful as the likelihood
pen.size <- 0.1

# Guess starting parameters and fit the model 

guess <- c( beta_m=1, beta_r=0.6, theta=0.1, stngcy=0.4,p=0.2) 

#guess <- c( beta_m=1, stngcy=0.4,beta_r=0.6, theta=0.1,beff=0.8) 



pen.fit_BC <- optim(fn=func_penloglik,  par=guess, lower=c(0,0,0.01, 0,0), 
                    upper = c(Inf,Inf,1,1,1), 
                    method = "L-BFGS-B", 
                    parameters = parameters, test_prop=test_prop, 
                    pen.size=pen.size, penalties = penalties, dat_omic=dat_omic, hessian=T)
pen.fit_BC



# check the fit 
func_loglik(pen.fit_BC$par, test_prop, dat_omic,parameters) 
parameters[names(guess)] <- pen.fit_BC$par


#check total cases and compare to seroprevalence data 


out_samp <- as.data.frame(deSolve::ode(y=init,time=times,func=sveirs,
                                       parms=parameters)) 

reportable = parameters[["p"]]*parameters[["sigma"]]*(out_samp$Er + out_samp$Erv + out_samp$Erw +
                                                        out_samp$Em + out_samp$Emv +
                                                        out_samp$Emw)
incidence = test_prop*reportable #(length differs)



incid = data.frame(date = intro_date+out_samp$time, 
                   inc = incidence, reportable = reportable)
tmp =  merge(dat, incid, by="date") %>% dplyr::select(date, cases, inc, reportable) %>% 
  pivot_longer(c(2,3,4), names_to = "type", values_to = "number")
ggplot(tmp, aes(x=date, y=number, color=type))+geom_line()        

###############sanity checks against seroprevalence data ############

#from data
#1:5.3 (18% ascertainment from Sept/Oct 2021 to March 2022)
#9% inf prevalence in Sept/Oct 2021, and 43% in March 2022 
#seroprevalence of vaccinated or infected 96% (from serprev)


true_incidence =  reportable/parameters[["p"]]
tot_true <- sum(true_incidence)
tot_reported <- sum(dat_omic$value)

#from model
tot_true/N #(27% infection from Nov 30, 2021 to March 30, 2022)
tot_reported/tot_true #(9% ascertainment from Nov 30, 2021 to March 30, 2022)

tot_inf_vax <- (out_samp$V+ out_samp$Erv+ out_samp$Emv+out_samp$Irv + out_samp$Imv+  out_samp$Rv + out_samp$R + 
                  out_samp$W+out_samp$Erw+out_samp$Emw+out_samp$Irw + out_samp$Imw + out_samp$Rw)

last(tot_inf_vax/N) #(95.5%, consistent)



################################# re-sampling CI


times <- 1:(nrow(dat_omic) + 0)
# Sample from asymptotic distribution of MLEs, 100 times
resampled <- mvrnorm(n = 100, mu = pen.fit_BC$par,Sigma = (1/nrow(dat_omic))*solve(pen.fit_BC$hessian)) # don't need -ve because we MINimised NEGloglike


#sample a new transmission rate from a distribution, then use that to calculate the selecytion coefficient 



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

project_dat_BC  =  as.data.frame(aaply(bound_mlesample_BC,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
  mutate(date=seq.Date(ymd(intro_date),ymd(intro_date)-1+length(times), 1))


# #last date of data   


prevelence <- data.frame(prev=(out_samp$Ir + out_samp$Im + out_samp$Irv + out_samp$Imv + out_samp$Irw + out_samp$Imw))
prevelence <- prevelence %>%  mutate(date=seq.Date(ymd(intro_date), ymd(intro_date)+length(times), 1))

ggplot() + geom_line(data=prevelence,aes(x=date,y=prev), col="green",size=1.5,alpha=0.4)

dat_reported <- dat_omic  %>% mutate(date=seq.Date(ymd(intro_date),ymd(intro_date)-1+length(dat_omic$day), 1))



ggplot() + geom_line(data=project_dat_BC,aes(x=date,y=`50%`), col="green",size=1.5,alpha=0.4) +
  geom_point(data=dat_reported,aes(x=date, y=value),color='grey48', alpha=0.8, size = 1.5) +
  geom_ribbon(data=project_dat_BC,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='darkgreen',alpha=0.1, size = 1.5) +
  labs(y="Reported cases",x="Date") + ylim(c(0,10000)) + #30000
  scale_x_date(date_breaks = "15 days", date_labels = "%b-%d-%y") +theme_light() +
  scale_color_manual(values = cols) +  theme(axis.text=element_text(size=12),
                                             plot.title = element_text(size=15, face="bold"),
                                             axis.text.x = element_text(angle = 45, hjust = 1),
                                             legend.position = "bottom", legend.title = element_text(size=15),
                                             legend.text = element_text(size=12),
                                             axis.title=element_text(size=12,face="bold")) 


#saveRDS(project_dat_BC, file.path("data/BC_fit_omicron.rds"))






