#source relevant scripts
source("analysis-new/functions.R") #functions
source("analysis-new/BC_data.R") # data
source("analysis-new/likelihood_func.R") # likelihood function
source("analysis-new/function-new.R")
########## Data setup  ######
#import data 
dat = readRDS("data/BC-dat.rds")


# citf data 
citf_bc = data.frame(date = c(ymd("2021-12-22"),ymd("2022-01-16"), ymd("2022-02-14"), ymd("2022-03-16"), 
                              ymd("2022-04-15"), ymd("2022-05-15")), 
                     cumulperc = c(5.83, 11.82, 25.52, 32.09, 35.01, 45.43)) 
#ba2 wave in BC is mid-march to late may. see notes in pipps_projection
sum(filter(dat_full, date< ymd("2021-12-22"))$cases)/N # number reported 
ascprop2021 = (sum(filter(dat_full, date< ymd("2021-12-22"))$cases)/N )/ (0.08) # reported cases / approx serology
#resize data to include Omicron wave only---stopping in March to
#enable sensible comparison between model output and seroprevalence data


forecasts_days <- 1 
intro_date <-  ymd("2021-11-30")
stop_date <- ymd("2022-03-10")#last(dat$date)# #last date of data   

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
tptest = 1 -0.9/(1+exp(-0.2*as.numeric((x - ymd("2021-12-20")))))
# tptest = tptest - (0.97-0.67)/(1+exp(-0.2*as.numeric((x - ymd("2022-03-01")))))
plot(x,tptest)
test_prop= tptest
plot(x,test_prop)


########## Parameters setup 
#set values to generate initial conditions with make_init()
N=5.07e6
# Neff = 2.7e6 # this made no difference 
N_pop= N
vaxlevel_in = 0.88 # portion of the pop vaccinated at start time 
port_wane_in = 0.07 # portion boosted at start time 
past_infection_in = 0.14  #increased this from 0.1 to 0.18 # total in R at start time 
#New-----changed to past_infection_in to 14% (by end of Nov 2021), seroprev was 9% in Sept/Oct 

incres_in = 389 #( 389, then was 389*1.78 reported cases on Nov 30) # resident strain (delta) incidence at start 
incmut_in = 30# new (omicron) inc at stat (#first cases of Omicron reported on 30 Nov)
simu_size = 1e5 # number of times to resample the negative binom (for ribbons)
#forecasts_days =30 # how long to forecast for 
times = 1:nrow(dat_omic)

#declaring  parameters 
eff_date <-   ymd("2021-12-31")  # intervention date 
intv_date <-  ymd("2024-12-08") #  increase due to BA.2 - turned off
fur_intv_date <- ymd("2024-12-09") #increase due to reopening  - turned off

########################################
# notes 
########################################
# first experiment - basically don't wane , don't reinfect w same variant; 
# ba2: don't have immune escape, but allow increase in beta
# method: set waning to 1 year (still 1/2 year for boosters). set cs to 0. 
# result: ba2 wave looks good, needs 55% increase in transmissibility, 
# still maybe slightly the wrong shape. hosp admissions look great, 
# census looks great but lags behind admissions (14 day lag from cases to census, 6 to admissions)
# but both look great 

# second experiment: reduce the infectiousness of vax ppl
# new ode func, sveirs.vred, with two new pars, vred and vredb (reduction in beta for vax and boosted)
# bigger ba1 wave than we think we had; fiddling a little (w vredb and ve at 0.3 instead of 0.45) 
# that results in a slow ba2 wave under teh same assns (good 1 to 2 cross immunity) 
# would require an even higher increase in transmissibility for ba2. 
# with higher ve in ba1, would likely be better 
# hosps not as good in this reduced beta model - basically yes, this makes the overall 
# transmissibility even lower, and so beta (baseline) has to be really high
# to make a quick ba2 wave. 
# HOWEVER - if vredb (booster effect) goes UP in the BA2 wave, it can look good, even with 
# waning immunity as in our default. hosps end up too high, cases too, but the shape is good.
# i'd say not as good as the longer-immunity model . see save.image("experiment2-vredb.Rdata")


parameters <-         c(sigma=1/1, # incubation period (days) 
                        gamma=1/4, #recovery rate 
                        nu =0.007, #vax rate: 0.7% per day 
                        mu=1/(82*365), # 1/life expectancy 
                        w_m= 1/(0.62*365),# waning rate from R to S 
                        w_r= 1/(0.62*365),
                        w_b= 1/(0.5*365), # waning rate from Rv to V 
                        ve=1, # I think this should be 1. it is not really efficacy  
                        beta_r=0.7, #transmission rate 
                        beta_m=1, #transmission rate 
                        c_m = 0.05,  # KEEP AT 0.05 for now . (1-protection) protection from mutant variants when individuals  just recovered  from it (made-up (for now)) 
                        c_r = 0.05, # KEEP AT 0.05 for now1- protection from resident variants when individuals  just recovered  from it (made-up (for now)) 
                        c_mr = 0.05, # CAN BE HIGHER than 0.05 1 -cross immunity of resident  from mutant 
                        c_rm = 0.05, #KEEP AT 0.05 for now1-cross immunity of mutant  from resident 
                        epsilon_r = (1-0.8), # % this should be 1-ve against delta - fine at 0.8
                        epsilon_m = (1-0.3), # % 1-ve omicron . ve omicron is at most 0.3 #DISCUSSED WITH JS
                        b= 0.015,#0.018, # booster rate DISCUSSED
                        beffr = 0.75, # booster efficacy, resident strain
                        beffm = 0.75, # booster efficacy, mutant strain. 
                        N=N_pop,
                        stngcy= 0.4, #(*%(reduction)) strength of intervention (reduction in beta's)
                        eff_t = as.numeric(eff_date - intro_date),
                        relx_level = 0.65,
                        fur_relx_level = 0,
                        rlx_t = as.numeric(intv_date - intro_date),
                        fur_rlx_t = as.numeric(fur_intv_date - intro_date),
                        p = ascprop2021, # ascertainment fraction from pre-Omicron seroprevalence estimates
                        theta = 0.1, #negative binomial dispersion
                        vred=0.8, 
                        vredb=0.5
                        
)
init <- make_init(N=N_pop)   #generate initial states

# NOTE HERE - REMOVE FOR ANOTHER EXPT 
sveirs = sveirs.vred



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

# check booster
vv = get_vax(outtest)
ggplot(vv, aes(x=time, y=boosted/N))+geom_line() 


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
known_growth <- 0.18 # daily rise of omicron on dec 12 
period_known_growth <- c("2021-12-05", "2021-12-15")
penalties <- list(known_prop = known_prop, date_known_prop = date_known_prop, 
                  known_growth = known_growth, period_known_growth = period_known_growth)

# Determine weight of penalty. 
# Qu: how strong should penalty be on scale of 0-1? 0 = no penalty. 1 = relatively as impactful as the likelihood
pen.size <- 0.1

# Guess starting parameters and fit the model 

guess <- c( beta_m=1, beta_r=0.6, theta=0.1, stngcy=0.4)  # NOTE removed fitting p 

#guess <- c( beta_m=1, stngcy=0.4,beta_r=0.6, theta=0.1,beff=0.8) 



exp.fit_BC <- optim(fn=func_penloglik,  par=guess, lower=c(0,0,0.01, 0,0.01), 
                    upper = c(Inf,Inf,1,1,1), 
                    method = "L-BFGS-B", 
                    parameters = parameters, test_prop=test_prop, 
                    pen.size=pen.size, penalties = penalties, dat_omic=dat_omic, hessian=T)
exp.fit_BC



# check the fit 
func_loglik(exp.fit_BC$par, test_prop, dat_omic,parameters) 
parameters[names(guess)] <- exp.fit_BC$par


#check total cases and compare to seroprevalence data 


out_samp <- as.data.frame(deSolve::ode(y=init,time=times,func=sveirs,
                                       parms=parameters)) 
get_growth_rate(out_samp, startoffset = 2, duration = 10)

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
tot_true/N # compare to Danuta's estimate of 34% of BC in the BA1 wave. 
tot_reported/tot_true # approx ascertainment 

tot_inf_vax <- (out_samp$V+ out_samp$Erv+ out_samp$Emv+out_samp$Irv + out_samp$Imv+   out_samp$Rrv+out_samp$Rmv + 
                  out_samp$W+out_samp$Erw+out_samp$Emw+out_samp$Irw + out_samp$Imw + out_samp$Rrw + out_samp$Rmw)

last(tot_inf_vax/N) #(95.5%, consistent)
# check delta and omicron growth rate 
get_growth_rate(out_samp, startoffset = 2, duration = 10)


# check vaccination / booster level 
vv = get_vax(out_samp)
ggplot(vv, aes(x=time, y=boosted/N))+geom_line()

# check two strain dynamics: did delta actually decline ? 
incdf= get_total_incidence(out_samp, parameters)
incdf$date = incdf$time+ymd("2021-12-01")-1 
ggplot(incdf, aes(x=date, y=inc_res/parameters["p"]))+geom_line() + 
  geom_line(data = incdf, aes(x=date, y = inc_mut/parameters["p"]), inherit.aes = F)

# look at our incident cases along with Tara Moriarty's estimates 
taram = readr::read_csv("data/Canadian COVID data_Infections_Line chart.csv")
taram = taram %>% mutate(chardate = date) %>% mutate(date = dmy(chardate))

# tara says broadly consistent with serology, and indeed: 
ba1tara = filter(taram, date > ymd("2021-12-01") & date < ymd("2022-03-30"))
sum(ba1tara$BC)/N # 38% -- high. Danuta had only 34% of BC in BA1 and citf has only 27%. 

ggplot(ba1tara, aes(x=date, y = BC))+geom_point(alpha=0.5)+
  geom_line(data = incdf, inherit.aes = F, aes(x=date, y = inc_mut/parameters["p"]), 
            color="blue") + 
  scale_x_date(date_breaks = "months", date_labels = "%b-%d") +theme_light()
# i think what we have now is a good enough interpolation between what we know from 



########################
# Quick check with hosps
########################
source("analysis-new/hosp-data.R")
hospdat <- get_hosp_data(intro_date, stop_date)
IHR <- get_IHR() # note - the old IHR from Nicola's fit was 0.0035 I think

# prediction
predict_hosps <- data.frame(inc=reportable/parameters[["p"]],
   date=seq.Date(ymd(intro_date), ymd(intro_date)+length(reportable)-1, 1)) %>%
  mutate(hosp=lag(inc*IHR,6))

# plot
ggplot(hospdat, aes(x=week_of, y=new/7))+#weekly to daily
  geom_point(size=2.5)+
  geom_point(col="grey")+
  geom_line(data=predict_hosps, aes(x=date, y=hosp), col="darkblue", size=1.5)+
  labs(x="Date", y="Predicted Hospital Admissions")+
  theme_light()

# save.image(file = "simulationscript_out.Rdata")

########################
# Nicola's ORIGINAL quick check with hosps
########################

# **use admission data once we get it
hosp_data <- get_can_covid_tracker_data("bc") %>%
  mutate(date=as.Date(date)) %>%
  dplyr::select("date", "total_hospitalizations") %>%
  filter(date <= stop_date & date >= intro_date) %>%
  rename(hosp_census = total_hospitalizations) %>%
  mutate(hosp_admit = as.numeric(hosp_census)/8) # approx. based on av stay in hosp

l <- 14# lag
ls_fit_hosp <- function(par, data){# super basic..
  predict <- par[["IHR"]]*lag(incid$reportable, l)/parameters[["p"]]
  predict <- predict[1:length(data)]
  return(sum((predict-data)^2, na.rm=TRUE))
}

guess <- c(IHR=0.001)
fit_BC_hosp <- optim(fn=ls_fit_hosp,  par=guess, lower=c(0),
                     upper = c(1), method = "L-BFGS-B", data = hosp_data$hosp_admit)

# get result
IHR <- fit_BC_hosp$par[["IHR"]]
predict <- data.frame(n=IHR*lag(incid$reportable, l)/parameters[["p"]],
                      date = incid$date)

ggplot(hosp_data, aes(x=date, y=hosp_admit))+
  geom_point(col="grey")+
  geom_line(data=predict, aes(x=date, y=n), col="blue")+
  labs(x="", y="Predicted Hospital Admissions")+
  theme_light()

########################################################################
##### now ba2 ####
########################################################################

dat_full = readRDS("data/BC-dat.rds")
tot_ba2_estimate = 1.26*1.15*0.13 # percent of BC infected in ba2

forecasts_days <- 1 
old_intro_date  = intro_date # Keep track of 'day 0'
intro_date <-   ymd("2022-02-27") # cc made this a little earlier. 
stop_date <- last(dat_full$date)

dat_rem <- dat_full %>% filter(date >= intro_date &  date <= stop_date)
# align 'day' to correct date, in order to continue counting on same scale as _simulation.R
dat_rem$day <- as.numeric(intro_date - old_intro_date):(as.numeric(stop_date - old_intro_date))

# eff_date <-   ymd("2021-12-31")

#dat_rem <- filter(dat_full, date >= intro_date) %>% select(c("day", "value", "date"))
dat_full$day <- 1:nrow(dat_full)

plot(dat_rem$date, dat_rem$value)

#now using Jessica's function to switch variants 

rem_parameters  <- parameters

times = dat_rem$day

# Set the desired characteristics of the new mutant.
# change booster rate etc as needed 
# You can include any of the named elements of rem_parameters here
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.55,#1.11
                        "gamma"=1/4,
                        "sigma"=1, 
                        "epsilon_m" = (1-0.4), # was (1-0.45) in pipps-simn
                        "c_m" = rem_parameters["c_m"]*1,#BA.2's protection against itself higher than BA.1's?
                        "c_mr" = rem_parameters["c_mr"]*1, # lowering this (wo other changes) slows it down. 
                        "c_rm" = rem_parameters["c_rm"]*1,
                        "w_m" =  rem_parameters["w_m"]*1,
                        "b"=1/(0.4*365), # lower booster rate 
                        "w_b" = 1/(0.5*365), # set booster waning to around booster rate 
                        "beffr" = 0.75, 
                        "beffm"=0.75,
                        "vredb"=0.75
)



# Swap resident and mutant, then set up new mutant. 
# This assumes that the new mutant 'arrives' with mut_prop% of current cases
out_samp$date = out_samp$time + old_intro_date -1 
out_samp = dplyr::filter(out_samp, date < intro_date)
new_model <- swap_strains(out_old = out_samp, params_old = rem_parameters, 
                          params_newmutant = params_newmutant, 
                          mut_prop = 0.6, res_to_s_prop =  0.3)
init_proj <- new_model$init_newm
proj_parameters <- new_model$newm_parameters

# Make projections
#forecasts_days <- nrow(dat_rem)
times <- dat_rem$day # Keep the same time count going
proj_out <- as.data.frame(deSolve::ode(y=init_proj, time=times,func= sveirs,
                                       parms=proj_parameters)) 
#check growth rate 
get_growth_rate(output= proj_out, startoffset = 20, duration = 7)

# check booster
#vv = get_vax(proj_out)
# ggplot(vv, aes(x=time, y=boosted/N))+geom_line()


# Simple plot of the projection
proj_out <- proj_out %>% mutate(Total=last(test_prop)*proj_parameters[["p"]]*
                                  proj_parameters[["sigma"]]*(proj_out$Er + proj_out$Erv + proj_out$Erw + 
                                                                proj_out$Em + proj_out$Emv + proj_out$Emw), 
                                Resident=last(test_prop)*proj_parameters[["p"]]*
                                  proj_parameters[["sigma"]]*(proj_out$Er + proj_out$Erv + proj_out$Erw), 
                                Mutant=last(test_prop)*proj_parameters[["p"]]*
                                  proj_parameters[["sigma"]]*(proj_out$Em + proj_out$Emv + proj_out$Emw)) %>% 
  mutate(date=seq.Date(ymd(intro_date),ymd(intro_date )-1+length(times), 1)) 

#pivot_longer(proj_out, c(Total, Resident, Mutant), names_to = "Strain", values_to = "count") %>%
ggplot(proj_out) + geom_line(aes( x=date, y=Resident), col="blue") + 
  geom_line(aes( x=date, y=Mutant), col="red") +  geom_line(aes( x=date, y=Total), col="green") + 
  geom_point(aes(x=dat_rem$date, y=dat_rem$value)) # this is really only useful for peak time

# and what do we know about ba2 in BC? 
# ww: peak approx 1/2 the ba1 peak, late april, trough in 
# mid-June
# serology with fudge - about 13-18% of BC got BA2. 
# hosp and ww: it is approx 1/2 the peak height of BA1 

# add plot showing BA1 so we can compare
incba1 = dplyr::filter(incdf, date < intro_date+2) %>% dplyr::select(date, inc_tot) 
incba2 = proj_out[3:nrow(proj_out),] %>% dplyr::select(date, Total) %>% rename(inc_tot = Total) %>% 
  mutate(inc_tot = inc_tot / last(test_prop))
ggplot(rbind(incba1, incba2), aes(x=date, y=inc_tot))+geom_line()+ 
  scale_x_date(date_breaks = "months", date_labels = "%b-%d")


# 
infected =get_total_infection(proj_out, from_date = ymd("2022-03-16"), 
                              to_date = ymd("2022-05-15"), 
                              parameters=proj_parameters)
infected/N # 




# -- Hospitalizations -------------
source("analysis-new/hosp-data.R")
hospdat <- get_hosp_data(intro_date, stop_date)
# IHR <- get_IHR()*1.1 # account for reporting change...

proj_out <- proj_out %>% 
  mutate(incid=proj_parameters[["sigma"]]*
           (Er + Erv + Erw + Em + Emv + Emw), 
         prev = Ir + Irv + Irw + Im + Imv + Imw) %>% 
  mutate(hosp = lag(incid,6)*IHR) %>%  
  mutate(date=seq.Date(ymd(intro_date),ymd(intro_date )-1+length(times), 1)) 

ggplot(hospdat, aes(x=week_of, y=new/7))+ #weekly to daily
  geom_point(size=2.5)+
  geom_point(col="grey")+
  geom_line(data=proj_out, aes(x=date, y=hosp), col="darkblue", size=1.5)+
  labs(x="Date", y="Predicted Hospital Admissions")+
  theme_light()


# --- check against census numbers (/8) ---- # 
hosp_data <- get_can_covid_tracker_data("bc") %>%
  mutate(date=as.Date(date)) %>%
  dplyr::select("date", "total_hospitalizations") %>%
  filter(date <= stop_date & date >= intro_date) %>%
  rename(hosp_census = total_hospitalizations) %>%
  mutate(hosp_admit = as.numeric(hosp_census)/8) # approx. based on av stay in hosp


ggplot(hosp_data, aes(x=date, y=hosp_admit))+
  geom_point(col="grey")+
  geom_line(data=proj_out, aes(x=date, y=lag(hosp,8)), col="blue")+
  labs(x="", y="Model compared to adjusted public census")+
  theme_light()




