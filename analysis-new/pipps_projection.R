#load("simulationscript_out.Rdata") - if you want to read in the saved output of simulation.R

#add dat to data for plotting 
#here I match model output to the remaining data point (i.e. from March 31, 2022 onward)
dat_full = readRDS("data/BC-dat.rds")

# approximately f1*f2*(13%) new infections in taht time period
# where f1, f2 account for  citf under-rep and the fact that there were some reinfections
# take f1 from the ba1 - ba2 difference. Danuta: 34% of BC in the BA1 wave
# CITF: 4.68 (sept-oct) to end ba1 (march) 32.09 -> 27%. 
# Difference is f1 = 1.26. 
# f2 depends on efficacy of first infection, but probably a lot of people getting 
# ba2 didn't get ba1. say f2 is about 1.15. 
# then the total infections in the ba2 wave should be 
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

#find how much intervention and relaxation have impacted transmission 
with(as.list( rem_parameters), {
  c <- (1 - stngcy/(1+ exp(-1.25*(times-eff_t)))) 
  rlx <- (1 + relx_level/(1+ exp(-1.25*(times-rlx_t)))) # relaxation 
  further_rlx <- (1 + fur_relx_level/(1+ exp(-1.25*(times-fur_rlx_t))))
  infectionfactor <- c*rlx*further_rlx
  plot(times,infectionfactor)})

#increase due to reopening 
#rem_parameters["beta_r"] <- rem_parameters["beta_r"]*(1.5) #rem_parameters["beta_r"]*0.4*(1.5) - JS: I think this is already being done by rlx_t?
#rem_parameters["beta_m"] <- rem_parameters["beta_m"]*(1.5) #rem_parameters["beta_m"]*0.4*(1.5)
#rem_parameters["epsilon_r"] <- (1-0.15) 
#rem_parameters[["stngcy"]] <- 0.35
#rem_parameters[["p"]] <- 0.5
#rem_parameters[["beff"]] <- 0.95
#rem_parameters[["wf"]] <- 0.01
#rem_parameters[["b"]] <- 0.018*1.2
#initialm data matching 
x = last(out_samp)[,2:ncol(out_samp)]
pie(as.numeric(x),  
    labels = colnames(out_samp)[2:ncol(out_samp)], radius = 1.5)

# Set the desired characteristics of the new mutant.
# change booster rate etc as needed 
# You can include any of the named elements of rem_parameters here
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.30,#1.11
                        "gamma"=1/4,
                        "sigma"=1, 
                        "epsilon_m" = (1-0.4), # was (1-0.3) where 0.3 is ve 
                        "c_m" = rem_parameters["c_m"]*1,#BA.2's protection against itself higher than BA.1's?
                        "c_mr" = rem_parameters["c_mr"]*2, # lowering this (wo other changes) slows it down. 
                        "c_rm" = rem_parameters["c_rm"]*1,
                        "w_m" =  rem_parameters["w_m"]*1,
                        "b"=1/(0.4*365), # lower booster rate 
                        "w_b" = 1/(0.5*365), # set booster waning to around booster rate 
                        "beffr" = 0.75, 
                        "beffm"=0.75
                        
                      #  "eff_t"= 137,
                      #  "stngcy"=0.21
)



# Swap resident and mutant, then set up new mutant. 
# This assumes that the new mutant 'arrives' with mut_prop% of current cases
out_samp$date = out_samp$time + old_intro_date -1 
out_samp = dplyr::filter(out_samp, date < intro_date)
new_model <- swap_strains(out_old = out_samp, params_old = rem_parameters, 
                          params_newmutant = params_newmutant, 
                          mut_prop = 0.6, res_to_s_prop =  0.6)
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
source("~/Omicron_projection/analysis-new/hosp-data.R")
hospdat <- get_hosp_data(intro_date, stop_date)
IHR <- get_IHR()*1.1 # account for reporting change...

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


# -- Split cases across HAs ------------
source("analysis-new/pipps_geographical.R")
# 'which_wave_match' tells this function whether to make a 'delta-like' wave, a 'ba.1-like wave' and so on 
#                                                               - you can currently provide any wave 1:7 (7 = ba.2)
project_HAs(total_out = proj_out, which_wave_match = 6)



save.image(file = "projectionscript_out.Rdata")

