#load("simulationscript_out.Rdata") - if you want to read in the saved output of simulation.R

#add dat to data for plotting 
#here I match model output to the remaining data point (i.e. from March 31, 2022 onward)
dat_full = readRDS("data/BC-dat.rds")



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


# Set the desired characteristics of the new mutant. You can include any of the named elements of rem_parameters here
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.55,#1.11
                        "epsilon_m" = 0.88,
                        "c_m" = rem_parameters["c_m"]*1,#BA.2's protection against itself higher than BA.1's?
                        "c_mr" = rem_parameters["c_mr"]*1, # lowering this (wo other changes) slows it down. 
                        "c_rm" = rem_parameters["c_rm"]*1,
                        "w_m" =  rem_parameters["w_m"]*1
                      #  "eff_t"= 137,
                      #  "stngcy"=0.21
)



# Swap resident and mutant, then set up new mutant. 
# This assumes that the new mutant 'arrives' with mut_prop% of current cases
new_model <- swap_strains(out_old = out_samp, params_old = rem_parameters, 
                          params_newmutant = params_newmutant, mut_prop = 0.35)
init_proj <- new_model$init_newm
proj_parameters <- new_model$newm_parameters

# Make projections
#forecasts_days <- nrow(dat_rem)
times <- dat_rem$day # Keep the same time count going
proj_out <- as.data.frame(deSolve::ode(y=init_proj, time=times,func= sveirs,
                                       parms=proj_parameters)) 
#check growth rate 
get_growth_rate(output= proj_out, startoffset = 20, duration = 7)


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
  geom_point(aes(x=dat_rem$date, y=dat_rem$value))



# -- Hospitalizations -------------
# IHR = 0.01 * (5/16) # see slack w nicola
IHR <- 0.00258 # from pipps_simulation
l <- 14 # from pipps_simulation. NOTE - if the age distribution changed, the lag would change too 
# seems reasonable for ba2 to have a higher IHR and lower lag, but of course not a negative lag... 

hosp_data <- get_can_covid_tracker_data("bc") %>%
  mutate(date=as.Date(date)) %>%
  dplyr::select("date", "total_hospitalizations") %>%
  filter(date %in% proj_out$date) %>%
  mutate(hosp_admit = as.numeric(total_hospitalizations)/8)

proj_out <- proj_out %>% 
  mutate(incid=proj_parameters[["sigma"]]*
           (Er + Erv + Erw + Em + Emv + Emw), 
         prev = Ir + Irv + Irw + Im + Imv + Imw) %>% 
  mutate(hosp = lag(incid,l)*IHR) %>%  
  mutate(date=seq.Date(ymd(intro_date),ymd(intro_date )-1+length(times), 1)) 

ggplot(proj_out) + 
  geom_line(aes(x=date, y=hosp), col="blue") +
  geom_point(data=hosp_data, aes(x=date, y=hosp_admit), col="grey")+
  labs(x="",y="Projected Hospital Admissions")


# -- Split cases across HAs ------------
source("analysis-new/pipps_geographical.R")
# 'which_wave_match' tells this function whether to make a 'delta-like' wave, a 'ba.1-like wave' and so on 
#                                                               - you can currently provide any wave 1:7 (7 = ba.2)
project_HAs(total_out = proj_out, which_wave_match = 5)





