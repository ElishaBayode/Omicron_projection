
####### Doing 'ba2' like projections for fall  - we project in March '22 under various scenarios and then shift the x-axis...
# This script follows on directly from pipps_simulation.r


# Adjust dates for projection
old_intro_date  = intro_date # Keep track of 'day 0'
intro_date <-   ymd("2022-03-10")
stop_date <- intro_date + 300 # 300 days projection
times = as.numeric(intro_date - old_intro_date):(as.numeric(stop_date - old_intro_date)) # Keep the same time count going




rem_parameters  <- parameters # ...to keep 'parameters' safe
# Set the desired characteristics of the new mutant. You can include any of the named elements of rem_parameters here -------------
# Choose which scenario you want to run:

# 1. Worst case scenario
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.5,
                        "epsilon_m" = 0.9, # Previous mutant epsilon was 0.7
                        "c_m" = rem_parameters["c_m"]*1,
                        "c_mr" = rem_parameters["c_mr"]*1.3,
                        "c_rm" = rem_parameters["c_rm"]*1,
                        "w_m" =  rem_parameters["w_m"]*1) 
IHR_factor <- 2 # multiplier for IHR (see below)

# 2. Intermediate case scenario
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.3,
                        "epsilon_m" = 0.8, # Previous mutant epsilon was 0.7
                        "c_m" = rem_parameters["c_m"]*1,
                        "c_mr" = rem_parameters["c_mr"]*1.3,
                        "c_rm" = rem_parameters["c_rm"]*1,
                        "w_m" =  rem_parameters["w_m"]*1) 
IHR_factor <- 2 # multiplier for IHR (see below)

# 3. Best case scenario
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.0,
                        "epsilon_m" = 0.8, # Previous mutant epsilon was 0.7
                        "c_m" = rem_parameters["c_m"]*1,
                        "c_mr" = rem_parameters["c_mr"]*1,
                        "c_rm" = rem_parameters["c_rm"]*1,
                        "w_m" =  rem_parameters["w_m"]*1) 
                        # May want to consider: "beff" =  rem_parameters["beff"]*0.7) 
IHR_factor <- 1 # multiplier for IHR (see below)




# Swap resident and mutant, then set up new mutant -------------
# This assumes that the new mutant 'arrives' with mut_prop% of current cases
new_model <- swap_strains(out_old = out_samp, params_old = rem_parameters, 
                          params_newmutant = params_newmutant, mut_prop = 0.3)
proj_parameters <- new_model$newm_parameters

# Make projections
proj_out <- as.data.frame(deSolve::ode(y=new_model$init_newm, time=times,func= sveirs,
                                       parms=proj_parameters)) 


# Plot of the projected cases -------------
proj_out <- proj_out %>% mutate(Total=
                                  proj_parameters[["sigma"]]*(proj_out$Er + proj_out$Erv + proj_out$Erw + 
                                                                proj_out$Em + proj_out$Emv + proj_out$Emw), 
                                "BA4/5"=
                                  proj_parameters[["sigma"]]*(proj_out$Er + proj_out$Erv + proj_out$Erw), 
                                "New variant X"=
                                  proj_parameters[["sigma"]]*(proj_out$Em + proj_out$Emv + proj_out$Emw)) %>% 
  mutate(date=seq.Date(ymd("2022-10-01"),ymd("2022-10-01")+300, 1)) 

pivot_longer(proj_out, c(Total,"BA4/5", "New variant X"), names_to = "Strain", values_to = "count") %>%
  ggplot(aes(x=date, y=count, colour=Strain)) + geom_line() + ylab("Incident cases") + xlab("Date") + theme_minimal()
  


# Plot of the projected  hospitalizations -------------
IHR = 0.01 * (5/16)*IHR_factor # see slack w nicola

proj_out <- proj_out %>% mutate(prev =Ir + Irv + Irw + Im + Imv +Imw) %>% mutate(hosp = Total*IHR*(1e5/N))  # note per 100K hosp
  
ggplot(proj_out) + geom_line(aes( x=date, y=hosp), col="blue") + theme_minimal() + ylab("Hospitalizations, per 100k population") + xlab("Date")





