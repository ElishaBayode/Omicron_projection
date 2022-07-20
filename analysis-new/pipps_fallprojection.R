
####### Doing 'ba2' like projections for fall  - we project in March '22 under various scenarios and then shift the x-axis...

# This script follows on directly from either pipps_simulation.r or pipps_projection.r
#load("simulationscript_out.Rdata")
# load("projectionscript_out.Rdata")


# Adjust dates for projection
old_intro_date  = intro_date # Keep track of 'day 0'
intro_date <-   ymd("2022-07-01")
stop_date <- intro_date + 300 # 300 days projection
times = as.numeric(intro_date - old_intro_date):(as.numeric(stop_date - old_intro_date)) # Keep the same time count going



rem_parameters  <- proj_parameters # ...to keep parameters safe
# Set the desired characteristics of the new mutant. You can include any of the named elements of rem_parameters here -------------
# Choose which scenario you want to run:

# 1. Worst case scenario: cc modified to explore what makes what worse.
# insights -- only get a really high hosp peak if beta is higher. 
# but the other parameters, unsurprisingly, make the eventual endemic level (jan ++ ) really bad
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.4,
                        "gamma"=1/4,
                        "sigma"=3, 
                     #    "epsilon_m" = (1-0.4), for ref this is what i did for ba2
                        "epsilon_m" = (1-0.25), # Previous mutant epsilon was 0.7
                        "c_m" = rem_parameters["c_m"]*1.3,
                        "c_mr" = rem_parameters["c_mr"]*1.3,
                        "c_rm" = rem_parameters["c_rm"]*1,
                        "w_m" =  rem_parameters["w_m"]*1) 
IHR_factor <- 2 # multiplier for IHR (see below)
res_swap = 0.7

# 2. Intermediate case scenario
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.3,
                        "epsilon_m" = 0.8, # Previous mutant epsilon was 0.7
                        "c_m" = rem_parameters["c_m"]*1.3,
                        "c_mr" = rem_parameters["c_mr"]*1.3,
                        "c_rm" = rem_parameters["c_rm"]*1,
                        "w_m" =  rem_parameters["w_m"]*1) 
IHR_factor <- 2 # multiplier for IHR (see below)
res_swap = 0.3

# 3. Best case scenario
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.1,
                        "epsilon_m" = 0.8, # Previous mutant epsilon was 0.7
                        "c_m" = rem_parameters["c_m"]*1.3,
                        "c_mr" = rem_parameters["c_mr"]*1.3,
                        "c_rm" = rem_parameters["c_rm"]*1,
                        "w_m" =  rem_parameters["w_m"]*1) 
                        # May want to consider: "beff" =  rem_parameters["beff"]*0.7) 
IHR_factor <- 1 # multiplier for IHR (see below)
res_swap = 0.3




# Swap resident and mutant, then set up new mutant -------------
# This assumes that the new mutant 'arrives' with mut_prop% of current cases
new_model <- swap_strains(out_old = proj_out, params_old = rem_parameters, 
                          params_newmutant = params_newmutant, mut_prop = 0.5, 
                          res_to_s_prop =  res_swap)
fproj_parameters <- new_model$newm_parameters

# Make projections
fproj_out <- as.data.frame(deSolve::ode(y=new_model$init_newm, time=times,func= sveirs,
                                       parms=fproj_parameters)) 

get_growth_rate(fproj_out, startoffset = 5, duration = 10)

# Plot of the projected cases -------------
fproj_out <- fproj_out %>% mutate(Total=
                                  fproj_parameters[["sigma"]]*(fproj_out$Er + fproj_out$Erv + fproj_out$Erw + 
                                                                fproj_out$Em + fproj_out$Emv + fproj_out$Emw), 
                                "BA4/5"=
                                  fproj_parameters[["sigma"]]*(fproj_out$Er + fproj_out$Erv + fproj_out$Erw), 
                                "New variant X"=
                                  fproj_parameters[["sigma"]]*(fproj_out$Em + fproj_out$Emv + fproj_out$Emw)) %>% 
  mutate(date=seq.Date(ymd("2022-10-01"),ymd("2022-10-01")+300, 1)) 

pivot_longer(fproj_out, c(Total,"BA4/5", "New variant X"), names_to = "Strain", values_to = "count") %>%
  ggplot(aes(x=date, y=count, colour=Strain)) + geom_line() + ylab("Incident cases") + xlab("Date") + theme_minimal()
 
# Split projected cases across HAs ------------
source("analysis-new/pipps_geographical.R")
# 'which_wave_match' tells this function whether to make a 'delta-like' wave, a 'ba.1-like wave' and so on 
#                                                               - you can currently provide any wave 1:7 (7 = ba.2)
project_HAs(total_out = fproj_out, which_wave_match = 5) 


# Plot of the projected  hospitalizations -------------
IHRx = IHR*IHR_factor # see slack w nicola

fproj_out <- fproj_out %>% mutate(prev =Ir + Irv + Irw + Im + Imv +Imw) %>% mutate(hosp = Total*IHRx*(1e5/N))  # note per 100K hosp
  
ggplot(fproj_out) + geom_line(aes( x=date, y=hosp), col="blue") + theme_minimal() + ylab("Hospitalizations, per 100k population") + xlab("Date")





