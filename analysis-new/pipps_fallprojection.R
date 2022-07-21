
####### Doing 'ba2' like projections for fall  - we project in March '22 under various scenarios and then shift the x-axis...

# This script follows on directly from pipps_projection.r
load("projectionscript_out.Rdata")



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
params_newmutant = list("beta_r" = rem_parameters["beta_r"]*0.9, # compensate for artificial boost
                        "beta_m" = rem_parameters["beta_m"]*1.3,
                        "gamma"=1/4,
                        "sigma"=1,
                        #    "epsilon_m" = (1-0.4), for ref this is what i did for ba2
                        "epsilon_m" = (1-0.3), #
                        "c_m" = rem_parameters["c_m"]*1,
                        "c_mr" = rem_parameters["c_mr"]*1.2,
                        "c_rm" = rem_parameters["c_rm"]*1,
                        "w_m" =  rem_parameters["w_m"]*1,
                        "beffr" = 0.7,
                        "beffm"=0.7)

IHR_factor <- 2 # multiplier for IHR (see below)
res_swap = 0.7
projfilename = "worstcase-projection.csv"
projfilename = "worstcase-projection.csv"
projfilename = "worstcase-projection.csv"


# 2. Intermediate case scenario
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.2,
                        "epsilon_m" = (1-0.3), # Previous mutant epsilon was 0.7
                        "c_m" = rem_parameters["c_m"]*1,
                        "c_mr" = rem_parameters["c_mr"]*1.1,
                        "c_rm" = rem_parameters["c_rm"]*1,
                        "w_m" =  rem_parameters["w_m"]*1,
                        "beffr" = 0.70,
                        "beffm"=0.70)
IHR_factor <- 1 # multiplier for IHR (see below)
res_swap = 0.55
projfilename = "mediumcase-projection.csv"
projfilename = "mediumcase-projection.csv"
projfilename = "mediumcase-projection.csv"


# 3. Best case scenario - limited immune evasion etc
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.1,
                        "epsilon_m" = (1-0.3), # Previous mutant epsilon was 0.7
                        "c_m" = rem_parameters["c_m"]*1.1,
                        "c_mr" = rem_parameters["c_mr"]*1.1,
                        "c_rm" = rem_parameters["c_rm"]*1,
                        "w_m" =  rem_parameters["w_m"]*1,
                        "beffr" = 0.75,
                        "beffm"=0.75)

IHR_factor <- 1 # multiplier for IHR (see below)
res_swap = 0.3
projfilename = "bestcase-projection.csv"
projfilename = "bestcase-projection.csv"
projfilename = "bestcase-projection.csv"




# Swap resident and mutant, then set up new mutant -------------
# This assumes that the new mutant 'arrives' with mut_prop% of current cases
new_model <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                          params_newmutant = params_newmutant, mut_prop = 0.1,
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
                                "Established type"=
                                  fproj_parameters[["sigma"]]*(fproj_out$Er + fproj_out$Erv + fproj_out$Erw),
                                "New variant X"=
                                  fproj_parameters[["sigma"]]*(fproj_out$Em + fproj_out$Emv + fproj_out$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-10-01"),ymd("2022-10-01")+300, 1))

pivot_longer(fproj_out, c(Total,"Established type", "New variant X"), names_to = "Strain", values_to = "count") %>%
  ggplot(aes(x=date, y=count, colour=Strain)) + geom_line() + ylab("Incident cases") + xlab("Date") + theme_minimal()

maxdate = ymd("2023-02-01")

# change name according to which scenario you run
plot.best = ggplot(filter(fproj_out, date < maxdate), aes(x=date, y=fproj_parameters["p"]*Total/50))+
  geom_line(color="blue") +  ylab("Incidence (symptomatic infection per 100K)") + xlab("")

plot.best

# CHANGE THIS LINE DEPENDING ON WHAT SCENARIO YOU RUN !!
bestdf = fproj_out %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters["p"]*Total,
         Scenario = "best case" )

readr::write_csv(fproj_out, file = projfilename)
readr::write_csv(as.data.frame(new_model$init_newm), file = initfilename)
readr::write_csv(as.data.frame(fproj_parameters), file = parfilename)


# add plot showing BA2 so we can compare
#incba2 = proj_out[3:nrow(proj_out),] %>% dplyr::select(date, incid) 
#incbaX  = fproj_out[3:nrow(fproj_out),] %>% dplyr::select(date, Total) %>% rename(incid = Total)
#ggplot(rbind(incba1, incba2, incbaX), aes(x=date, y=incid))+geom_line()+ 
#  scale_x_date(date_breaks = "months", date_labels = "%b-%d")





# having run those , make one figure with all of them , for infections and hosps

alldf  = rbind(bestdf, mediumdf, worstdf)
glimpse(alldf)
ggplot(filter(alldf, date<maxdate), aes(x=date, y = sympinfect/50, fill=Scenario))+
   ylab("Incidence (symptomatic infection per 100K)") +
  geom_ribbon(aes(x=date, ymin = 0.85*sympinfect/50, ymax=1.25*sympinfect/50,
                  fill=Scenario),  alpha=0.5)+
  theme(legend.position = "bottom") + xlab("")

# now make one with hospitalizations

worstdf = worstdf %>% mutate(hosp = Total*IHR*2) # 2x increase
mediumdf = mediumdf %>% mutate(hosp = Total*IHR) # no increase
bestdf = bestdf %>% mutate(hosp = Total*IHR) # no increase
alldf  = rbind(bestdf, mediumdf, worstdf)

ggplot(filter(alldf, date<maxdate), aes(x=date, y = hosp, fill=Scenario))+
  ylab("COVID admissions") +
  geom_ribbon(aes(x=date, ymin = 0.75*hosp, ymax=1.25*hosp,
                  fill=Scenario),  alpha=0.5)+
  theme(legend.position = "bottom") + xlab("")

# Split projected cases across HAs ------------
source("analysis-new/functions_splittingwaves.R")
# 'which_wave_match' tells this function whether to make a 'delta-like' wave, a 'ba.1-like wave' and so on
#                                                               - you can currently provide any wave 1:7 (7 = ba.2)
project_HAs(total_out = fproj_out, which_wave_match = 6)


# Plot of the projected  hospitalizations -------------
IHRx = IHR*IHR_factor # see slack w nicola

fproj_out <- fproj_out %>% mutate(prev =Ir + Irv + Irw + Im + Imv +Imw) %>% mutate(hosp = Total*IHRx*(1e5/N))  # note per 100K hosp

ggplot(fproj_out) + geom_line(aes( x=date, y=hosp), col="blue") + theme_minimal() + ylab("Hospitalizations, per 100k population") + xlab("Date")
