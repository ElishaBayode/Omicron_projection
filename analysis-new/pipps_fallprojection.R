
####### Fall scenarios - we project in summer under various scenarios and then shift the x-axis to a hypothetical fall date

# This script follows on directly from pipps_projection.r
load("projectionscript_out.Rdata")


# Adjust dates for projection
old_intro_date  = intro_date # Keep track of 'day 0'
intro_date <-   ymd("2022-07-01")
stop_date <- intro_date + 300 # 300 days projection
times = as.numeric(intro_date - old_intro_date):(as.numeric(stop_date - old_intro_date)) # Keep the same time count going



rem_parameters  <- proj_parameters # ...to keep parameters safe
# Set the desired characteristics of the new mutant. You can include any of the named elements of rem_parameters 
# Choose which scenario you want to run:

# Which date to plot projections out until
maxdate = ymd("2023-02-01")


#### Non-omicron scenarios ######################################
#### 1. Worst case scenario:  ---------------------------
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
projfilename = "NO_worstcase-projection.csv"
initfilename = "NO_worstcase-init.csv"
parfilename = "NO_worstcase-par.csv"

# Swap resident and mutant, then set up new mutant 
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

plot.worst = ggplot(filter(fproj_out, date < maxdate), aes(x=date, y=fproj_parameters["p"]*Total/50))+
  geom_line(color="blue") +  ylab("Incidence (symptomatic infection per 100K)") + xlab("")
plot.worst

NO_worstdf = fproj_out %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters["p"]*Total,
         Scenario = "Worst case")

readr::write_csv(fproj_out, file = projfilename)
readr::write_csv(as.data.frame(new_model$init_newm), file = initfilename)
readr::write_csv(as.data.frame(fproj_parameters), file = parfilename)




##########################################################
#### 2. Intermediate case scenario   ---------------------------
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.2,
                        "epsilon_m" = (1-0.3), # Previous mutant epsilon was 0.6
                        "c_m" = rem_parameters["c_m"]*1,
                        "c_mr" = rem_parameters["c_mr"]*1.1,
                        "c_rm" = rem_parameters["c_rm"]*1,
                        "w_m" =  rem_parameters["w_m"]*1,
                        "beffr" = 0.70,
                        "beffm"=0.70)
IHR_factor <- 1 # multiplier for IHR (see below)
res_swap = 0.55
projfilename = "NO_mediumcase-projection.csv"
initfilename = "NO_mediumcase-init.csv"
parfilename = "NO_mediumcase-par.csv"


# Swap resident and mutant, then set up new mutant 
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

plot.medium = ggplot(filter(fproj_out, date < maxdate), aes(x=date, y=fproj_parameters["p"]*Total/50))+
  geom_line(color="blue") +  ylab("Incidence (symptomatic infection per 100K)") + xlab("")
plot.medium

NO_mediumdf = fproj_out %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters["p"]*Total,
         Scenario = "Medium case" )

readr::write_csv(fproj_out, file = projfilename)
readr::write_csv(as.data.frame(new_model$init_newm), file = initfilename)
readr::write_csv(as.data.frame(fproj_parameters), file = parfilename)





#### Omicron scenarios ######################################
#### 1. Worst case scenario ---------------------------
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.1, # this one diff to inter
                        "epsilon_m" = (1-0.3), 
                        "c_m" = rem_parameters["c_m"]*1.0,
                        "c_mr" = rem_parameters["c_mr"]*1.0, 
                        "c_rm" = rem_parameters["c_rm"]*1.0,
                        "w_m" =  rem_parameters["w_m"]*1.0,
                        "beffr" = 0.7,
                        "beffm"=0.7)

IHR_factor <- 1 # multiplier for IHR (see below)
res_swap = 0.3
projfilename = "O_worstcase-projection.csv"
initfilename = "O_worstcase-init.csv"
parfilename = "O_worstcase-par.csv"

# Swap resident and mutant, then set up new mutant -------------
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
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

pivot_longer(fproj_out, c(Total,"Established type", "New variant X"), names_to = "Strain", values_to = "count") %>%
  ggplot(aes(x=date, y=count, colour=Strain)) + geom_line() + ylab("Incident cases") + xlab("Date") + theme_minimal()

plot.worst = ggplot(filter(fproj_out, date < maxdate), aes(x=date, y=fproj_parameters["p"]*Total/50))+
  geom_line(color="blue") +  ylab("Incidence (symptomatic infection per 100K)") + xlab("")
plot.worst

O_worstdf = fproj_out %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters["p"]*Total,
         Scenario = "Worst case" )

readr::write_csv(fproj_out, file = projfilename)
readr::write_csv(as.data.frame(new_model$init_newm), file = initfilename)
readr::write_csv(as.data.frame(fproj_parameters), file = parfilename)

#### 2. Intermediate case scenario (Increased beta & epsilon, gives ba.2 like wave) --------
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.1,
                        "epsilon_m" = (1-0.3), # Previous mutant epsilon was 0.6
                        "c_m" = rem_parameters["c_m"]*1.0,
                        "c_mr" = rem_parameters["c_mr"]*1.0,
                        "c_rm" = rem_parameters["c_rm"]*1.0,
                        "w_m" =  rem_parameters["w_m"]*1.0,
                        "beffr" = 0.75,
                        "beffm"=0.75)

IHR_factor <- 1 # multiplier for IHR (see below)
res_swap = 0.3
projfilename = "O_mediumcase-projection.csv"
initfilename = "O_mediumcase-init.csv"
parfilename = "O_mediumcase-par.csv"

# Swap resident and mutant, then set up new mutant -------------
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
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

pivot_longer(fproj_out, c(Total,"Established type", "New variant X"), names_to = "Strain", values_to = "count") %>%
  ggplot(aes(x=date, y=count, colour=Strain)) + geom_line() + ylab("Incident cases") + xlab("Date") + theme_minimal()

plot.medium = ggplot(filter(fproj_out, date < maxdate), aes(x=date, y=fproj_parameters["p"]*Total/50))+
  geom_line(color="blue") +  ylab("Incidence (symptomatic infection per 100K)") + xlab("")
plot.medium

O_mediumdf = fproj_out %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters["p"]*Total,
         Scenario = "Medium case" )

readr::write_csv(fproj_out, file = projfilename)
readr::write_csv(as.data.frame(new_model$init_newm), file = initfilename)
readr::write_csv(as.data.frame(fproj_parameters), file = parfilename)

#### 3. Best case scenario (continued decreased trend, from ba.2-like parameters) --------
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.0,
                        "epsilon_m" = (1-0.4), # Previous mutant epsilon was 0.6
                        "c_m" = rem_parameters["c_m"]*1.0,
                        "c_mr" = rem_parameters["c_mr"]*1.0,
                        "c_rm" = rem_parameters["c_rm"]*1.0,
                        "w_m" =  rem_parameters["w_m"]*1.0,
                        "beffr" = 0.75,
                        "beffm"=0.75)

IHR_factor <- 1 # multiplier for IHR (see below)
res_swap = 0.3
projfilename = "O_bestcase-projection.csv"
initfilename = "O_bestcase-init.csv"
parfilename = "O_bestcase-par.csv"

# Swap resident and mutant, then set up new mutant -------------
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
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

pivot_longer(fproj_out, c(Total,"Established type", "New variant X"), names_to = "Strain", values_to = "count") %>%
  ggplot(aes(x=date, y=count, colour=Strain)) + geom_line() + ylab("Incident cases") + xlab("Date") + theme_minimal()

plot.best = ggplot(filter(fproj_out, date < maxdate), aes(x=date, y=fproj_parameters["p"]*Total/50))+
  geom_line(color="blue") +  ylab("Incidence (symptomatic infection per 100K)") + xlab("")
plot.best

O_bestdf = fproj_out %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters["p"]*Total,
         Scenario = "Best case" )

readr::write_csv(fproj_out, file = projfilename)
readr::write_csv(as.data.frame(new_model$init_newm), file = initfilename)
readr::write_csv(as.data.frame(fproj_parameters), file = parfilename)



#### Make figures with (i) non-omicron and (ii) omicron scenarios, 
# for infections and hosps
# 


# (i) Non-omicron scenarios
# Added grey box so we can highlight that likely we would DO something about the increased cases
alldf  = rbind(NO_mediumdf, NO_worstdf)
glimpse(alldf)
ggplot(filter(alldf, date<maxdate), aes(x=date, y = sympinfect/50, fill=Scenario))+
  ylab("Incidence (symptomatic infection per 100K)") + 
  annotate("rect", xmin = as.Date("2022-11-01"), xmax = as.Date("2023-02-01"), ymin = 0, ymax = Inf, alpha=0.2, fill="grey") +
  geom_ribbon(aes(x=date, ymin = 0.85*sympinfect/50, ymax=1.25*sympinfect/50,
                  fill=Scenario),  alpha=0.5) + theme_minimal() +
  theme(legend.position = "bottom") + xlab("") +  
  scale_x_continuous(breaks= c(as.Date("2022-10-01"),as.Date("2022-11-01"),as.Date("2022-12-01"),as.Date("2023-01-01"),as.Date("2023-02-01")), labels = c("Month 1","Month 2","Month 3","Month 4","Month 5"))

# now make one with hospitalizations
NO_worstdf = NO_worstdf %>% mutate(hosp = Total*IHR_BA2*2, census = lag(Total,8)*IHR_BA2*2*9) # 2x increase
NO_mediumdf = NO_mediumdf %>% mutate(hosp = Total*IHR_BA2, census = lag(Total,8)*IHR_BA2*9) # no increase
alldf  = rbind(NO_mediumdf, NO_worstdf)

ggplot(filter(alldf, date<maxdate), aes(x=date, y = hosp, fill=Scenario))+
  ylab("COVID- all reported admissions") +
  annotate("rect", xmin = as.Date("2022-11-01"), xmax = as.Date("2023-02-01"), ymin = 0, ymax = Inf, alpha=0.2, fill="grey") +
  geom_ribbon(aes(x=date, ymin = 0.75*hosp, ymax=1.25*hosp,
                  fill=Scenario),  alpha=0.5) + theme_minimal() +
  theme(legend.position = "bottom") +  
  #geom_hline(yintercept=73, linetype="dashed") + geom_text(aes(as.Date("2022-10-03"),73,label = "Approx. BA.2 peak", vjust = -1))+ 
  #geom_hline(yintercept=110, linetype="dashed") + geom_text(aes(as.Date("2022-10-03"),110,label = "Approx. BA.1 peak", vjust = -1))+
  xlab("") +  
  scale_x_continuous(breaks= c(as.Date("2022-10-01"),as.Date("2022-11-01"),as.Date("2022-12-01"),as.Date("2023-01-01"),as.Date("2023-02-01")), labels = c("Month 1","Month 2","Month 3","Month 4","Month 5"))


ggplot(filter(alldf, date<maxdate), aes(x=date, y = census, fill=Scenario))+
  ylab("COVID- census hospitalizations") +
  annotate("rect", xmin = as.Date("2022-11-01"), xmax = as.Date("2023-02-01"), ymin = 0, ymax = Inf, alpha=0.2, fill="grey") +
  geom_ribbon(aes(x=date, ymin = 0.75*census, ymax=1.25*census,
                  fill=Scenario),  alpha=0.5) + theme_minimal() +
  theme(legend.position = "bottom") + 
  #geom_hline(yintercept=590, linetype="dashed") + geom_text(aes(as.Date("2022-10-03"),590,label = "Approx. BA.2 peak", vjust = -1))+ 
  #geom_hline(yintercept=1000, linetype="dashed") + geom_text(aes(as.Date("2022-10-03"),1000,label = "Approx. BA.1 peak", vjust = -1)) +
  xlab("")+  
  scale_x_continuous(breaks= c(as.Date("2022-10-01"),as.Date("2022-11-01"),as.Date("2022-12-01"),as.Date("2023-01-01"),as.Date("2023-02-01")), labels = c("Month 1","Month 2","Month 3","Month 4","Month 5"))




# (ii) omicron scenarios
alldf  = rbind(O_bestdf, O_mediumdf, O_worstdf)
glimpse(alldf)
ggplot(filter(alldf, date<maxdate), aes(x=date, y = sympinfect/50, fill=Scenario))+
   ylab("Incidence (symptomatic infection per 100K)") +
  geom_ribbon(aes(x=date, ymin = 0.85*sympinfect/50, ymax=1.25*sympinfect/50,
                  fill=Scenario),  alpha=0.5) + theme_minimal() +
  theme(legend.position = "bottom") + xlab("")

# now make one with hospitalizations
O_worstdf = O_worstdf %>% mutate(hosp = Total*IHR_BA2, census = lag(Total,8)*IHR_BA2*9) # no increase
O_mediumdf = O_mediumdf %>% mutate(hosp = Total*IHR_BA2, census = lag(Total,8)*IHR_BA2*9) # no increase
O_bestdf = O_bestdf %>% mutate(hosp = Total*IHR_BA2, census = lag(Total,8)*IHR_BA2*9) # no increase
alldf  = rbind(O_bestdf, O_mediumdf, O_worstdf)

ggplot(filter(alldf, date<maxdate), aes(x=date, y = hosp, fill=Scenario))+
  ylab("COVID- all reported admissions") +
  geom_ribbon(aes(x=date, ymin = 0.75*hosp, ymax=1.25*hosp,
                  fill=Scenario),  alpha=0.5) + theme_minimal() +
  theme(legend.position = "bottom") + 
  #geom_hline(yintercept=73, linetype="dashed") + geom_text(aes(as.Date("2022-09-05"),73,label = "Approx. BA.2 peak", vjust = -1))+ 
  #geom_hline(yintercept=110, linetype="dashed") + geom_text(aes(as.Date("2022-09-05"),110,label = "Approx. BA.1 peak", vjust = -1))+
  xlab("")
  
ggplot(filter(alldf, date<maxdate), aes(x=date, y = census, fill=Scenario))+
  ylab("COVID- census hospitalizations") +
  geom_ribbon(aes(x=date, ymin = 0.75*census, ymax=1.25*census,
                  fill=Scenario),  alpha=0.5) + theme_minimal() +
  theme(legend.position = "bottom") +  
 # geom_hline(yintercept=590, linetype="dashed") + geom_text(aes(as.Date("2022-09-05"),590,label = "Approx. BA.2 peak", vjust = -1))+ 
 # geom_hline(yintercept=1000, linetype="dashed") + geom_text(aes(as.Date("2022-09-05"),1000,label = "Approx. BA.1 peak", vjust = -1))+
  xlab("")




# Split projected cases across age groups and HAs ------------
source("analysis-new/functions_splittingwaves.R")
# 'which_wave_match' tells these functions whether to make a 'delta-like' wave, a 'ba.1-like wave' and so on
#                                                               - you can currently provide any wave 1:7 (7 = ba.2)
NO_worst.HAs <- project_HAs(total_out = NO_worstdf[NO_worstdf$date<=maxdate,], 
                         which_wave_match = 5, facets = TRUE)
NO_medium.HAs <-project_HAs(total_out = NO_mediumdf[NO_mediumdf$date<=maxdate,],
                         which_wave_match = 5, facets = TRUE)
O_worst.HAs <- project_HAs(total_out = O_worstdf[O_worstdf$date<=maxdate,], 
                         which_wave_match = 5, facets = TRUE)
O_medium.HAs <-project_HAs(total_out = O_mediumdf[O_mediumdf$date<=maxdate,],
                         which_wave_match = 5, facets = TRUE)
O_best.HAs <-project_HAs(total_out = O_bestdf[O_bestdf$date<=maxdate,], 
                       which_wave_match = 5, facets = TRUE)
NO_worst.HAs$plot + ggtitle("Worst case")
NO_medium.HAs$plot + ggtitle("Medium case")
O_worst.HAs$plot + ggtitle("Worst case")
O_medium.HAs$plot + ggtitle("Medium case")
O_best.HAs$plot + ggtitle("Best case")
readr::write_csv(NO_worst.HAs$df, file = "NO_worstcase_byHA.csv")
readr::write_csv(NO_medium.HAs$df, file = "NO_mediumcase_byHA.csv")
readr::write_csv(O_worst.HAs$df, file = "O_worstcase_byHA.csv")
readr::write_csv(O_medium.HAs$df, file = "O_mediumcase_byHA.csv")
readr::write_csv(O_best.HAs$df, file = "O_bestcase_byHA.csv")

NO_worst.ages <- project_ages(total_out = NO_worstdf[NO_worstdf$date<=maxdate,], 
                            which_wave_match = 5, facets = TRUE)
NO_medium.ages <-project_ages(total_out = NO_mediumdf[NO_mediumdf$date<=maxdate,],
                            which_wave_match = 5, facets = TRUE)
O_worst.ages <- project_ages(total_out = O_worstdf[O_worstdf$date<=maxdate,], 
                           which_wave_match = 5, facets = TRUE)
O_medium.ages <-project_ages(total_out = O_mediumdf[O_mediumdf$date<=maxdate,],
                           which_wave_match = 5, facets = TRUE)
O_best.ages <-project_ages(total_out = O_bestdf[O_bestdf$date<=maxdate,], 
                         which_wave_match = 5, facets = TRUE)

NO_worst.ages$plot + ggtitle("Worst case")
NO_medium.ages$plot + ggtitle("Medium case")
O_worst.ages$plot + ggtitle("Worst case")
O_medium.ages$plot + ggtitle("Medium case")
O_best.ages$plot + ggtitle("Best case")
readr::write_csv(NO_worst.ages$df, file = "NO_worstcase_byage.csv")
readr::write_csv(NO_medium.ages$df, file = "NO_mediumcase_byage.csv")
readr::write_csv(O_worst.ages$df, file = "O_worstcase_byage.csv")
readr::write_csv(O_medium.ages$df, file = "O_mediumcase_byage.csv")
readr::write_csv(O_best.ages$df, file = "O_bestcase_byage.csv")






# Plot of the projected  hospitalizations -------------
#IHRx = IHR_BA2*IHR_factor # see slack w nicola
#fproj_out <- fproj_out %>% mutate(prev =Ir + Irv + Irw + Im + Imv +Imw) %>% mutate(hosp = Total*IHRx*(1e5/N))  # note per 100K hosp
#ggplot(fproj_out) + geom_line(aes( x=date, y=hosp), col="blue") + theme_minimal() + ylab("Hospitalizations, per 100k population") + xlab("Date")







