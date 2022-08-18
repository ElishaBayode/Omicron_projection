
####### Sensitivity analysis for fall scenarios and ba1 fits


### Fall scenarios - vary those params that define the scenarios
load("projectionscript_out.Rdata")

# Adjust dates for projection
old_intro_date  = intro_date # Keep track of 'day 0'
intro_date <-   ymd("2022-07-01")
stop_date <- intro_date + 300 # 300 days projection
times = as.numeric(intro_date - old_intro_date):(as.numeric(stop_date - old_intro_date)) # Keep the same time count going

rem_parameters  <- proj_parameters # ...to keep parameters safe
# Which date to plot projections out until
maxdate = ymd("2023-02-01")


##### 1. Omicron - vary beta_m

par_to_vary <- "beta_m"
par_vals <- c(1.1, 1.05, 1.15, 1.1, 1.05, 1.15, 1.0, 0.95, 1.05) # worst, med, best (within = base, down, up)

# Worst case scenario - baseline ---------------------------
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.1, 
                        "epsilon_m" = (1-0.3), 
                        "c_m" = rem_parameters["c_m"]*1.0,
                        "c_mr" = rem_parameters["c_mr"]*1.0, 
                        "c_rm" = rem_parameters["c_rm"]*1.0,
                        "w_m" =  rem_parameters["w_m"]*1.0,
                        "beffr" = 0.7,
                        "beffm"=0.7)
params_newmutant[[par_to_vary]] = rem_parameters[par_to_vary]*par_vals[1]
IHR_factor <- 1 # multiplier for IHR (see below)
res_swap = 0.3

# Swap resident and mutant, then set up new mutant -------------
new_model_base <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                          params_newmutant = params_newmutant, mut_prop = 0.1,
                          res_to_s_prop =  res_swap)
fproj_parameters_base <- new_model_base$newm_parameters

params_newmutant[[par_to_vary]] = rem_parameters[par_to_vary]*par_vals[2]
new_model_down <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                               params_newmutant = params_newmutant, mut_prop = 0.1,
                               res_to_s_prop =  res_swap)
fproj_parameters_down <- new_model_down$newm_parameters

params_newmutant[[par_to_vary]] = rem_parameters[par_to_vary]*par_vals[3]
new_model_up <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                               params_newmutant = params_newmutant, mut_prop = 0.1,
                               res_to_s_prop =  res_swap)
fproj_parameters_up <- new_model_up$newm_parameters

# Make projections
fproj_out_base <- as.data.frame(deSolve::ode(y=new_model_base$init_newm, time=times,func= sveirs,
                                        parms=fproj_parameters_base))
fproj_out_down <- as.data.frame(deSolve::ode(y=new_model_down$init_newm, time=times,func= sveirs,
                                             parms=fproj_parameters_down))
fproj_out_up <- as.data.frame(deSolve::ode(y=new_model_up$init_newm, time=times,func= sveirs,
                                             parms=fproj_parameters_up))

fproj_out_base <- fproj_out_base %>% mutate(Total=
                                    fproj_parameters_base[["sigma"]]*(fproj_out_base$Er + fproj_out_base$Erv + fproj_out_base$Erw +
                                                                   fproj_out_base$Em + fproj_out_base$Emv + fproj_out_base$Emw),
                                  "Established type"=
                                    fproj_parameters_base[["sigma"]]*(fproj_out_base$Er + fproj_out_base$Erv + fproj_out_base$Erw),
                                  "New variant X"=
                                    fproj_parameters_base[["sigma"]]*(fproj_out_base$Em + fproj_out_base$Emv + fproj_out_base$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

fproj_out_down <- fproj_out_down %>% mutate(Total=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Er + fproj_out_down$Erv + fproj_out_down$Erw +
                                                                                  fproj_out_down$Em + fproj_out_down$Emv + fproj_out_down$Emw),
                                            "Established type"=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Er + fproj_out_down$Erv + fproj_out_down$Erw),
                                            "New variant X"=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Em + fproj_out_down$Emv + fproj_out_down$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

fproj_out_up <- fproj_out_up %>% mutate(Total=
                                              fproj_parameters_up[["sigma"]]*(fproj_out_up$Er + fproj_out_up$Erv + fproj_out_up$Erw +
                                                                                  fproj_out_up$Em + fproj_out_up$Emv + fproj_out_up$Emw),
                                            "Established type"=
                                              fproj_parameters_up[["sigma"]]*(fproj_out_up$Er + fproj_out_up$Erv + fproj_out_up$Erw),
                                            "New variant X"=
                                              fproj_parameters_up[["sigma"]]*(fproj_out_up$Em + fproj_out_up$Emv + fproj_out_up$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

O_worstdf = fproj_out_base %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters_base["p"]*Total,
         Scenario = "Pessimistic, beta=1.1 (baseline)" )
O_worstdf <- rbind(O_worstdf, fproj_out_down %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters_down["p"]*Total,
         Scenario = "Pessimistic, beta=1.05" ))
O_worstdf <- rbind(O_worstdf, fproj_out_up %>% dplyr::select(date, Total) %>%
                     mutate(sympinfect = fproj_parameters_up["p"]*Total,
                            Scenario = "Pessimistic, beta=1.15" ))

#### 2. Intermediate case scenario (Increased beta & epsilon, gives ba.2 like wave) --------
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.1,
                        "epsilon_m" = (1-0.3), # Previous mutant epsilon was 0.6
                        "c_m" = rem_parameters["c_m"]*1.0,
                        "c_mr" = rem_parameters["c_mr"]*1.0,
                        "c_rm" = rem_parameters["c_rm"]*1.0,
                        "w_m" =  rem_parameters["w_m"]*1.0,
                        "beffr" = 0.75,
                        "beffm"=0.75)
params_newmutant[[par_to_vary]] = rem_parameters[par_to_vary]*par_vals[4]
IHR_factor <- 1 # multiplier for IHR (see below)
res_swap = 0.3


# Swap resident and mutant, then set up new mutant -------------
new_model_base <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                               params_newmutant = params_newmutant, mut_prop = 0.1,
                               res_to_s_prop =  res_swap)
fproj_parameters_base <- new_model_base$newm_parameters

params_newmutant[[par_to_vary]] = rem_parameters[par_to_vary]*par_vals[5]
new_model_down <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                               params_newmutant = params_newmutant, mut_prop = 0.1,
                               res_to_s_prop =  res_swap)
fproj_parameters_down <- new_model_down$newm_parameters

params_newmutant[[par_to_vary]] = rem_parameters[par_to_vary]*par_vals[6]
new_model_up <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                             params_newmutant = params_newmutant, mut_prop = 0.1,
                             res_to_s_prop =  res_swap)
fproj_parameters_up <- new_model_up$newm_parameters

# Make projections
fproj_out_base <- as.data.frame(deSolve::ode(y=new_model_base$init_newm, time=times,func= sveirs,
                                             parms=fproj_parameters_base))
fproj_out_down <- as.data.frame(deSolve::ode(y=new_model_down$init_newm, time=times,func= sveirs,
                                             parms=fproj_parameters_down))
fproj_out_up <- as.data.frame(deSolve::ode(y=new_model_up$init_newm, time=times,func= sveirs,
                                           parms=fproj_parameters_up))

fproj_out_base <- fproj_out_base %>% mutate(Total=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Er + fproj_out_base$Erv + fproj_out_base$Erw +
                                                                                  fproj_out_base$Em + fproj_out_base$Emv + fproj_out_base$Emw),
                                            "Established type"=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Er + fproj_out_base$Erv + fproj_out_base$Erw),
                                            "New variant X"=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Em + fproj_out_base$Emv + fproj_out_base$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

fproj_out_down <- fproj_out_down %>% mutate(Total=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Er + fproj_out_down$Erv + fproj_out_down$Erw +
                                                                                  fproj_out_down$Em + fproj_out_down$Emv + fproj_out_down$Emw),
                                            "Established type"=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Er + fproj_out_down$Erv + fproj_out_down$Erw),
                                            "New variant X"=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Em + fproj_out_down$Emv + fproj_out_down$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

fproj_out_up <- fproj_out_up %>% mutate(Total=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Er + fproj_out_up$Erv + fproj_out_up$Erw +
                                                                            fproj_out_up$Em + fproj_out_up$Emv + fproj_out_up$Emw),
                                        "Established type"=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Er + fproj_out_up$Erv + fproj_out_up$Erw),
                                        "New variant X"=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Em + fproj_out_up$Emv + fproj_out_up$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

O_mediumdf = fproj_out_base %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters_base["p"]*Total,
         Scenario = "Intermediate, beta=1.1  (baseline)" )
O_mediumdf <- rbind(O_mediumdf, fproj_out_down %>% dplyr::select(date, Total) %>%
                     mutate(sympinfect = fproj_parameters_down["p"]*Total,
                            Scenario = "Intermediate, beta=1.05" ))
O_mediumdf <- rbind(O_mediumdf, fproj_out_up %>% dplyr::select(date, Total) %>%
                     mutate(sympinfect = fproj_parameters_up["p"]*Total,
                            Scenario = "Intermediate, beta=1.15" ))


#### 3. Best case scenario (continued decreased trend, from ba.2-like parameters) --------
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.0,
                        "epsilon_m" = (1-0.4), # Previous mutant epsilon was 0.6
                        "c_m" = rem_parameters["c_m"]*1.0,
                        "c_mr" = rem_parameters["c_mr"]*1.0,
                        "c_rm" = rem_parameters["c_rm"]*1.0,
                        "w_m" =  rem_parameters["w_m"]*1.0,
                        "beffr" = 0.75,
                        "beffm"=0.75)
params_newmutant[[par_to_vary]] = rem_parameters[par_to_vary]*par_vals[7]
IHR_factor <- 1 # multiplier for IHR (see below)
res_swap = 0.3


# Swap resident and mutant, then set up new mutant -------------
new_model_base <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                               params_newmutant = params_newmutant, mut_prop = 0.1,
                               res_to_s_prop =  res_swap)
fproj_parameters_base <- new_model_base$newm_parameters

params_newmutant[[par_to_vary]] = rem_parameters[par_to_vary]*par_vals[8]
new_model_down <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                               params_newmutant = params_newmutant, mut_prop = 0.1,
                               res_to_s_prop =  res_swap)
fproj_parameters_down <- new_model_down$newm_parameters

params_newmutant[[par_to_vary]] = rem_parameters[par_to_vary]*par_vals[9]
new_model_up <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                             params_newmutant = params_newmutant, mut_prop = 0.1,
                             res_to_s_prop =  res_swap)
fproj_parameters_up <- new_model_up$newm_parameters

# Make projections
fproj_out_base <- as.data.frame(deSolve::ode(y=new_model_base$init_newm, time=times,func= sveirs,
                                             parms=fproj_parameters_base))
fproj_out_down <- as.data.frame(deSolve::ode(y=new_model_down$init_newm, time=times,func= sveirs,
                                             parms=fproj_parameters_down))
fproj_out_up <- as.data.frame(deSolve::ode(y=new_model_up$init_newm, time=times,func= sveirs,
                                           parms=fproj_parameters_up))

fproj_out_base <- fproj_out_base %>% mutate(Total=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Er + fproj_out_base$Erv + fproj_out_base$Erw +
                                                                                  fproj_out_base$Em + fproj_out_base$Emv + fproj_out_base$Emw),
                                            "Established type"=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Er + fproj_out_base$Erv + fproj_out_base$Erw),
                                            "New variant X"=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Em + fproj_out_base$Emv + fproj_out_base$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

fproj_out_down <- fproj_out_down %>% mutate(Total=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Er + fproj_out_down$Erv + fproj_out_down$Erw +
                                                                                  fproj_out_down$Em + fproj_out_down$Emv + fproj_out_down$Emw),
                                            "Established type"=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Er + fproj_out_down$Erv + fproj_out_down$Erw),
                                            "New variant X"=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Em + fproj_out_down$Emv + fproj_out_down$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

fproj_out_up <- fproj_out_up %>% mutate(Total=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Er + fproj_out_up$Erv + fproj_out_up$Erw +
                                                                            fproj_out_up$Em + fproj_out_up$Emv + fproj_out_up$Emw),
                                        "Established type"=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Er + fproj_out_up$Erv + fproj_out_up$Erw),
                                        "New variant X"=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Em + fproj_out_up$Emv + fproj_out_up$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

O_bestdf = fproj_out_base %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters_base["p"]*Total,
         Scenario = "Optimistic, beta=1.0 (baseline)" )
O_bestdf <- rbind(O_bestdf, fproj_out_down %>% dplyr::select(date, Total) %>%
                      mutate(sympinfect = fproj_parameters_down["p"]*Total,
                             Scenario = "Optimistic, beta=0.95" ))
O_bestdf <- rbind(O_bestdf, fproj_out_up %>% dplyr::select(date, Total) %>%
                      mutate(sympinfect = fproj_parameters_up["p"]*Total,
                             Scenario = "Optimistic, beta=1.05" ))



#### Figures ######################################
alldf  = rbind(O_bestdf, O_mediumdf, O_worstdf)
glimpse(alldf)
ggplot(filter(alldf, date<maxdate), aes(x=date, y = sympinfect/50, colour=Scenario, linetype=Scenario))+
  ylab("Incidence (symptomatic infection per 100K)") +
  geom_line(lwd=2,  alpha=0.5) + theme_minimal() +
  theme(legend.position = "bottom", legend.key.width=unit(3.8,"cm")) + xlab("") + guides(colour=guide_legend(nrow=3)) +
  scale_colour_manual(values=c("#34568B", "#92A8D1", "#6B5B95", "#9B2335", "#E15D44", "#BC243C", "#55B4B0", "#009B77", "#45B8AC")) +
  scale_linetype_manual(values=rep(c("dashed", "solid","twodash"),3)) + ggtitle(expression(paste("Sensitivity to ", beta[m])))

# now make one with hospitalizations
O_worstdf = O_worstdf %>% mutate(hosp = Total*IHR_BA2, census = lag(Total,8)*IHR_BA2*9) # no increase
O_mediumdf = O_mediumdf %>% mutate(hosp = Total*IHR_BA2, census = lag(Total,8)*IHR_BA2*9) # no increase
O_bestdf = O_bestdf %>% mutate(hosp = Total*IHR_BA2, census = lag(Total,8)*IHR_BA2*9) # no increase
alldf  = rbind(O_bestdf, O_mediumdf, O_worstdf)

ggplot(filter(alldf, date<maxdate), aes(x=date, y = hosp, colour=Scenario, linetype=Scenario))+
  ylab("COVID- all reported admissions") +
  geom_line(lwd=2,  alpha=0.5) + theme_minimal() +
  theme(legend.position = "bottom", legend.key.width=unit(3.8,"cm")) + xlab("") + guides(colour=guide_legend(nrow=3)) +
  scale_colour_manual(values=c("#34568B", "#92A8D1", "#6B5B95", "#9B2335", "#E15D44", "#BC243C", "#55B4B0", "#009B77", "#45B8AC")) +
  scale_linetype_manual(values=rep(c("solid", "dashed", "twodash"),3)) + ggtitle(expression(paste("Sensitivity to ", beta[m])))


ggplot(filter(alldf, date<maxdate), aes(x=date, y = census, colour=Scenario, linetype=Scenario))+
  ylab("COVID- census hospitalizations") +
  geom_line(lwd=2,  alpha=0.5) + theme_minimal() +
  theme(legend.position = "bottom", legend.key.width=unit(3.8,"cm")) + xlab("") + guides(colour=guide_legend(nrow=3)) +
  scale_colour_manual(values=c("#34568B", "#92A8D1", "#6B5B95", "#9B2335", "#E15D44", "#BC243C", "#55B4B0", "#009B77", "#45B8AC")) +
  scale_linetype_manual(values=rep(c("solid", "dashed", "twodash"),3)) + ggtitle(expression(paste("Sensitivity to ", beta[m])))







##### 1. Omicron - vary epsilon_m

par_to_vary <- "epsilon_m"
par_vals <- c(0.7, 0.65, 0.75, 0.7, 0.65, 0.75, 0.6, 0.55, 0.65) # worst, med, best (within = base, down, up)

# Worst case scenario - baseline ---------------------------
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.1, 
                        "epsilon_m" = (1-0.3), 
                        "c_m" = rem_parameters["c_m"]*1.0,
                        "c_mr" = rem_parameters["c_mr"]*1.0, 
                        "c_rm" = rem_parameters["c_rm"]*1.0,
                        "w_m" =  rem_parameters["w_m"]*1.0,
                        "beffr" = 0.7,
                        "beffm"=0.7)
params_newmutant[[par_to_vary]] = par_vals[1]
IHR_factor <- 1 # multiplier for IHR (see below)
res_swap = 0.3

# Swap resident and mutant, then set up new mutant -------------
new_model_base <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                               params_newmutant = params_newmutant, mut_prop = 0.1,
                               res_to_s_prop =  res_swap)
fproj_parameters_base <- new_model_base$newm_parameters

params_newmutant[[par_to_vary]] = par_vals[2]
new_model_down <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                               params_newmutant = params_newmutant, mut_prop = 0.1,
                               res_to_s_prop =  res_swap)
fproj_parameters_down <- new_model_down$newm_parameters

params_newmutant[[par_to_vary]] = par_vals[3]
new_model_up <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                             params_newmutant = params_newmutant, mut_prop = 0.1,
                             res_to_s_prop =  res_swap)
fproj_parameters_up <- new_model_up$newm_parameters

# Make projections
fproj_out_base <- as.data.frame(deSolve::ode(y=new_model_base$init_newm, time=times,func= sveirs,
                                             parms=fproj_parameters_base))
fproj_out_down <- as.data.frame(deSolve::ode(y=new_model_down$init_newm, time=times,func= sveirs,
                                             parms=fproj_parameters_down))
fproj_out_up <- as.data.frame(deSolve::ode(y=new_model_up$init_newm, time=times,func= sveirs,
                                           parms=fproj_parameters_up))

fproj_out_base <- fproj_out_base %>% mutate(Total=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Er + fproj_out_base$Erv + fproj_out_base$Erw +
                                                                                  fproj_out_base$Em + fproj_out_base$Emv + fproj_out_base$Emw),
                                            "Established type"=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Er + fproj_out_base$Erv + fproj_out_base$Erw),
                                            "New variant X"=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Em + fproj_out_base$Emv + fproj_out_base$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

fproj_out_down <- fproj_out_down %>% mutate(Total=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Er + fproj_out_down$Erv + fproj_out_down$Erw +
                                                                                  fproj_out_down$Em + fproj_out_down$Emv + fproj_out_down$Emw),
                                            "Established type"=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Er + fproj_out_down$Erv + fproj_out_down$Erw),
                                            "New variant X"=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Em + fproj_out_down$Emv + fproj_out_down$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

fproj_out_up <- fproj_out_up %>% mutate(Total=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Er + fproj_out_up$Erv + fproj_out_up$Erw +
                                                                            fproj_out_up$Em + fproj_out_up$Emv + fproj_out_up$Emw),
                                        "Established type"=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Er + fproj_out_up$Erv + fproj_out_up$Erw),
                                        "New variant X"=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Em + fproj_out_up$Emv + fproj_out_up$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

O_worstdf = fproj_out_base %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters_base["p"]*Total,
         Scenario = "Pessimistic, epsilon=0.7 (baseline)" )
O_worstdf <- rbind(O_worstdf, fproj_out_down %>% dplyr::select(date, Total) %>%
                     mutate(sympinfect = fproj_parameters_down["p"]*Total,
                            Scenario = "Pessimistic, epsilon=0.65" ))
O_worstdf <- rbind(O_worstdf, fproj_out_up %>% dplyr::select(date, Total) %>%
                     mutate(sympinfect = fproj_parameters_up["p"]*Total,
                            Scenario = "Pessimistic, epsilon=0.75" ))

#### 2. Intermediate case scenario (Increased beta & epsilon, gives ba.2 like wave) --------
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.1,
                        "epsilon_m" = (1-0.3), # Previous mutant epsilon was 0.6
                        "c_m" = rem_parameters["c_m"]*1.0,
                        "c_mr" = rem_parameters["c_mr"]*1.0,
                        "c_rm" = rem_parameters["c_rm"]*1.0,
                        "w_m" =  rem_parameters["w_m"]*1.0,
                        "beffr" = 0.75,
                        "beffm"=0.75)
params_newmutant[[par_to_vary]] = par_vals[4]
IHR_factor <- 1 # multiplier for IHR (see below)
res_swap = 0.3


# Swap resident and mutant, then set up new mutant -------------
new_model_base <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                               params_newmutant = params_newmutant, mut_prop = 0.1,
                               res_to_s_prop =  res_swap)
fproj_parameters_base <- new_model_base$newm_parameters

params_newmutant[[par_to_vary]] = par_vals[5]
new_model_down <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                               params_newmutant = params_newmutant, mut_prop = 0.1,
                               res_to_s_prop =  res_swap)
fproj_parameters_down <- new_model_down$newm_parameters

params_newmutant[[par_to_vary]] = par_vals[6]
new_model_up <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                             params_newmutant = params_newmutant, mut_prop = 0.1,
                             res_to_s_prop =  res_swap)
fproj_parameters_up <- new_model_up$newm_parameters

# Make projections
fproj_out_base <- as.data.frame(deSolve::ode(y=new_model_base$init_newm, time=times,func= sveirs,
                                             parms=fproj_parameters_base))
fproj_out_down <- as.data.frame(deSolve::ode(y=new_model_down$init_newm, time=times,func= sveirs,
                                             parms=fproj_parameters_down))
fproj_out_up <- as.data.frame(deSolve::ode(y=new_model_up$init_newm, time=times,func= sveirs,
                                           parms=fproj_parameters_up))

fproj_out_base <- fproj_out_base %>% mutate(Total=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Er + fproj_out_base$Erv + fproj_out_base$Erw +
                                                                                  fproj_out_base$Em + fproj_out_base$Emv + fproj_out_base$Emw),
                                            "Established type"=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Er + fproj_out_base$Erv + fproj_out_base$Erw),
                                            "New variant X"=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Em + fproj_out_base$Emv + fproj_out_base$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

fproj_out_down <- fproj_out_down %>% mutate(Total=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Er + fproj_out_down$Erv + fproj_out_down$Erw +
                                                                                  fproj_out_down$Em + fproj_out_down$Emv + fproj_out_down$Emw),
                                            "Established type"=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Er + fproj_out_down$Erv + fproj_out_down$Erw),
                                            "New variant X"=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Em + fproj_out_down$Emv + fproj_out_down$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

fproj_out_up <- fproj_out_up %>% mutate(Total=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Er + fproj_out_up$Erv + fproj_out_up$Erw +
                                                                            fproj_out_up$Em + fproj_out_up$Emv + fproj_out_up$Emw),
                                        "Established type"=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Er + fproj_out_up$Erv + fproj_out_up$Erw),
                                        "New variant X"=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Em + fproj_out_up$Emv + fproj_out_up$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

O_mediumdf = fproj_out_base %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters_base["p"]*Total,
         Scenario = "Intermediate, epsilon=0.7  (baseline)" )
O_mediumdf <- rbind(O_mediumdf, fproj_out_down %>% dplyr::select(date, Total) %>%
                      mutate(sympinfect = fproj_parameters_down["p"]*Total,
                             Scenario = "Intermediate, epsilon=0.65" ))
O_mediumdf <- rbind(O_mediumdf, fproj_out_up %>% dplyr::select(date, Total) %>%
                      mutate(sympinfect = fproj_parameters_up["p"]*Total,
                             Scenario = "Intermediate, epsilon=0.75" ))


#### 3. Best case scenario (continued decreased trend, from ba.2-like parameters) --------
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.0,
                        "epsilon_m" = (1-0.4), # Previous mutant epsilon was 0.6
                        "c_m" = rem_parameters["c_m"]*1.0,
                        "c_mr" = rem_parameters["c_mr"]*1.0,
                        "c_rm" = rem_parameters["c_rm"]*1.0,
                        "w_m" =  rem_parameters["w_m"]*1.0,
                        "beffr" = 0.75,
                        "beffm"=0.75)
params_newmutant[[par_to_vary]] = par_vals[7]
IHR_factor <- 1 # multiplier for IHR (see below)
res_swap = 0.3


# Swap resident and mutant, then set up new mutant -------------
new_model_base <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                               params_newmutant = params_newmutant, mut_prop = 0.1,
                               res_to_s_prop =  res_swap)
fproj_parameters_base <- new_model_base$newm_parameters

params_newmutant[[par_to_vary]] = par_vals[8]
new_model_down <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                               params_newmutant = params_newmutant, mut_prop = 0.1,
                               res_to_s_prop =  res_swap)
fproj_parameters_down <- new_model_down$newm_parameters

params_newmutant[[par_to_vary]] = par_vals[9]
new_model_up <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                             params_newmutant = params_newmutant, mut_prop = 0.1,
                             res_to_s_prop =  res_swap)
fproj_parameters_up <- new_model_up$newm_parameters

# Make projections
fproj_out_base <- as.data.frame(deSolve::ode(y=new_model_base$init_newm, time=times,func= sveirs,
                                             parms=fproj_parameters_base))
fproj_out_down <- as.data.frame(deSolve::ode(y=new_model_down$init_newm, time=times,func= sveirs,
                                             parms=fproj_parameters_down))
fproj_out_up <- as.data.frame(deSolve::ode(y=new_model_up$init_newm, time=times,func= sveirs,
                                           parms=fproj_parameters_up))

fproj_out_base <- fproj_out_base %>% mutate(Total=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Er + fproj_out_base$Erv + fproj_out_base$Erw +
                                                                                  fproj_out_base$Em + fproj_out_base$Emv + fproj_out_base$Emw),
                                            "Established type"=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Er + fproj_out_base$Erv + fproj_out_base$Erw),
                                            "New variant X"=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Em + fproj_out_base$Emv + fproj_out_base$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

fproj_out_down <- fproj_out_down %>% mutate(Total=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Er + fproj_out_down$Erv + fproj_out_down$Erw +
                                                                                  fproj_out_down$Em + fproj_out_down$Emv + fproj_out_down$Emw),
                                            "Established type"=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Er + fproj_out_down$Erv + fproj_out_down$Erw),
                                            "New variant X"=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Em + fproj_out_down$Emv + fproj_out_down$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

fproj_out_up <- fproj_out_up %>% mutate(Total=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Er + fproj_out_up$Erv + fproj_out_up$Erw +
                                                                            fproj_out_up$Em + fproj_out_up$Emv + fproj_out_up$Emw),
                                        "Established type"=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Er + fproj_out_up$Erv + fproj_out_up$Erw),
                                        "New variant X"=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Em + fproj_out_up$Emv + fproj_out_up$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

O_bestdf = fproj_out_base %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters_base["p"]*Total,
         Scenario = "Optimistic, epsilon=0.6 (baseline)" )
O_bestdf <- rbind(O_bestdf, fproj_out_down %>% dplyr::select(date, Total) %>%
                    mutate(sympinfect = fproj_parameters_down["p"]*Total,
                           Scenario = "Optimistic, epsilon=0.55" ))
O_bestdf <- rbind(O_bestdf, fproj_out_up %>% dplyr::select(date, Total) %>%
                    mutate(sympinfect = fproj_parameters_up["p"]*Total,
                           Scenario = "Optimistic, epsilon=0.65" ))



#### Figures ######################################
alldf  = rbind(O_bestdf, O_mediumdf, O_worstdf)
glimpse(alldf)
ggplot(filter(alldf, date<maxdate), aes(x=date, y = sympinfect/50, colour=Scenario, linetype=Scenario))+
  ylab("Incidence (symptomatic infection per 100K)") +
  geom_line(lwd=2,  alpha=0.5) + theme_minimal() +
  theme(legend.position = "bottom", legend.key.width=unit(3.8,"cm")) + xlab("") + guides(colour=guide_legend(nrow=3)) +
  scale_colour_manual(values=c("#34568B", "#92A8D1", "#6B5B95", "#9B2335", "#E15D44", "#BC243C", "#55B4B0", "#009B77", "#45B8AC")) +
  scale_linetype_manual(values=rep(c("dashed", "solid","twodash"),3)) + ggtitle(expression(paste("Sensitivity to ", epsilon[m])))











##### 1. Omicron - vary res_swap

par_to_vary <- "res_swap"
par_vals <- c(0.3, 0.25, 0.35 , 0.3, 0.25, 0.35 , 0.3, 0.25, 0.35) # worst, med, best (within = base, down, up)

# Worst case scenario - baseline ---------------------------
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.1, 
                        "epsilon_m" = (1-0.3), 
                        "c_m" = rem_parameters["c_m"]*1.0,
                        "c_mr" = rem_parameters["c_mr"]*1.0, 
                        "c_rm" = rem_parameters["c_rm"]*1.0,
                        "w_m" =  rem_parameters["w_m"]*1.0,
                        "beffr" = 0.7,
                        "beffm"=0.7)
IHR_factor <- 1 # multiplier for IHR (see below)
res_swap = par_vals[1]

# Swap resident and mutant, then set up new mutant -------------
new_model_base <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                               params_newmutant = params_newmutant, mut_prop = 0.1,
                               res_to_s_prop =  res_swap)
fproj_parameters_base <- new_model_base$newm_parameters

res_swap = par_vals[2]
new_model_down <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                               params_newmutant = params_newmutant, mut_prop = 0.1,
                               res_to_s_prop =  res_swap)
fproj_parameters_down <- new_model_down$newm_parameters

res_swap = par_vals[3]
new_model_up <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                             params_newmutant = params_newmutant, mut_prop = 0.1,
                             res_to_s_prop =  res_swap)
fproj_parameters_up <- new_model_up$newm_parameters

# Make projections
fproj_out_base <- as.data.frame(deSolve::ode(y=new_model_base$init_newm, time=times,func= sveirs,
                                             parms=fproj_parameters_base))
fproj_out_down <- as.data.frame(deSolve::ode(y=new_model_down$init_newm, time=times,func= sveirs,
                                             parms=fproj_parameters_down))
fproj_out_up <- as.data.frame(deSolve::ode(y=new_model_up$init_newm, time=times,func= sveirs,
                                           parms=fproj_parameters_up))

fproj_out_base <- fproj_out_base %>% mutate(Total=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Er + fproj_out_base$Erv + fproj_out_base$Erw +
                                                                                  fproj_out_base$Em + fproj_out_base$Emv + fproj_out_base$Emw),
                                            "Established type"=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Er + fproj_out_base$Erv + fproj_out_base$Erw),
                                            "New variant X"=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Em + fproj_out_base$Emv + fproj_out_base$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

fproj_out_down <- fproj_out_down %>% mutate(Total=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Er + fproj_out_down$Erv + fproj_out_down$Erw +
                                                                                  fproj_out_down$Em + fproj_out_down$Emv + fproj_out_down$Emw),
                                            "Established type"=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Er + fproj_out_down$Erv + fproj_out_down$Erw),
                                            "New variant X"=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Em + fproj_out_down$Emv + fproj_out_down$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

fproj_out_up <- fproj_out_up %>% mutate(Total=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Er + fproj_out_up$Erv + fproj_out_up$Erw +
                                                                            fproj_out_up$Em + fproj_out_up$Emv + fproj_out_up$Emw),
                                        "Established type"=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Er + fproj_out_up$Erv + fproj_out_up$Erw),
                                        "New variant X"=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Em + fproj_out_up$Emv + fproj_out_up$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

O_worstdf = fproj_out_base %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters_base["p"]*Total,
         Scenario = "Pessimistic, rho=0.3 (baseline)" )
O_worstdf <- rbind(O_worstdf, fproj_out_down %>% dplyr::select(date, Total) %>%
                     mutate(sympinfect = fproj_parameters_down["p"]*Total,
                            Scenario = "Pessimistic, rho=0.25" ))
O_worstdf <- rbind(O_worstdf, fproj_out_up %>% dplyr::select(date, Total) %>%
                     mutate(sympinfect = fproj_parameters_up["p"]*Total,
                            Scenario = "Pessimistic, rho=0.35" ))

#### 2. Intermediate case scenario (Increased beta & epsilon, gives ba.2 like wave) --------
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.1,
                        "epsilon_m" = (1-0.3), # Previous mutant epsilon was 0.6
                        "c_m" = rem_parameters["c_m"]*1.0,
                        "c_mr" = rem_parameters["c_mr"]*1.0,
                        "c_rm" = rem_parameters["c_rm"]*1.0,
                        "w_m" =  rem_parameters["w_m"]*1.0,
                        "beffr" = 0.75,
                        "beffm"=0.75)
res_swap  = par_vals[4]
IHR_factor <- 1 # multiplier for IHR (see below)


# Swap resident and mutant, then set up new mutant -------------
new_model_base <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                               params_newmutant = params_newmutant, mut_prop = 0.1,
                               res_to_s_prop =  res_swap)
fproj_parameters_base <- new_model_base$newm_parameters

res_swap = par_vals[5]
new_model_down <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                               params_newmutant = params_newmutant, mut_prop = 0.1,
                               res_to_s_prop =  res_swap)
fproj_parameters_down <- new_model_down$newm_parameters

res_swap = par_vals[6]
new_model_up <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                             params_newmutant = params_newmutant, mut_prop = 0.1,
                             res_to_s_prop =  res_swap)
fproj_parameters_up <- new_model_up$newm_parameters

# Make projections
fproj_out_base <- as.data.frame(deSolve::ode(y=new_model_base$init_newm, time=times,func= sveirs,
                                             parms=fproj_parameters_base))
fproj_out_down <- as.data.frame(deSolve::ode(y=new_model_down$init_newm, time=times,func= sveirs,
                                             parms=fproj_parameters_down))
fproj_out_up <- as.data.frame(deSolve::ode(y=new_model_up$init_newm, time=times,func= sveirs,
                                           parms=fproj_parameters_up))

fproj_out_base <- fproj_out_base %>% mutate(Total=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Er + fproj_out_base$Erv + fproj_out_base$Erw +
                                                                                  fproj_out_base$Em + fproj_out_base$Emv + fproj_out_base$Emw),
                                            "Established type"=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Er + fproj_out_base$Erv + fproj_out_base$Erw),
                                            "New variant X"=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Em + fproj_out_base$Emv + fproj_out_base$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

fproj_out_down <- fproj_out_down %>% mutate(Total=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Er + fproj_out_down$Erv + fproj_out_down$Erw +
                                                                                  fproj_out_down$Em + fproj_out_down$Emv + fproj_out_down$Emw),
                                            "Established type"=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Er + fproj_out_down$Erv + fproj_out_down$Erw),
                                            "New variant X"=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Em + fproj_out_down$Emv + fproj_out_down$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

fproj_out_up <- fproj_out_up %>% mutate(Total=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Er + fproj_out_up$Erv + fproj_out_up$Erw +
                                                                            fproj_out_up$Em + fproj_out_up$Emv + fproj_out_up$Emw),
                                        "Established type"=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Er + fproj_out_up$Erv + fproj_out_up$Erw),
                                        "New variant X"=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Em + fproj_out_up$Emv + fproj_out_up$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

O_mediumdf = fproj_out_base %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters_base["p"]*Total,
         Scenario = "Intermediate, epsilon=0.3  (baseline)" )
O_mediumdf <- rbind(O_mediumdf, fproj_out_down %>% dplyr::select(date, Total) %>%
                      mutate(sympinfect = fproj_parameters_down["p"]*Total,
                             Scenario = "Intermediate, epsilon=0.25" ))
O_mediumdf <- rbind(O_mediumdf, fproj_out_up %>% dplyr::select(date, Total) %>%
                      mutate(sympinfect = fproj_parameters_up["p"]*Total,
                             Scenario = "Intermediate, epsilon=0.35" ))


#### 3. Best case scenario (continued decreased trend, from ba.2-like parameters) --------
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.0,
                        "epsilon_m" = (1-0.4), # Previous mutant epsilon was 0.6
                        "c_m" = rem_parameters["c_m"]*1.0,
                        "c_mr" = rem_parameters["c_mr"]*1.0,
                        "c_rm" = rem_parameters["c_rm"]*1.0,
                        "w_m" =  rem_parameters["w_m"]*1.0,
                        "beffr" = 0.75,
                        "beffm"=0.75)
res_swap = par_vals[7]
IHR_factor <- 1 # multiplier for IHR (see below)


# Swap resident and mutant, then set up new mutant -------------
new_model_base <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                               params_newmutant = params_newmutant, mut_prop = 0.1,
                               res_to_s_prop =  res_swap)
fproj_parameters_base <- new_model_base$newm_parameters

res_swap  = par_vals[8]
new_model_down <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                               params_newmutant = params_newmutant, mut_prop = 0.1,
                               res_to_s_prop =  res_swap)
fproj_parameters_down <- new_model_down$newm_parameters

res_swap  = par_vals[9]
new_model_up <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                             params_newmutant = params_newmutant, mut_prop = 0.1,
                             res_to_s_prop =  res_swap)
fproj_parameters_up <- new_model_up$newm_parameters

# Make projections
fproj_out_base <- as.data.frame(deSolve::ode(y=new_model_base$init_newm, time=times,func= sveirs,
                                             parms=fproj_parameters_base))
fproj_out_down <- as.data.frame(deSolve::ode(y=new_model_down$init_newm, time=times,func= sveirs,
                                             parms=fproj_parameters_down))
fproj_out_up <- as.data.frame(deSolve::ode(y=new_model_up$init_newm, time=times,func= sveirs,
                                           parms=fproj_parameters_up))

fproj_out_base <- fproj_out_base %>% mutate(Total=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Er + fproj_out_base$Erv + fproj_out_base$Erw +
                                                                                  fproj_out_base$Em + fproj_out_base$Emv + fproj_out_base$Emw),
                                            "Established type"=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Er + fproj_out_base$Erv + fproj_out_base$Erw),
                                            "New variant X"=
                                              fproj_parameters_base[["sigma"]]*(fproj_out_base$Em + fproj_out_base$Emv + fproj_out_base$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

fproj_out_down <- fproj_out_down %>% mutate(Total=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Er + fproj_out_down$Erv + fproj_out_down$Erw +
                                                                                  fproj_out_down$Em + fproj_out_down$Emv + fproj_out_down$Emw),
                                            "Established type"=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Er + fproj_out_down$Erv + fproj_out_down$Erw),
                                            "New variant X"=
                                              fproj_parameters_down[["sigma"]]*(fproj_out_down$Em + fproj_out_down$Emv + fproj_out_down$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

fproj_out_up <- fproj_out_up %>% mutate(Total=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Er + fproj_out_up$Erv + fproj_out_up$Erw +
                                                                            fproj_out_up$Em + fproj_out_up$Emv + fproj_out_up$Emw),
                                        "Established type"=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Er + fproj_out_up$Erv + fproj_out_up$Erw),
                                        "New variant X"=
                                          fproj_parameters_up[["sigma"]]*(fproj_out_up$Em + fproj_out_up$Emv + fproj_out_up$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-09-01"),ymd("2022-09-01")+300, 1))

O_bestdf = fproj_out_base %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters_base["p"]*Total,
         Scenario = "Optimistic, epsilon=0.3 (baseline)" )
O_bestdf <- rbind(O_bestdf, fproj_out_down %>% dplyr::select(date, Total) %>%
                    mutate(sympinfect = fproj_parameters_down["p"]*Total,
                           Scenario = "Optimistic, epsilon=0.25" ))
O_bestdf <- rbind(O_bestdf, fproj_out_up %>% dplyr::select(date, Total) %>%
                    mutate(sympinfect = fproj_parameters_up["p"]*Total,
                           Scenario = "Optimistic, epsilon=0.35" ))



#### Figures ######################################
alldf  = rbind(O_bestdf, O_mediumdf, O_worstdf)
glimpse(alldf)
ggplot(filter(alldf, date<maxdate), aes(x=date, y = sympinfect/50, colour=Scenario, linetype=Scenario))+
  ylab("Incidence (symptomatic infection per 100K)") +
  geom_line(lwd=2,  alpha=0.5) + theme_minimal() +
  theme(legend.position = "bottom", legend.key.width=unit(3.8,"cm")) + xlab("") + guides(colour=guide_legend(nrow=3)) +
  scale_colour_manual(values=c("#34568B", "#92A8D1", "#6B5B95", "#9B2335", "#E15D44", "#BC243C", "#55B4B0", "#009B77", "#45B8AC")) +
  scale_linetype_manual(values=rep(c("dashed", "solid","twodash"),3)) + ggtitle(expression(paste("Sensitivity to ", rho)))















#### Find various parameter sets that give the pessimistic scenario

#### 1. Baseline worst case scenario:  ---------------------------
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

# Swap resident and mutant, then set up new mutant 
new_model <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                          params_newmutant = params_newmutant, mut_prop = 0.1,
                          res_to_s_prop =  res_swap)
fproj_parameters <- new_model$newm_parameters
# Make projections
fproj_out <- as.data.frame(deSolve::ode(y=new_model$init_newm, time=times,func= sveirs,
                                        parms=fproj_parameters))

# Plot of the projected cases -------------
fproj_out <- fproj_out %>% mutate(Total=
                                    fproj_parameters[["sigma"]]*(fproj_out$Er + fproj_out$Erv + fproj_out$Erw +
                                                                   fproj_out$Em + fproj_out$Emv + fproj_out$Emw),
                                  "Established type"=
                                    fproj_parameters[["sigma"]]*(fproj_out$Er + fproj_out$Erv + fproj_out$Erw),
                                  "New variant X"=
                                    fproj_parameters[["sigma"]]*(fproj_out$Em + fproj_out$Emv + fproj_out$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-10-01"),ymd("2022-10-01")+300, 1))

NO_worstdf = fproj_out %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters["p"]*Total,
         Scenario = "Baseline pessimistic")



##########################################################
#### 2. A first match to pessimistic scenario  ---------------------------
params_newmutant = list("beta_r" = rem_parameters["beta_r"]*0.9, # compensate for artificial boost
                        "beta_m" = rem_parameters["beta_m"]*1.4,
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
res_swap = 0.5

# Swap resident and mutant, then set up new mutant 
new_model <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                          params_newmutant = params_newmutant, mut_prop = 0.1,
                          res_to_s_prop =  res_swap)
fproj_parameters <- new_model$newm_parameters

# Make projections
fproj_out <- as.data.frame(deSolve::ode(y=new_model$init_newm, time=times,func= sveirs,
                                        parms=fproj_parameters))

# Plot of the projected cases -------------
fproj_out <- fproj_out %>% mutate(Total=
                                    fproj_parameters[["sigma"]]*(fproj_out$Er + fproj_out$Erv + fproj_out$Erw +
                                                                   fproj_out$Em + fproj_out$Emv + fproj_out$Emw),
                                  "Established type"=
                                    fproj_parameters[["sigma"]]*(fproj_out$Er + fproj_out$Erv + fproj_out$Erw),
                                  "New variant X"=
                                    fproj_parameters[["sigma"]]*(fproj_out$Em + fproj_out$Emv + fproj_out$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-10-01"),ymd("2022-10-01")+300, 1))

NO_alt1df = fproj_out %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters["p"]*Total,
         Scenario = "Alt 1" )


#### 3. A second match to pessimistic scenario  ---------------------------
params_newmutant = list("beta_r" = rem_parameters["beta_r"]*0.9, # compensate for artificial boost
                        "beta_m" = rem_parameters["beta_m"]*1.5,
                        "gamma"=1/4,
                        "sigma"=1,
                        #    "epsilon_m" = (1-0.4), for ref this is what i did for ba2
                        "epsilon_m" = (1-0.4), #
                        "c_m" = rem_parameters["c_m"]*1,
                        "c_mr" = rem_parameters["c_mr"]*1.0,
                        "c_rm" = rem_parameters["c_rm"]*1,
                        "w_m" =  rem_parameters["w_m"]*1,
                        "beffr" = 0.75,
                        "beffm"=0.75)
IHR_factor <- 2 # multiplier for IHR (see below)
res_swap = 0.7

# Swap resident and mutant, then set up new mutant -------------
new_model <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                          params_newmutant = params_newmutant, mut_prop = 0.1,
                          res_to_s_prop =  res_swap)
fproj_parameters <- new_model$newm_parameters

# Make projections
fproj_out <- as.data.frame(deSolve::ode(y=new_model$init_newm, time=times,func= sveirs,
                                        parms=fproj_parameters))

# Plot of the projected cases -------------
fproj_out <- fproj_out %>% mutate(Total=
                                    fproj_parameters[["sigma"]]*(fproj_out$Er + fproj_out$Erv + fproj_out$Erw +
                                                                   fproj_out$Em + fproj_out$Emv + fproj_out$Emw),
                                  "Established type"=
                                    fproj_parameters[["sigma"]]*(fproj_out$Er + fproj_out$Erv + fproj_out$Erw),
                                  "New variant X"=
                                    fproj_parameters[["sigma"]]*(fproj_out$Em + fproj_out$Emv + fproj_out$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-10-01"),ymd("2022-10-01")+300, 1))

NO_alt2df = fproj_out %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters["p"]*Total,
         Scenario = "Alt 2" )


#### 4. A third match to pessimistic scenario  --------
params_newmutant = list("beta_r" = rem_parameters["beta_r"]*0.9, # compensate for artificial boost
                        "beta_m" = rem_parameters["beta_m"]*1.8,
                        "gamma"=1/4,
                        "sigma"=1,
                        #    "epsilon_m" = (1-0.4), for ref this is what i did for ba2
                        "epsilon_m" = (1-0.3), #
                        "c_m" = rem_parameters["c_m"]*1,
                        "c_mr" = rem_parameters["c_mr"]*1.1,
                        "c_rm" = rem_parameters["c_rm"]*1,
                        "w_m" =  rem_parameters["w_m"]*1,
                        "beffr" = 0.7,
                        "beffm"=0.7)
IHR_factor <- 1 # multiplier for IHR (see below)
res_swap = 0.7

# Swap resident and mutant, then set up new mutant -------------
new_model <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                          params_newmutant = params_newmutant, mut_prop = 0.03,
                          res_to_s_prop =  res_swap)
fproj_parameters <- new_model$newm_parameters
# Make projections
fproj_out <- as.data.frame(deSolve::ode(y=new_model$init_newm, time=times,func= sveirs,
                                        parms=fproj_parameters))

# Plot of the projected cases -------------
fproj_out <- fproj_out %>% mutate(Total=
                                    fproj_parameters[["sigma"]]*(fproj_out$Er + fproj_out$Erv + fproj_out$Erw +
                                                                   fproj_out$Em + fproj_out$Emv + fproj_out$Emw),
                                  "Established type"=
                                    fproj_parameters[["sigma"]]*(fproj_out$Er + fproj_out$Erv + fproj_out$Erw),
                                  "New variant X"=
                                    fproj_parameters[["sigma"]]*(fproj_out$Em + fproj_out$Emv + fproj_out$Emw)) %>%
  mutate(date=seq.Date(ymd("2022-10-01"),ymd("2022-10-01")+300, 1))

NO_alt3df = fproj_out %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters["p"]*Total,
         Scenario = "Alt 3" )





#### Make figure(s)

NO_worstdf = NO_worstdf %>% mutate(hosp = Total*IHR_BA2*2, census = lag(Total,8)*IHR_BA2*2*9) # 2x increase
NO_alt1df = NO_alt1df %>% mutate(hosp = Total*IHR_BA2*2, census = lag(Total,8)*IHR_BA2*2*9) # 2x increase
NO_alt2df = NO_alt2df %>% mutate(hosp = Total*IHR_BA2*2, census = lag(Total,8)*IHR_BA2*2*9) # 2x increase
NO_alt3df = NO_alt3df %>% mutate(hosp = Total*IHR_BA2*1, census = lag(Total,8)*IHR_BA2*1*9) # no increase
alldf  = rbind(NO_worstdf, NO_alt1df, NO_alt2df, NO_alt3df)

ggplot(filter(alldf, date<maxdate), aes(x=date, y = hosp, fill=Scenario))+
  ylab("COVID- all reported admissions") +
  geom_line(data = filter(alldf, date<maxdate & Scenario=="Baseline pessimistic"), lty = "dashed", lwd=1.3) +
  geom_ribbon(aes(x=date, ymin = 0.75*hosp, ymax=1.25*hosp,
                  fill=Scenario),  alpha=0.5) + theme_minimal() +
  theme(legend.position = "bottom") +  
  xlab("") +  
  scale_x_continuous(breaks= c(as.Date("2022-10-01"),as.Date("2022-11-01"),as.Date("2022-12-01"),as.Date("2023-01-01"),as.Date("2023-02-01")), labels = c("Month 1","Month 2","Month 3","Month 4","Month 5")) 


