
# This script follows on directly from pipps_projection.r
load("projectionscript_out_JS.Rdata")

# Adjust dates for projection
old_intro_date  = intro_date # Keep track of 'day 0'
intro_date <-   ymd("2022-06-20") # NOTE MOVED TO JUNE 
proj_out = filter(proj_out, date<= intro_date) # sets last part of trajectory to the intro date to get the last state, for the IC of the new run 
initscale=1.5
stop_date <- intro_date + 300 # 300 days projection
times = as.numeric(intro_date - old_intro_date):(as.numeric(stop_date - old_intro_date)) # Keep the same time count going

rem_parameters  <- proj_parameters # ...to keep parameters safe
# Set the desired characteristics of the new mutant. You can include any of the named elements of rem_parameters 
# Choose which scenario you want to run:

# Which date to plot projections out until
maxdate = ymd("2022-12-30")
# when do you want to start the date for the model? 
startdate = ymd("2022-09-15")
IHRcdn = 0.00484 # note: Canada-wide IHR
#COVID-19 resources Canada (run by Guillaume Bourque and Tara Moriarty. https://covid19resources.ca/covid-hazard-index/) provide:
#  - Data: reported Canadian hospitalizations since Dec 2nd 2021: 116115
#- Estimates: Canadian cumulative % population infected since Dec 2nd 2021: 62.04% (as of 16/09/22. Calculated from blood donation seroprevalence)
#Canada population: 38,654,738 on Apr 1 2022 https://www150.statcan.gc.ca/n1/daily-quotidien/220622/dq220622d-eng.htm
# 116115/(38654738*0.6204) = 0.00484

#### 1. BASELINE Omicron scenario 
params_newmutant = list("beta_r" = rem_parameters["beta_r"]*1.15,
                        "beta_m" = rem_parameters["beta_m"]*1.18,
                        "epsilon_m" = (1-0.325), # Previous mutant epsilon was 0.6
                        "epsilon_r" = (1-0.325),
                        "c_m" = rem_parameters["c_m"]*1.0,
                        "c_mr" = rem_parameters["c_mr"]*1.0,
                        "c_rm" = rem_parameters["c_rm"]*1.0,
                        "w_m" =  rem_parameters["w_m"]*1.0,
                      #  "b" = 0.015, # orig was 0.015 
                      #  "w_b" = 1/(0.5*365),
                        "beffr" = 0.7,
                        "beffm"=0.7) # boosted NUMBER stays the same at about 0.5

IHR_factor <- 1 # multiplier for IHR (see below)
res_swap = 0.1

# swap pars
new_model <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                          params_newmutant = params_newmutant, mut_prop = 0.1,
                          res_to_s_prop =  res_swap)
fproj_parameters <- new_model$newm_parameters

# visualize and fiddle with IC
x =new_model$init_newm
pie(as.numeric(x),  
    labels = colnames(out_samp)[2:ncol(out_samp)], radius = 1.0)

newstate = tweak_init(oldstate=x, myfactor=initscale) 
pie(as.numeric(newstate),  
    labels = colnames(out_samp)[2:ncol(out_samp)], radius = 1.0)


# Make projection
fproj_out <- as.data.frame(deSolve::ode(y=newstate, time=times,func= sveirs,
                                        parms=fproj_parameters))
get_growth_rate(fproj_out, startoffset = 5, duration = 10)



# check some features of the trajectory 

# check trajectory
fproj_out <- fproj_out %>% mutate(Total=
                                    fproj_parameters[["sigma"]]*(fproj_out$Er + fproj_out$Erv + fproj_out$Erw +
                                                                   fproj_out$Em + fproj_out$Emv + fproj_out$Emw),
                                  "Established type"=
                                    fproj_parameters[["sigma"]]*(fproj_out$Er + fproj_out$Erv + fproj_out$Erw),
                                  "New variant X"=
                                    fproj_parameters[["sigma"]]*(fproj_out$Em + fproj_out$Emv + fproj_out$Emw)) %>%
  mutate(date=seq.Date(startdate,startdate+300, 1))

baseline.df = fproj_out %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters["p"]*Total,
         Scenario = "Baseline")


#### 2. MORE BOOSTING than the baseline scenario. raise beffm i think. 
params_newmutant = list("beta_r" = rem_parameters["beta_r"]*1.15,
                        "beta_m" = rem_parameters["beta_m"]*1.18,
                        "epsilon_m" = (1-0.325), # Previous mutant epsilon was 0.6
                        "epsilon_r" = (1-0.325),
                        "c_m" = rem_parameters["c_m"]*1.0,
                        "c_mr" = rem_parameters["c_mr"]*1.0,
                        "c_rm" = rem_parameters["c_rm"]*1.0,
                        "w_m" =  rem_parameters["w_m"]*1.0,
                        #  "b" = 0.015, # orig was 0.015 
                        #  "w_b" = 1/(0.5*365),
                        "beffr" = 0.7,
                        "beffm"=0.8) # boosted NUMBER stays the same at about 0.5

IHR_factor <- 1 # multiplier for IHR (see below)
res_swap = 0.1

# swap pars
new_model <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                          params_newmutant = params_newmutant, mut_prop = 0.1,
                          res_to_s_prop =  res_swap)
fproj_parameters <- new_model$newm_parameters

# visualize and fiddle with IC
x =new_model$init_newm
pie(as.numeric(x),  
    labels = colnames(out_samp)[2:ncol(out_samp)], radius = 1.0)

newstate = tweak_init(oldstate=x,  myfactor=initscale) 
pie(as.numeric(newstate),  
    labels = colnames(out_samp)[2:ncol(out_samp)], radius = 1.0)


# Make projection
fproj_out <- as.data.frame(deSolve::ode(y=newstate, time=times,func= sveirs,
                                        parms=fproj_parameters))
get_growth_rate(fproj_out, startoffset = 5, duration = 10)



# check some features of the trajectory 

# check trajectory
fproj_out <- fproj_out %>% mutate(Total=
                                    fproj_parameters[["sigma"]]*(fproj_out$Er + fproj_out$Erv + fproj_out$Erw +
                                                                   fproj_out$Em + fproj_out$Emv + fproj_out$Emw),
                                  "Established type"=
                                    fproj_parameters[["sigma"]]*(fproj_out$Er + fproj_out$Erv + fproj_out$Erw),
                                  "New variant X"=
                                    fproj_parameters[["sigma"]]*(fproj_out$Em + fproj_out$Emv + fproj_out$Emw)) %>%
  mutate(date=seq.Date(startdate,startdate+300, 1))


moreboost.df = fproj_out %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters["p"]*Total,
         Scenario = "More boosting")


#### 2. LESS BOOSTING than the baseline scenario. raise beffm i think. 
params_newmutant = list("beta_r" = rem_parameters["beta_r"]*1.15,
                        "beta_m" = rem_parameters["beta_m"]*1.18,
                        "epsilon_m" = (1-0.325), # Previous mutant epsilon was 0.6
                        "epsilon_r" = (1-0.325),
                        "c_m" = rem_parameters["c_m"]*1.0,
                        "c_mr" = rem_parameters["c_mr"]*1.0,
                        "c_rm" = rem_parameters["c_rm"]*1.0,
                        "w_m" =  rem_parameters["w_m"]*1.0,
                        #  "b" = 0.015, # orig was 0.015 
                        #  "w_b" = 1/(0.5*365),
                        "beffr" = 0.7,
                        "beffm"=0.6) # boosted NUMBER stays the same at about 0.5

IHR_factor <- 1 # multiplier for IHR (see below)
res_swap = 0.1

# swap pars
new_model <- swap_strains(out_old = proj_out, params_old = rem_parameters,
                          params_newmutant = params_newmutant, mut_prop = 0.1,
                          res_to_s_prop =  res_swap)
fproj_parameters <- new_model$newm_parameters

# visualize and fiddle with IC
x =new_model$init_newm
pie(as.numeric(x),  
    labels = colnames(out_samp)[2:ncol(out_samp)], radius = 1.0)

newstate = tweak_init(oldstate=x,  myfactor=initscale) 
pie(as.numeric(newstate),  
    labels = colnames(out_samp)[2:ncol(out_samp)], radius = 1.0)


# Make projection
fproj_out <- as.data.frame(deSolve::ode(y=newstate, time=times,func= sveirs,
                                        parms=fproj_parameters))
get_growth_rate(fproj_out, startoffset = 5, duration = 10)



# check some features of the trajectory 

# check trajectory
fproj_out <- fproj_out %>% mutate(Total=
                                    fproj_parameters[["sigma"]]*(fproj_out$Er + fproj_out$Erv + fproj_out$Erw +
                                                                   fproj_out$Em + fproj_out$Emv + fproj_out$Emw),
                                  "Established type"=
                                    fproj_parameters[["sigma"]]*(fproj_out$Er + fproj_out$Erv + fproj_out$Erw),
                                  "New variant X"=
                                    fproj_parameters[["sigma"]]*(fproj_out$Em + fproj_out$Emv + fproj_out$Emw)) %>%
  mutate(date=seq.Date(startdate+5,startdate+5+300, 1))


lessboost.df = fproj_out %>% dplyr::select(date, Total) %>%
  mutate(sympinfect = fproj_parameters["p"]*Total,
         Scenario = "Less boosting")

alldf = rbind(baseline.df, moreboost.df, lessboost.df)


# ---- plotting section 

# plot all projections together 

# incident infections in Canada 
ggplot(filter(alldf, date <= maxdate), aes(x=date, y=Total, fill=Scenario))+
  geom_ribbon(aes(x=date, ymin = 0.75*(38/5)*(100000/38e6)*Total, ymax=1.25*(38/5)*(100000/38e6)*Total), 
               alpha=0.45)+  scale_fill_brewer(palette="Set1")+ xlab("")+
  ylab("Daily new infections per 100K") + theme_minimal()
ggsave("boost-scenarios-infections.pdf") 

# total hosps in Canada based on an 8 day stay and the IHR from public data, above
ggplot(filter(alldf, date <= maxdate), aes(x=date, y=Total, fill=Scenario))+
  geom_ribbon(aes(x=date, ymin = IHRcdn*8*0.75*(38/5)*Total, 
                  ymax=IHRcdn*8*1.25*(38/5)*Total), 
              alpha=0.45)+  scale_fill_brewer(palette="Set1")+ xlab("")+
  ylab("Total hospitalizations") + 
  geom_point(data=filter(hosp, Date > min(fproj_out$date-250)), 
             aes(x=Date,y=COVID_HOSP),inherit.aes = F)+xlab("") + theme_minimal()
theme_minimal()
ggsave("boost-scenarios-hosp.pdf") 


# plots of one projection at a time
pivot_longer(filter(fproj_out, date <= maxdate), c(Total,"Established type", "New variant X"), names_to = "Strain", values_to = "count") %>%
  ggplot(aes(x=date, y=count, colour=Strain)) + geom_line() + ylab("Incident infections") + xlab("Date") + theme_minimal()

# first go read in the hosp data, below 
pivot_longer(filter(fproj_out, date <= maxdate), c(Total,"Established type", "New variant X"), names_to = "Strain", values_to = "count") %>%
  ggplot(aes(x=date, y=count*IHRcdn*8*(38/5), colour=Strain)) + geom_line() + ylab("Total hospitalizations") +
  geom_point(data=filter(hosp, Date > min(fproj_out$date-250)), aes(x=Date,y=COVID_HOSP),inherit.aes = F)+xlab("Date") + theme_minimal()



# --- hosp data - public, from https://health-infobase.canada.ca/covid-19/
hosp = readr::read_csv("data/covid19-epiSummary-hospVentICU.csv")
glimpse(hosp)



