# now use Jessica's new code to add a different kind of BA.2 wave to the end of the BA.1 wave 
# and see if we can get it doing the right thing without swapping all the BA.1
# infections to BA.2 infections.

# instead, JS code will make the new m strain BA2 and will let us add 
# its own parameters

# i learned that I can fix p = .2 and get just as good a fit as we had before

# this follws on from line 242 of pipps_simulation.R 
out_samp$date = intro_date+out_samp$time

params_newmutant = list("beta_m" = parameters["beta_m"]*1.7, "eff_t" = -10) # eff_t is just pushed back beyond the time horizon of the projection
new_model <- swap_strains(out_old = out_samp, params_old = parameters, 
                          params_newmutant = params_newmutant, mut_prop = 0.3)
init_proj <- new_model$init_newm
proj_parameters <- new_model$newm_parameters
forecasts_days <- 200
times <- 1:(forecasts_days)
proj_out <- as.data.frame(deSolve::ode(y=init_proj, time=times,func= sveirs,
                                       parms=proj_parameters)) 

#check growth rate 
get_growth_rate(output= proj_out, startoffset = 20, duration = 7)

proj_out <- proj_out %>% mutate(Total=last(test_prop)*proj_parameters[["p"]]*
                                  proj_parameters[["sigma"]]*(proj_out$Er + proj_out$Erv + proj_out$Erw + 
                                                                proj_out$Em + proj_out$Emv + proj_out$Emw), 
                                Resident=last(test_prop)*proj_parameters[["p"]]*
                                  proj_parameters[["sigma"]]*(proj_out$Er + proj_out$Erv + proj_out$Erw), 
                                Mutant=last(test_prop)*proj_parameters[["p"]]*
                                  proj_parameters[["sigma"]]*(proj_out$Em + proj_out$Emv + proj_out$Emw)) %>% 
  mutate(date=seq.Date(ymd("2022-03-31"),ymd("2022-03-31")-1+length(times), 1)) 


pivot_longer(proj_out, c(Total, Resident, Mutant), names_to = "Strain", values_to = "count") %>%
  ggplot(aes(x=date, y=count, colour=Strain)) + geom_line()

tmp = pivot_longer(proj_out, c(Total, Resident, Mutant), 
                   names_to = "Strain", values_to = "count") %>%
  dplyr::select(date, Strain, count)

dat_full = readRDS("data/BC-dat.rds")

ggplot() + geom_line(data=project_dat_BC,aes(x=date,y=`50%`), col="green",size=1.5,alpha=0.4) + 
  geom_point(data=filter(dat_full,date>= intro_date),aes(x=date, y=value),color='grey48', alpha=0.8, size = 1.5) + 
  geom_line(data =tmp, aes(x=date, y=count, colour=Strain), inherit.aes = F)
