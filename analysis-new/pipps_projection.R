#add dat to data for plotting 
#here I match model output to the remaining data point (i.e. from March 31, 2022 onward)
dat_full = readRDS("data/BC-dat.rds")

forecasts_days <- 1 
intro_date <-   ymd("2022-03-13")
stop_date <- last(dat_full$date)

dat_rem <- dat_full %>% filter(date >= intro_date &  date <= stop_date)

#dat_rem <- filter(dat_full, date >= intro_date) %>% select(c("day", "value", "date"))
dat_full$day <- 1:nrow(dat_full)

plot(dat_rem$value)

#initiating with the last point on fit from pipps_simulation.R script 
init_rem <-     c(S=last(out_samp$S),
                Er=1,Em=last(out_samp$Em)-1,
                Ir=2,Im=last(out_samp$Im)-2,
                R=last(out_samp$R),V=last(out_samp$V),
                Erv=10,Emv=last(out_samp$Emv)-10,
                Irv=1,Imv=last(out_samp$Imv)-1,
                Rv=last(out_samp$Rv),W=last(out_samp$W),
                Erw=0,Emw=last(out_samp$Emw),
                Irw=0,Imw=last(out_samp$Imw)-0
                ,Rw=last(out_samp$Rw)) 




rem_parameters  <- parameters

rem_parameters["beta_r"] <- rem_parameters["beta_m"]*(1-rem_parameters[["stngcy"]])*1.5*1.35
rem_parameters["beta_m"] <- rem_parameters["beta_m"]*(1-rem_parameters[["stngcy"]])*1.5
rem_parameters["eff_t"]  <- 1000 #set to  some time in the future beyond  the projection period 
rem_parameters["epsilon_r"] <- (1-0.15) 
#rem_parameters[["stngcy"]] <- 0.35
#rem_parameters[["p"]] <- 0.5
#rem_parameters[["beff"]] <- 0.95
#rem_parameters[["wf"]] <- 0.01
rem_parameters[["b"]] <- 0.018*1.3
#initialm data matching 




forecasts_days <- 1
times <- 1:(nrow(dat_rem) + forecasts_days)
rem_outtest <- as.data.frame(deSolve::ode(y=init_rem,time=times,func= sveirs,
                             parms=rem_parameters))  
rem_outtest <- rem_outtest %>% mutate(incid=last(test_prop)*rem_parameters[["p"]]*
                             rem_parameters[["sigma"]]*(rem_outtest$Er + rem_outtest$Erv + rem_outtest$Erw + 
                              rem_outtest$Em + rem_outtest$Emv + rem_outtest$Emw)) %>% 
                             mutate(date=seq.Date(ymd(intro_date),ymd(intro_date)-1+length(times), 1))

ggplot(data =rem_outtest, aes(x=date, y = incid))+geom_line() +
geom_line(data=dat_rem, aes(x=date, y= value), color = "blue") #a bit off at the onset 


# 
 guess <- c( beta_m=0.1, b=0.03) 
# 
# #the parameters are constrained  accordingly (lower and upper)
# 
rem_fit <- optim(fn=func_loglik_2,  par=guess, lower=c(0,0), 
              upper = c( Inf,1), method = "L-BFGS-B", 
          parameters = rem_parameters,dat = dat_rem, hessian=T)
#  #check the values:
    rem_fit
# 
 func_loglik_2(rem_fit$par, dat_rem,rem_parameters) 
# 
 rem_parameters[names(guess)] <- rem_fit$par
#rem_parameters["beta_m"] <- rem_parameters["beta_m"]*(1-rem_parameters[["stngcy"]])



times <- 1:(nrow(dat_rem) + forecasts_days)
rem_out <- as.data.frame(deSolve::ode(y=init_rem,time=times,func= sveirs,
                                          parms=rem_parameters))  
rem_out <- rem_out %>% mutate(incid=last(test_prop)*rem_parameters[["p"]]*
                                        rem_parameters[["sigma"]]*(rem_out$Er + rem_out$Erv + rem_out$Erw + 
                                        rem_out$Em + rem_out$Emv + rem_out$Emw)) %>% 
  mutate(date=seq.Date(ymd(intro_date),ymd(intro_date)-1+length(times), 1))

#with test_prop 





uncert_bound_fit = raply(simu_size,rnbinom(n=length(rem_out$incid),
                                           mu=rem_out$incid,
                                           size=1/rem_parameters[["theta"]]))

project_dat_fit =  as.data.frame(aaply(uncert_bound_fit ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
mutate(date=seq.Date(ymd(intro_date),ymd(intro_date)-1+length(times), 1))



 

ggplot() + geom_line(data=project_dat_BC,aes(x=date,y=`50%`), col="green",size=1.5,alpha=0.4) +
  geom_point(data=dat_reported,aes(x=date, y=value),color='grey48', alpha=0.8, size = 1.5) +
  geom_ribbon(data=project_dat_BC,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='darkgreen',alpha=0.1, size = 1.5) +
  labs(y="Reported cases",x="Date") + ylim(c(0,7000)) + #30000
  geom_line(data=project_dat_fit,aes(x=date,y=`50%`), col="darkgreen",size=1.5,alpha=0.4) +
  geom_point(data=dat_rem,aes(x=date, y=value),color='grey48', alpha=0.8, size = 2) +
  geom_ribbon(data=project_dat_fit,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='darkgreen',alpha=0.1, size = 1.5)+  
  labs(y="Reported cases", x="Date") +
  scale_x_date(date_breaks = "15 days", date_labels = "%b-%d-%y") +theme_light() +
  scale_color_manual(values = cols) +  theme(axis.text=element_text(size=15),
                                             plot.title = element_text(size=15, face="bold"),
                                             axis.text.x = element_text(angle = 45, hjust = 1),
                                             legend.position = "bottom", legend.title = element_text(size=15),
                                             legend.text = element_text(size=15),
                                             axis.title=element_text(size=15,face="bold")) +  
  annotate("rect", xmin = as.Date("2021-11-30"), xmax = as.Date("2022-04-01"), ymin = 0, ymax = Inf,alpha = .2)





############### Scenario projections ##################


#initiate the model with a new mutant strain - e.g. introducing BA.5

# Set the desired characteristics of the new mutant. You can include any of the named elements of rem_parameters to be adjusted. 
# Old mutant will be swapped to resident automatically, and does not need to be defined here
params_newmutant = list("beta_m" = rem_parameters["beta_m"]*1.7, "eff_t" = 600) # eff_t is just pushed back beyond the time horizon of the projection

# Swap resident and mutant, then set up new mutant. 
# This assumes that the new mutant 'arrives' with mut_prop% of current cases
new_model <- swap_strains(out_old = rem_out, params_old = rem_parameters, params_newmutant = params_newmutant, mut_prop = 0.01)
init_proj <- new_model$init_newm
proj_parameters <- new_model$newm_parameters


# Make projections
forecasts_days <- 200
times <- 1:(forecasts_days)
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
  mutate(date=seq.Date(ymd(last(rem_out$date)),ymd(last(rem_out$date))-1+length(times), 1)) 

pivot_longer(proj_out, c(Total, Resident, Mutant), names_to = "Strain", values_to = "count") %>%
ggplot(aes(x=date, y=count, colour=Strain)) + geom_line()






