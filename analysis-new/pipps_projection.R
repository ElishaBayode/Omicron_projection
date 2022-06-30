#add dat to data for plotting 
#here I match model output to the remaining data point (i.e. from March 31, 2022 onward)
dat_full = readRDS("data/BC-dat.rds")

forecasts_days <- 1 
intro_date <-  ymd("2022-03-31")
stop_date <- last(dat_full$date)

dat_rem <- dat_full %>% filter(date >= intro_date &  date <= stop_date)

#dat_rem <- filter(dat_full, date >= intro_date) %>% select(c("day", "value", "date"))
dat_full$day <- 1:nrow(dat_full)

plot(dat_rem$value)

#initiating with the last point on fit from pipps_simulation.R script 
init_rem <-  c(S=last(out_samp$S),
                Er=last(out_samp$Er),Em=last(out_samp$Em),
                Ir=last(out_samp$Ir),Im=last(out_samp$Im),
                R=last(out_samp$R),V=last(out_samp$V),
                Erv=last(out_samp$Erv),Emv=last(out_samp$Emv),
                Irv=last(out_samp$Irv),Imv=last(out_samp$Imv),
                Rv=last(out_samp$Rv),W=last(out_samp$W),Erw=last(out_samp$Erw)
                ,Emw=last(out_samp$Emw),Irw=last(out_samp$Irw),Imw=last(out_samp$Imw)
                ,Rw=last(out_samp$Rw)) 

rem_parameters  <- parameters

rem_parameters["beta_r"] <- rem_parameters["beta_r"]*(1-rem_parameters[["stngcy"]])
rem_parameters["beta_m"] <- rem_parameters["beta_m"]*(1-rem_parameters[["stngcy"]])*1.5
rem_parameters["eff_t"]  <- 100
rem_parameters[["stngcy"]] <- 0.35
#rem_parameters[["p"]] <- 0.5
rem_parameters[["beff"]] <- 0.5

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



guess <- c( beta_m=0.27) 

#the parameters are constrained  accordingly (lower and upper)

rem_fit <- optim(fn=func_loglik_2,  par=guess, lower=c(0), 
             upper = c(Inf), method = "L-BFGS-B", 
           parameters = rem_parameters,dat = dat_rem, hessian=T)
 #check the values:
   rem_fit

func_loglik_2(rem_fit$par, dat_rem,rem_parameters) 

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


#initiate the model with a new mutant starin 
#switch r for m ########### 

init_proj <-  c(S=last(rem_out$S),Er=5, Em=last(rem_out$Er) + (last(rem_out$Em)-5),
               Ir=4,Im=last(rem_out$Ir) + (last(rem_out$Im)-4), R=last(rem_out$R),
               V=last(rem_out$V),
               Erv=5,   Emv= last(rem_out$Erv) + (last(rem_out$Emv)-5),
               Irv=1 ,Imv= last(rem_out$Irv) + (last(rem_out$Imv)-1),
               Rv=last(rem_out$Rv),
               W=last(rem_out$W),  
               Erw=1,Emw=last(rem_out$Erw) + (last(rem_out$Emw)-1),
               Irw=5,Imw= last(rem_out$Irw) + (last(rem_out$Imw)- 5)
               ,Rw=last(rem_out$Rw)) 


proj_parameters <- rem_parameters

proj_parameters["beta_r"] <- proj_parameters["beta_m"]*3
proj_parameters["eff_t"]  <- 600
proj_parameters[["stngcy"]] <- 0.35

forecasts_days <- 100

times <- 1:(forecasts_days)
proj_out <- as.data.frame(deSolve::ode(y=init_proj, time=times,func= sveirs,
                                      parms=proj_parameters))  
proj_out <- proj_out %>% mutate(incid=last(test_prop)*proj_parameters[["p"]]*
                              proj_parameters[["sigma"]]*(proj_out$Er + proj_out$Erv + proj_out$Erw + 
                              proj_out$Em + proj_out$Emv + proj_out$Emw)) %>% 
  mutate(date=seq.Date(ymd(last(rem_out$date)),ymd(last(rem_out$date))-1+length(times), 1))

plot(proj_out$incid)




