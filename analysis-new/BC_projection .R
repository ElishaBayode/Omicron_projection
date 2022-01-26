#include Omicron wave only 

intro_date <-  ymd("2021-12-07")
stop_date <- ymd("2022-01-23") # make it uniform 
#import data 
dat = readRDS("data/BC-dat.rds")


dat <- dat %>% filter(date >= intro_date)
dat_omic <- filter(dat, date <= stop_date)
dat_omic <- filter(dat_omic, date >= intro_date) %>% select(c("day", "value"))
dat_omic$day <- 1:nrow(dat_omic)


#this assumes we've run catch_data.R (to load testpropdf)
test_prop_BC <- filter(testpropdf, date >= intro_date)$test_prop
test_prop <- test_prop_BC


#set values to generate initial conditions with make_init()

N=5.07e6
N_pop=N
#ascFrac <- 0.5
vaxlevel_in = 0.88
port_wane_in = 0.04 
past_infection_in = 0.2  #increase this from 0.1 to 0.2
incres_in = 300
incmut_in = 70
simu_size = 1e5
forecasts_days =30
times = 1:(nrow(dat_omic) + forecasts_days)
 



#declaring fixed parameters 

parameters <-         c(sigma=1/3, # incubation period (3 days) (to fixed)
                        gamma=1/(4), #recovery rate (fixed)
                        nu =0.007, #vax rate: 0.7% per day (fixed)
                        mu=1/(82*365), # 1/life expectancy (fixed)
                        w1= 1/(3*365),# waning rate from R to S (fixed)
                        w2= 1/(3*365), # waning rate from Rv to V (fixed)
                        w3= 1/(3*365),# waning rate Rw to W (fixed)
                        ve=1, # I think this should be 1. it is not really efficacy  ( fixed)
                        #beta_r=0.72, #transmission rate (to estimate) (0.35)
                        #beta_m=0.8*2.2, #transmission rate (to estimate)(*1.9)
                        epsilon_r = (1-0.8), # % this should be 1-ve 
                        epsilon_m = (1-0.6), # % escape capacity #(fixed)
                        b= 0.006, # booster rate  (fixed)
                        beff = 0.7, # booster efficacy
                        wf=0.2, # protection for newly recoverd #0.2
                        N=5e6,
                        stngcy= 0,#0.78, #(*%(reduction)) strength of intervention (reduction in beta's)
                        eff_t = as.numeric(eff_date - intro_date),
                        c=1
)

init <- make_init() #generate initial state


#load pomp object for fitting 
source("analysis-new/mod_fitting_setup.R")



#guess <- c(log(1.03), logit(0.65), log(1.34), log(0.01)) #c(log(10),log(15),log(1))
guess <- c(log(0.81), logit(0.8), log(2.1), log(0.01))

#the parameters are constrained  accordingly (lower and upper)

fit_BC <- optim(fn=f_loglik,par=guess, lower=c(log(0.6), 0.0001, log(2), log(0.1)), 
                upper = c(log(0.8), 0.001, log(2.5), log(0)), method = "L-BFGS-B")


#this catches estimated parameter values from MLE 
mle_est_BC <- c(beta_r=exp(fit_BC$par[1]),p=expit(fit_BC$par[2]), 
                beta_m=exp(fit_BC$par[3]),theta=exp(fit_BC$par[4]))

signif(mle_est_BC,3)




#estimated parameter values will now be used to check the fit 

coef(pomp_obj ) <- c(c(S_0=init[[1]],Er_0=init[[2]],Em_0=init[[3]],Ir_0=init[[4]],
                       Im_0=init[[5]],R_0=init[[6]],V_0=init[[7]],Erv_0=init[[8]], 
                       Emv_0=init[[9]],Irv_0=init[[10]],Imv_0=init[[11]],Rv_0=init[[12]],
                       W_0=init[[13]],Erw_0=init[[14]],Emw_0=init[[15]],Irw_0=init[[16]],
                       Imw_0=init[[17]],Rw_0=init[[18]] ,c(parameters,mle_est_BC)))



out_state <-  trajectory(pomp_obj) #Computes trajectories of the deterministic skeleton of a MP.


model.pred <- parameters[[1]]*(out_state["Er",,]+ 
                                 out_state["Erv",,] + out_state["Erw",,]+
                                 out_state["Em",,]+ out_state["Emv",,] + 
                                 out_state["Emw",,])*test_prop

raply(simu_size,rnbinom(n=length(model.pred),
                      mu=coef(pomp_obj ,"p")*model.pred,
                      size=1/coef(pomp_obj ,"theta"))) -> simdat



aaply(simdat,2,quantile,probs=c(0.025,0.5,0.975)) -> quantiles
dat$day <- 1:nrow(dat)
dat_sim <- cbind(select(dat, c("day", "value")),
                 quantiles) 
dat_sim = dat_sim %>% mutate(date=seq.Date(ymd(intro_date),
                                           ymd(intro_date)-1+nrow(dat), 1))


#check fit, initial rises need to align   

ggplot() + geom_line(data=dat_sim,aes(x=date,y=`50%`),color='blue',size=1.2,alpha=0.4) +
geom_ribbon(data=dat_sim,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='blue',alpha=0.1) +
geom_point(data=dat_sim,aes(x=date, y=value),color='grey48', alpha=0.8)  + ylim(c(0,4000)) 


########now we make prediction and projections with estimated parameters 

#load initial state    
init_BC <- c(S=init[[1]],Er=init[[2]],Em=init[[3]],Ir=init[[4]],
           Im=init[[5]],R=init[[6]],V=init[[7]],Erv=init[[8]], 
           Emv=init[[9]],Irv=init[[10]],Imv=init[[11]],Rv=init[[12]],
           W=init[[13]],Erw=init[[14]],Emw=init[[15]],Irw=init[[16]],
           Imw=init[[17]],Rw=init[[18]])

eff_date <-   ymd("2021-12-30")  # intervention date 

parameters_BC <- c(sigma=1/3, # incubation period (3 days) (to fixed)
                  gamma=1/(4), #recovery rate (fixed)
                  nu =0.007, #vax rate: 0.7% per day (fixed)
                  mu=1/(82*365), # 1/life expectancy (fixed)
                  w1= 1/(3*365),# waning rate from R to S (fixed)
                  w2= 1/(3*365), # waning rate from Rv to V (fixed)
                  w3= 1/(3*365),# waning rate Rw to W (fixed)
                  ve=1, # I think this should be 1. it is not really efficacy  ( fixed)
                  epsilon_r = (1-0.8), # % this should be 1-ve 
                  epsilon_m = (1-0.6), # % escape capacity #(fixed)
                  b= 0.006, # booster rate  (fixed) orig 0.006 
                  beff = 0.7, # booster efficacy
                  wf=0.2, # protection for newly recoverd #0.2
                  N=5e6,
                  stngcy= 0.5,#0.78, #(*%(reduction)) strength of intervention (reduction in beta's)
                  eff_t = as.numeric(eff_date - intro_date)  # time to 50% intervention effectiveness
)


#use estimated parameters 

parameters <- c(parameters_BC,mle_est_BC)

#use estimated parameters to make projectons 

out_BC <- as.data.frame(deSolve::ode(y=init_BC,time=times,func= sveirs,
                           parms=parameters))   

#test_prop is shorter than projection, so we'll use the last value of test_prop for the rest of the simulation 
test_prop_BC <- c(test_prop, rep(last(test_prop),forecasts_days))

#with test_prop 
incidence_BC =  parameters[[1]]*(out_BC$Er + out_BC$Erv + out_BC$Erw +
                                                       out_BC$Em + out_BC$Emv +
                                                       out_BC$Emw)*test_prop_BC 

#without test_prop 
incidence_BC_rel =  parameters[[1]]*(out_BC$Er + out_BC$Erv + out_BC$Erw +
                                   out_BC$Em + out_BC$Emv +
                                   out_BC$Emw)

uncert_bound_BC = raply(simu_size,rnbinom(n=length(incidence_BC),
                                       mu=coef(pomp_obj ,"p")*incidence_BC,
                                       size=1/coef(pomp_obj ,"theta")))

quantiles_proj_BC =  as.data.frame(aaply(uncert_bound_BC 
                      ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975)))

project_dat_BC = quantiles_proj_BC %>% mutate(date=seq.Date(ymd(intro_date),
                                       ymd(intro_date)-1+length(times), 1))



uncert_bound_BC_rel = raply(simu_size,rnbinom(n=length(incidence_BC_rel),
                                  mu=coef(pomp_obj ,"p")*incidence_BC_rel,
                                  size=1/coef(pomp_obj ,"theta")))

quantiles_proj_BC_rel =  as.data.frame(aaply(uncert_bound_BC_rel,2
                        ,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975)))

project_dat_BC_rel = quantiles_proj_BC_rel %>% mutate(date=seq.Date(ymd(intro_date),
                                              ymd(intro_date)-1+length(times), 1))


saveRDS(project_dat_BC, file.path("data/test_constraints.rds"))
project_dat_BC =readRDS(file.path("data/test_constraints.rds"))

saveRDS(project_dat_BC_rel, file.path("data/no_constraints.rds"))
project_dat_BC_rel =readRDS(file.path("data/no_constraints.rds"))



#make figures 

cols <- c("Current, with testing constraints" = "darkgreen", "Without testing constraints"="orange")

gg_BC <- ggplot() + geom_line(data=project_dat_BC,aes(x=date,y=`50%`, colour = "Current, with testing constraints"),size=1.2,alpha=0.4) +
  geom_ribbon(data=project_dat_BC,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='darkgreen',alpha=0.1)+
  geom_line(data=project_dat_BC_rel,aes(x=date,y=`50%`, color="Without testing constraints"),size=1.2,alpha=0.4) +
  geom_ribbon(data=project_dat_BC_rel,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='orange',alpha=0.1)+
  geom_point(data=dat_sim,aes(x=date, y=value),color='grey48', alpha=0.8) + 
  #geom_line(aes(y=typical),color='blue') +
  labs(y="Reported cases",x="Date") + ylim(c(0,30000)) + 
  scale_x_date(date_breaks = "8 days", date_labels = "%b-%d-%y") +theme_light() +
  scale_color_manual(values = cols) +  theme(axis.text=element_text(size=15),
                                             plot.title = element_text(size=15, face="bold"),
                                             axis.text.x = element_text(angle = 45, hjust = 1),
                                             legend.position = "bottom", legend.title = element_text(size=15),
                                             legend.text = element_text(size=15),
                                             axis.title=element_text(size=15,face="bold")) +
  
  labs(color = " ",title="BC") +  geom_vline(xintercept=eff_date, linetype="dashed", 
                                              color = "grey", size=1)


gg_BC

ggsave(file="figs/BC_proj.png", gg_BC, width = 10, height = 8)
saveRDS(gg_BC, file.path("figs/BC-fig.rds"))






#plot(NA,NA, xlim=c(0,78), ylim=c(0,3e4))
#lines(dat_sim$value)
#lines(dat_sim$`50%`)
#lines(project_dat_BC$`50%`)
#lines(project_dat_BC$`50%`/c(test_prop, rep(last(test_prop),forecasts_days)))
#dat_sim$`50%` - project_dat_BC$`50%`









  






















