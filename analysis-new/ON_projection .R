
intro_date <-  ymd("2021-12-07")
stop_date <- ymd("2022-01-20") # make it uniform 
#import data 
dat = readRDS("data/ON-dat.rds")


dat <- dat %>% filter(date >= intro_date & date <= stop_date)
dat_omic <- filter(dat, date <= stop_date)
dat_omic <- filter(dat_omic, date >= intro_date) %>% select(c("day", "value"))
dat_omic$day <- 1:nrow(dat_omic)


tail(dat)
#this assumes we've run catch_data.R (to load testpropdf)
test_prop_ON <- filter(mytest_ON, date >= intro_date)$test_prop


#use a lower test_prop with more flexibility (and smooth)
test_prop_ON <- c(test_prop_ON[1:length(dat_omic$value)], rep(last(test_prop_ON),forecasts_days))
fake_test_prop_ON <- (1 - (1-0.05)/(1 + exp(-0.2*(1:length(test_prop_ON)-32))))

plot(fake_test_prop_ON)
lines(test_prop_ON)

test_prop_ON <- fake_test_prop_ON[1:length(dat_omic$day)]

test_prop <- test_prop_ON 


#set values to generate initial conditions with make_init()
N=14.57e6
vaxlevel_in = 0.78
N_pop = N 
port_wane_in = 0.04 
past_infection_in = 0.28
incres_in = 1300
incmut_in = 150


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
                        N=14.57e6,
                        stngcy= 0,#0.78, #(*%(reduction)) strength of intervention (reduction in beta's)
                        eff_t = as.numeric(eff_date - intro_date),
                        c=1
)

init <- make_init() #generate initial state


#load pomp object for fitting 
source("analysis-new/mod_fitting_setup.R")




#guess <- c(log(1.03), logit(0.65), log(1.34), log(0.01)) #c(log(10),log(15),log(1))
guess <- c(log(1.5), logit(0.51), log(2.6), log(0.1))

#the parameters are constrained  accordingly (lower and upper)

fit_ON <- optim(fn=f_loglik,par=guess, lower=c(log(0.8), 0.2, log(2.2), log(0.15)), 
                upper = c(log(1.7), 0.7, log(2.7), log(0.2)), method = "L-BFGS-B")


#this catches estimated parameter values from MLE 
mle_est_ON <- c(beta_r=exp(fit_ON$par[1]),p=expit(fit_ON$par[2]), 
              beta_m=exp(fit_ON$par[3]),theta=exp(fit_ON$par[4]))

signif(mle_est_ON,3)



#estimated parameter values will now be used to check the fit 

coef(pomp_obj ) <- c(c(S_0=init[[1]],Er_0=init[[2]],Em_0=init[[3]],Ir_0=init[[4]],
                       Im_0=init[[5]],R_0=init[[6]],V_0=init[[7]],Erv_0=init[[8]], 
                       Emv_0=init[[9]],Irv_0=init[[10]],Imv_0=init[[11]],Rv_0=init[[12]],
                       W_0=init[[13]],Erw_0=init[[14]],Emw_0=init[[15]],Irw_0=init[[16]],
                       Imw_0=init[[17]],Rw_0=init[[18]] ,c(parameters,mle_est_ON)))



out_state <-  trajectory(pomp_obj) #Computes trajectories of the deterministic skeleton of a MP.


model.pred <- parameters[[1]]*(out_state["Er",,]+ 
                                 out_state["Erv",,] + out_state["Erw",,]+
                                 out_state["Em",,]+ out_state["Emv",,] + 
                                 out_state["Emw",,])*test_prop[1:length(dat_omic$day)]

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
  geom_point(data=dat_sim,aes(x=date, y=value),color='grey48', alpha=0.8)  + ylim(c(0,10000)) 


########now we make prediction and projections with estimated parameters 

#load initial state    
init_ON <- c(S=init[[1]],Er=init[[2]],Em=init[[3]],Ir=init[[4]],
             Im=init[[5]],R=init[[6]],V=init[[7]],Erv=init[[8]], 
             Emv=init[[9]],Irv=init[[10]],Imv=init[[11]],Rv=init[[12]],
             W=init[[13]],Erw=init[[14]],Emw=init[[15]],Irw=init[[16]],
             Imw=init[[17]],Rw=init[[18]])

eff_date <-   ymd("2021-12-26")  # intervention date 

parameters_ON <- c(sigma=1/3, # incubation period (3 days) (to fixed)
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
                   N=14.57e6,
                   stngcy= 0.48,#0.78, #(*%(reduction)) strength of intervention (reduction in beta's)
                   eff_t = as.numeric(eff_date - intro_date)  # time to 50% intervention effectiveness
)


#use estimated parameters 

parameters <- c(parameters_ON,mle_est_ON)

#use estimated parameters to make projectons 

out_ON <- as.data.frame(deSolve::ode(y=init_ON,time=times,func= sveirs,
                                     parms=parameters))   

#test_prop is shorter than projection, so we'll use the last value of test_prop for the rest of the simulation 




#with test_prop 
incidence_ON =  parameters[[1]]*(out_ON$Er + out_ON$Erv + out_ON$Erw +
                                   out_ON$Em + out_ON$Emv +
                                   out_ON$Emw)*fake_test_prop_ON 

#without test_prop 
incidence_ON_rel =  parameters[[1]]*(out_ON$Er + out_ON$Erv + out_ON$Erw +
                                       out_ON$Em + out_ON$Emv +
                                       out_ON$Emw)

uncert_bound_ON = raply(simu_size,rnbinom(n=length(incidence_ON),
                                          mu=coef(pomp_obj ,"p")*incidence_ON,
                                          size=1/coef(pomp_obj ,"theta")))

quantiles_proj_ON =  as.data.frame(aaply(uncert_bound_ON 
                                         ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975)))

project_dat_ON = quantiles_proj_ON %>% mutate(date=seq.Date(ymd(intro_date),
                                                            ymd(intro_date)-1+length(times), 1))



uncert_bound_ON_rel = raply(simu_size,rnbinom(n=length(incidence_ON_rel),
                                              mu=coef(pomp_obj ,"p")*incidence_ON_rel,
                                              size=1/coef(pomp_obj ,"theta")))

quantiles_proj_ON_rel =  as.data.frame(aaply(uncert_bound_ON_rel,2
                                             ,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975)))

project_dat_ON_rel = quantiles_proj_ON_rel %>% mutate(date=seq.Date(ymd(intro_date),
                                                                    ymd(intro_date)-1+length(times), 1))


saveRDS(project_dat_ON, file.path("data/ON_test_constraints.rds"))
project_dat_ON =readRDS(file.path("data/ON_test_constraints.rds"))

saveRDS(project_dat_ON_rel, file.path("data/ON_no_constraints.rds"))
project_dat_ON_rel =readRDS(file.path("data/ON_no_constraints.rds"))



#make figures 

cols <- c("Current, with testing constraints" = "darkgreen", "Without testing constraints"="orange")

gg_ON <- ggplot() + geom_line(data=project_dat_ON,aes(x=date,y=`50%`, colour = "Current, with testing constraints"),size=1.2,alpha=0.4) +
  geom_ribbon(data=project_dat_ON,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='darkgreen',alpha=0.1)+
  geom_line(data=project_dat_ON_rel,aes(x=date,y=`50%`, color="Without testing constraints"),size=1.2,alpha=0.4) +
  geom_ribbon(data=project_dat_ON_rel,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='orange',alpha=0.1)+
  geom_point(data=dat_sim,aes(x=date, y=value),color='grey48', alpha=0.8) + 
  #geom_line(aes(y=typical),color='blue') +
  labs(y="Reported cases",x="Date") + ylim(c(0,140000)) + 
  scale_x_date(date_breaks = "8 days", date_labels = "%b-%d-%y") +theme_light() +
  scale_color_manual(values = cols) +  theme(axis.text=element_text(size=15),
                                             plot.title = element_text(size=15, face="bold"),
                                             axis.text.x = element_text(angle = 45, hjust = 1),
                                             legend.position = "bottom", legend.title = element_text(size=15),
                                             legend.text = element_text(size=15),
                                             axis.title=element_text(size=15,face="bold")) +
  
  labs(color = " ",title="ON") +  geom_vline(xintercept=eff_date, linetype="dashed", 
                                             color = "grey", size=1)


gg_ON

ggsave(file="figs/ON_proj.png", gg_ON, width = 10, height = 8)
saveRDS(gg_ON, file.path("figs/ON-fig.rds"))

