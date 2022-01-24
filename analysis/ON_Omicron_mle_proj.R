#including Omicron wave only 
intro_date <-  ymd("2021-12-07")
restart_date <- ymd("2022-01-01")
#importing data 
dat = readRDS("data/ON-dat.rds")

dat <- dat %>% filter(date >= intro_date)




dat_omic <- filter(dat, date <= restart_date)
dat_omic <- filter(dat_omic, date >= intro_date) %>% select(c("day", "value"))
dat_omic$day <- 1:nrow(dat_omic)


test_prop_ON <- filter(mytest_ON, date >= intro_date)$test_prop
test_prop <- test_prop_ON

source("analysis/mod_fitting_setup.R")
#setting values to generate initial conditions with make_init()
N=14.57e6
vaxlevel_in = 0.78
N_pop = N 
port_wane_in = 0.04 
past_infection_in = 0.28
incres_in = 1300
incmut_in = 60

intro_date <- ymd("2021-12-07")


#to generatee initial state 



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
                        c=1
)


init <- make_init()





#guess <- c(log(1.03), logit(0.65), log(1.34), log(0.01)) #c(log(10),log(15),log(1))
guess <- c(log(1.5), logit(0.51), log(2.6), log(0.1))

#the parameters are constrained  accordingly (lower and upper)

fit_ON <- optim(fn=f_loglik,par=guess, lower=c(log(0.8), 0.2, log(2.2), log(0.15)), 
                upper = c(log(1.7), 0.7, log(2.7), log(0.2)), method = "L-BFGS-B")


#this catches estimated parameter values from MLE 
mle_est_ON <- c(beta_r=exp(fit_ON$par[1]),p=expit(fit_ON$par[2]), beta_m=exp(fit_ON$par[3]),theta=exp(fit_ON$par[4]))

signif(mle_est_ON,3)


#estimated parameter values will now be used for for prediction 
coef(pomp_obj ) <- c(c(S_0=init[[1]],Er_0=init[[2]],Em_0=init[[3]],Ir_0=init[[4]],
                       Im_0=init[[5]],R_0=init[[6]],V_0=init[[7]],Erv_0=init[[8]], 
                       Emv_0=init[[9]],Irv_0=init[[10]],Imv_0=init[[11]],Rv_0=init[[12]],
                       W_0=init[[13]],Erw_0=init[[14]],Emw_0=init[[15]],Irw_0=init[[16]],
                       Imw_0=init[[17]],Rw_0=init[[18]] ,c(parameters,mle_est_ON)))






forecast_days=30 # This should match forecasts_days in case_projection()
eff_date <-   ymd("2022-01-01") - forecast_days
parameters_1 <- c(sigma=1/3, # incubation period (3 days) (to fixed)
                  gamma=1/(5), #recovery rate (fixed)
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
                  stngcy= 0,#0.78, #(*%(reduction)) strength of intervention (reduction in beta's)
                  eff_t = as.numeric(eff_date - intro_date)  # time to 50% intervention effectiveness
) 



# note 19% boosted as of Dec

estim_parameters <- c(parameters_1,mle_est_ON)


#an attempt to model changes in ascertainment probability 
#test_prop <- 1 - 0.55/(1+ exp(-1.25*(1:nrow(dat_omic)-20)))
#test_prop1 <- 1 - 0.75/(1+ exp(-1.25*(1:nrow(dat_omic)-100)))
#test_prop2 <- 1 - 0.65/(1+ exp(-1.25*(1:nrow(dat_omic)-100)))
#test_prop <- test_prop*test_prop1*test_prop2

# model predictions (incidence, sigma*(Er +Em + Erv + Emv + Erw + Emw)) 

out_state <-  trajectory(pomp_obj )

model.pred <- parameters[[1]]*(out_state["Er",,]+ 
                                 out_state["Erv",,] + out_state["Erw",,]+
                                 out_state["Em",,]+ out_state["Emv",,] + 
                                 out_state["Emw",,])



##################--------------------fitting the remaining data points 


dat_omic_end <- filter(dat, date > restart_date)
dat_omic  <- dat_omic_end 
dat_omic <- dat_omic_end  %>% select(c("day", "value"))
dat_omic$day <- 1:nrow(dat_omic)
source("analysis/mod_fitting_setup.R")


guess_part2 <- c(log(0.2), logit(0.51), log(0.8),log(0.01))

#the parameters are constrained  accordingly (lower and upper)

fit_ON2 <- optim(fn=f_loglik_2,par=guess_part2, lower=c(log(0.5), 0.2, log(0.8),  log(0.1)), 
                 upper = c(log(1.2), 0.7, log(5), log(0.2)), method = "L-BFGS-B")



#this catches estimated parameter values from MLE 
mle_est_ON2 <- c(beta_r=exp(fit_ON2$par[1]),p=expit(fit_ON$par[2]), beta_m=exp(fit_ON2$par[3]),theta=exp(fit_ON$par[4]))

signif(mle_est_ON2,3)

coef(pomp_obj) <- c(c(S_0=last(out_state["S",,]),Er_0=last(out_state["Er",,]),Em_0=last(out_state["Em",,]),
                      Ir_0=last(out_state["Ir",,]),Im_0=last(out_state["Im",,]), R_0=last(out_state["R",,]),
                      V_0=last(out_state["V",,]),Erv_0=last(out_state["Erv",,]), Emv_0=last(out_state["Emv",,]), 
                      Irv_0=last(out_state["Irv",,]),Imv_0=last(out_state["Imv",,]), Rv_0=last(out_state["Rv",,]),
                      W_0=last(out_state["W",,]),Erw_0=last(out_state["Erw",,]), Emw_0=last(out_state["Emw",,]),
                      Irw_0=last(out_state["Irw",,]),Imw_0=last(out_state["Imw",,]), Rw_0=last(out_state["Rw",,]),
                      c(parameters,mle_est_ON2)))


#


estim_parameters <- c(parameters_1,mle_est_ON2)


out_state_2 <-  trajectory(pomp_obj)

model.pred_2 <-   parameters[[1]]*(out_state_2["Er",,]+ 
                                     out_state_2["Erv",,] + out_state_2["Erw",,]+
                                     out_state_2["Em",,]+ out_state_2["Emv",,] + 
                                     out_state_2["Emw",,])

model.pred_fake  <- model.pred
model.pred_rel <- c(model.pred_fake,model.pred_2)
model.pred <- c(model.pred_fake,model.pred_2)*test_prop[1:length(model.pred_rel)]







# (out_state["Ir",,]+ 
#                 out_state["Irv",,] + out_state["Irw",,]+
#                out_state["Im",,]+ out_state["Imv",,] + 
#                 out_state["Imw",,]) + (out_state["Er",,]+ 
#                                         out_state["Erv",,] + out_state["Erw",,]+
#                                         out_state["Em",,]+ out_state["Emv",,] + 
#                                        out_state["Emw",,])

#raply(10000,rnbinom(n=length(pred_incid),
#                    mu=pred_incid,
#                   size=1/coef(pomp_obj ,"theta"))) -> sim_true_inc

#aaply(sim_true_inc,2,quantile,probs=c(0.025,0.5,0.975)) -> quantiles_true 
#this generates multiple  model realizations 

raply(1000000,rnbinom(n=length(model.pred),
                      mu=coef(pomp_obj ,"p")*model.pred,
                      size=1/coef(pomp_obj ,"theta"))) -> simdat

raply(100000,rnbinom(n=length(model.pred_rel),
                     mu=coef(pomp_obj ,"p")*model.pred_rel,
                     size=1/coef(pomp_obj ,"theta"))) -> simdat_rel

#this generates appriopriate quantiles

aaply(simdat,2,quantile,probs=c(0.025,0.5,0.975)) -> quantiles

aaply(simdat_rel,2,quantile,probs=c(0.025,0.5,0.975)) -> quantiles_rel
#sample a typical model realization (sanity check)
#typ <- sample(nrow(simdat),1) 

dat$day <- 1:nrow(dat)
dat_sim <- cbind(select(dat, c("day", "value")),
                 quantiles) 
dat_sim_rel <- cbind(select(dat, c("day", "value")),
                     quantiles_rel) 

#true_inc <-  cbind(as.data.frame(pomp_obj ),
#                  quantiles_true) 

#create dates 
dat_sim = dat_sim %>% mutate(date=seq.Date(ymd(intro_date),
                                           ymd(intro_date)-1+nrow(dat), 1))

dat_sim_rel = dat_sim_rel %>% mutate(date=seq.Date(ymd(intro_date),
                                                   ymd(intro_date)-1+nrow(dat), 1))



#true_inc  = true_inc  %>% mutate(date=seq.Date(ymd(intro_date),
#                                           ymd(intro_date)-1+nrow(dat_omic), 1))




######################

#### get uncorrected cases directly ##############3

# init_dr <- c(c(S=init[[1]],Er=init[[2]],Em=init[[3]],Ir=init[[4]],
#                Im=init[[5]],R=init[[6]],V=init[[7]],Erv=init[[8]], 
#                Emv=init[[9]],Irv=init[[10]],Imv=init[[11]],Rv=init[[12]],
#                W=init[[13]],Erw=init[[14]],Emw=init[[15]],Irw=init[[16]],
#                Imw=init[[17]],Rw=init[[18]]))
# 
# output_dr = as.data.frame(deSolve::ode(y=init_dr,times=1:nrow(dat),func=sveirs,
#                                        parms=c(parameters_1,mle_est_ON)))       
# 
# incidence_dr =  (estim_parameters[[1]]*(output_dr$Er + output_dr$Erv+ output_dr$Erw +
#                                          output_dr$Em + output_dr$Emv + output_dr$Emw))
# 
# 
# raply(100000,rnbinom(n=length(incidence_dr),
#                         mu=coef(pomp_obj ,"p")*incidence_dr,
#                         size=1/coef(pomp_obj ,"theta"))) -> simdat_dr
# 
# 
# 
# aaply(simdat_dr,2,quantile,probs=c(0.025,0.5,0.975)) -> quantiles_dr
# 
# dat_sim_dr <- as.data.frame( cbind(day=1:(nrow(dat)),
#                                    quantiles_dr) )
# 
# dat_sim_dr = dat_sim_dr %>% mutate(date=seq.Date(ymd(intro_date),
#                  ymd(intro_date)-1+nrow(dat), 1))

tail(dat_sim)
tail(dat_sim_rel)


######## forecast cases ############### 
ON_forecast <- case_projection(out_state=out_state_2,
                               forecast_days=30, 
                               test_prop = test_prop[1:length(model.pred_rel)],
                               lag=0,
                               parameters_estim = estim_parameters, 
                               simu_size = 100000,
                               pomp_obj = pomp_obj , 
                               dat_sim = dat_sim
                               
)



parameters_2 <- c(sigma=1/3, # incubation period (3 days) (to fixed)
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
                  stngcy= 0,#0.78, #(*%(reduction)) strength of intervention (reduction in beta's)
                  eff_t = as.numeric(eff_date - intro_date)  # time to 50% intervention effectiveness
)

estim_parameters <- c(parameters_2,mle_est_ON2)

ON_forecast_int <- case_projection_rel(out_state=out_state,
                                       forecast_days=30, 
                                       lag=0,
                                       parameters_estim = estim_parameters, 
                                       simu_size = 100000,
                                       pomp_obj = pomp_obj , 
                                       dat_sim = dat_sim_rel
                                       
)


saveRDS(dat_sim, file.path("data/ON_dat_sim.rds"))
dat_sim =readRDS(file.path("data/ON_dat_sim.rds"))

saveRDS(dat_sim_rel, file.path("data/ON_dat_sim_rel.rds"))
dat_sim_rel =readRDS(file.path("data/ON_dat_sim_rel.rds"))

saveRDS(ON_forecast, file.path("data/ON_forecast.rds"))
ON_forecast =readRDS(file.path("data/ON_forecast.rds"))

saveRDS(ON_forecast_int, file.path("data/ON_forecast_int.rds"))
ON_forecast_int =readRDS(file.path("data/ON_forecast_int.rds"))


cols <- c("Current, with testing constraints" = "darkgreen", "Without testing constraints"="orange")

gg_ON <- ggplot() + geom_line(data=dat_sim,aes(x=date,y=`50%`),color='blue',size=1.2,alpha=0.4) +
  geom_ribbon(data=dat_sim,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='blue',alpha=0.1)+
  geom_line(data=dat_sim_rel,aes(x=date,y=`50%`),color='orange',size=1.2,alpha=0.4) +
  geom_ribbon(data=dat_sim_rel,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='orange',alpha=0.1)+
  geom_point(data=dat_sim,aes(x=date, y=value),color='grey28', alpha=0.5) +
  geom_line(data=ON_forecast,aes(x=date, y=`50%`,color="Current, with testing constraints"),size=1.2,alpha=0.5) +
  geom_ribbon(data=ON_forecast,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='darkgreen',alpha=0.1)+
  geom_line(data=ON_forecast_int,aes(x=date, y=`50%`, color="Without testing constraints"),size=1.2,alpha=0.5) +
  geom_ribbon(data=ON_forecast_int,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='orange',alpha=0.1)+
  #geom_line(aes(y=typical),color='blue') +
  labs(y="Reported cases",x="Date") + ylim(c(0,120000)) + 
  scale_x_date(date_breaks = "8 days", date_labels = "%b-%d-%y") +theme_light() +
  scale_color_manual(values = cols) +  theme(axis.text=element_text(size=15),
                                             plot.title = element_text(size=15, face="bold"),
                                             axis.text.x = element_text(angle = 45, hjust = 1),
                                             legend.position = "bottom", legend.title = element_text(size=15),
                                             legend.text = element_text(size=15),
                                             axis.title=element_text(size=15,face="bold")) +
  
  labs(color = " ",title="ON")


gg_ON

ggsave(file="figs/ON_proj.png", gg_ON, width = 10, height = 8)
saveRDS(gg_ON, file.path("figs/ON-fig.rds"))