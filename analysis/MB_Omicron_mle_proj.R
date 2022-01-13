
dat = readRDS("data/MB-dat.rds")

#including Omicron wave only 
intro_date <- ymd("2021-12-07")
dat_omic <- filter(dat, date <= ymd("2022-01-09"))
dat_omic <- filter(dat_omic, date >= intro_date) %>% select(c("day", "value"))
dat_omic$day <- 1:nrow(dat_omic)
source("analysis/mod_fitting_setup.R")
#setting values to generate initial conditions with make_init()
N = 1.36e6
N_pop = N
vaxlevel_in = 0.78
port_wane_in = 0.1
past_infection_in = 0.1
incres_in = 200
incmut_in = 5 

parameters <-         c(sigma=1/3, # incubation period (3 days) (to fixed)
                        gamma=1/(5), #recovery rate (fixed)
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
                        N=1.36e6,
                        c=1
)

init <- make_init() #to generatee initial state 
#importing data 



#declaring fixed parameters 




#guess <- c(log(1.03), logit(0.65), log(1.34), log(0.01)) #c(log(10),log(15),log(1))
guess <- c(log(0.8), logit(0.52), log(2.1), log(0.15))

#the parameters are constrained  accordingly (lower and upper)

fit_MB <- optim(fn=f_loglik,par=guess, lower=c(log(0.3), 0.1, log(1.7), log(0.1)), 
                upper = c(log(0.0), 0.6, log(2.9), log(0.3)), method = "L-BFGS-B")


#this catches estimated parameter values from MLE 
mle_est_MB <- c(beta_r=exp(fit_MB$par[1]),p=expit(fit_MB$par[2]), beta_m=exp(fit_MB$par[3]),theta=exp(fit_MB$par[4]))

signif(mle_est_MB,3)


#estimated parameter values will now be used for for prediction 
coef(pomp_obj ) <- c(c(S_0=init[[1]],Er_0=init[[2]],Em_0=init[[3]],Ir_0=init[[4]],
                       Im_0=init[[5]],R_0=init[[6]],V_0=init[[7]],Erv_0=init[[8]], 
                       Emv_0=init[[9]],Irv_0=init[[10]],Imv_0=init[[11]],Rv_0=init[[12]],
                       W_0=init[[13]],Erw_0=init[[14]],Emw_0=init[[15]],Irw_0=init[[16]],
                       Imw_0=init[[17]],Rw_0=init[[18]] ,c(parameters,mle_est_MB)))

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






estim_parameters <- c(parameters_1,mle_est_MB)


#an attempt to model changes in ascertainment probability 
test_prop <- 1 - 0.4/(1+ exp(-1.25*(1:nrow(dat_omic)-19)))
test_prop1 <- 1 - 0.75/(1+ exp(-1.25*(1:nrow(dat_omic)-100)))
test_prop2 <- 1 - 0.65/(1+ exp(-1.25*(1:nrow(dat_omic)-100)))
test_prop <- test_prop*test_prop1*test_prop2

# model predictions (incidence, sigma*(Er +Em + Erv + Emv + Erw + Emw)) 

out_state <-  trajectory(pomp_obj )

model.pred <- test_prop*parameters[[1]]*(out_state["Er",,]+ 
                                           out_state["Erv",,] + out_state["Erw",,]+
                                           out_state["Em",,]+ out_state["Emv",,] + 
                                           out_state["Emw",,])

# pred_incid <- (out_state["Ir",,]+ 
#                 out_state["Irv",,] + out_state["Irw",,]+
#                out_state["Im",,]+ out_state["Imv",,] + 
#                 out_state["Imw",,]) + (out_state["Er",,]+ 
#                                         out_state["Erv",,] + out_state["Erw",,]+
#                                         out_state["Em",,]+ out_state["Emv",,] + 
#                                        out_state["Emw",,])

#raply(10000,rnbinom(n=length(pred_incid),
 #                   mu=pred_incid,
#                    size=1/coef(pomp_obj ,"theta"))) -> sim_true_inc

#aaply(sim_true_inc,2,quantile,probs=c(0.025,0.5,0.975)) -> quantiles_true 
#this generates multiple  model realizations 

raply(10000,rnbinom(n=length(model.pred),
                    mu=coef(pomp_obj ,"p")*model.pred,
                    size=1/coef(pomp_obj ,"theta"))) -> simdat

#this generates appriopriate quantiles

aaply(simdat,2,quantile,probs=c(0.025,0.5,0.975)) -> quantiles


#sample a typical model realization (sanity check)
typ <- sample(nrow(simdat),1) 


dat_sim <- cbind(as.data.frame(pomp_obj ),
                 quantiles,
                 typical=simdat[typ,]) 

#true_inc <-  cbind(as.data.frame(pomp_obj ),
#                  quantiles_true) 

#create dates 
dat_sim = dat_sim %>% mutate(date=seq.Date(ymd(intro_date),
                                           ymd(intro_date)-1+nrow(dat_omic), 1))


#true_inc  = true_inc  %>% mutate(date=seq.Date(ymd(intro_date),
#                                           ymd(intro_date)-1+nrow(dat_omic), 1))

tail(dat_sim)



######## forecast cases ############### 
MB_forecast <- case_projection(out_state=out_state,
                               forecast_days=30, 
                               lag=0,
                               parameters_estim = estim_parameters, 
                               simu_size = 100000,
                               pomp_obj = pomp_obj , 
                               dat_sim = dat_sim
                               
)



parameters_2 <- c(sigma=1/3, # incubation period (3 days) (to fixed)
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
                  stngcy= 0.5,#0.78, #(*%(reduction)) strength of intervention (reduction in beta's)
                  eff_t = as.numeric(eff_date - intro_date)  # time to 50% intervention effectiveness
)

estim_parameters <- c(parameters_2,mle_est_MB)

MB_forecast_int <- case_projection(out_state=out_state,
                                   forecast_days=30, 
                                   lag=0,
                                   parameters_estim = estim_parameters, 
                                   simu_size = 100000,
                                   pomp_obj = pomp_obj , 
                                   dat_sim = dat_sim
                                   
)


saveRDS(dat_sim, file.path("data/MB_dat_sim.rds"))
dat_sim =readRDS(file.path("data/MB_dat_sim.rds"))



saveRDS(MB_forecast, file.path("data/MB_forecast.rds"))
MB_forecast =readRDS(file.path("data/MB_forecast.rds"))

saveRDS(MB_forecast_int, file.path("data/MB_forecast_int.rds"))
MB_forecast_int =readRDS(file.path("data/MB_forecast_int.rds"))


cols <- c("Current" = "orange", "50% reduction in transmission"="green")

gg_MB <- ggplot() + geom_line(data=dat_sim,aes(x=date,y=`50%`),color='blue',size=1.2,alpha=0.4) +
  geom_ribbon(data=dat_sim,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='blue',alpha=0.1)+
  geom_point(data=dat_sim,aes(x=date, y=value),color='grey28', alpha=0.5) +
  geom_line(data=MB_forecast,aes(x=date, y=`50%`,color="Current"),size=1.2,alpha=0.5) +
  geom_ribbon(data=MB_forecast,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='orange',alpha=0.1)+
  geom_line(data=MB_forecast_int,aes(x=date, y=`50%`, color="50% reduction in transmission"),size=1.2,alpha=0.5) +
  geom_ribbon(data=MB_forecast_int,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='green',alpha=0.1)+
  #geom_line(aes(y=typical),color='blue') +
  labs(y="Reported cases",x="Date") + ylim(c(0,50000)) + 
  scale_x_date(date_breaks = "10 days", date_labels = "%b-%d-%y") +theme_light() +
  scale_color_manual(values = cols) +  theme(axis.text=element_text(size=15),
                                             plot.title = element_text(size=15, face="bold"),
                                             axis.text.x = element_text(angle = 45, hjust = 1),
                                             legend.position = "bottom", legend.title = element_text(size=15),
                                             legend.text = element_text(size=15),
                                             axis.title=element_text(size=15,face="bold")) +
  
  labs(color = " ",title="MB")


gg_MB 

ggsave(file="figs/MB_proj.png", gg_MB, width = 10, height = 8)
saveRDS(gg_MB, file.path("figs/MB-fig.rds"))

#ggplot() + geom_line(data=true_inc,aes(x=date,y=`50%`*2e-05),color='blue',size=1.2,alpha=0.5) +
#  geom_ribbon(data=true_inc,aes(x=date,ymin=`2.5%`*2e-05,ymax=`97.5%`*2e-05),fill='blue',alpha=0.1) +
#  labs(y="% infectious",x="date") + ylim(c(0,6)) + 
#  scale_x_date(date_breaks = "10 days", date_labels = "%d-%m-%Y") 





