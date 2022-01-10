source("analysis/mod_fitting_setup.R")
#setting values to generate initial conditions with make_init()
N=5.07e6
N_pop=N
#ascFrac <- 0.5
vaxlevel_in = 0.8
port_wane_in = 0.04 
past_infection_in = 0.1
incres_in = 300
incmut_in = 10
init <- make_init() #to generatee initial state 
#importing data 
dat = readRDS("data/BC-dat.rds")

#including Omicron wave only 
intro_date <-  ymd("2021-11-28")
dat_omic <- filter(dat, date >= intro_date) %>% select(c("day", "value"))
dat_omic$day <- 1:nrow(dat_omic)


#declaring fixed parameters 

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
                        wf=0.2, # protection for newly recoverd #0.2
                        N=5e6,
                        c=1
)




guess <- c(log(1.03), logit(0.65), log(1.34), log(0.01)) #c(log(10),log(15),log(1))

#the parameters are constrained  accordingly (lower and upper)

fit_BC <- optim(fn=f_loglik,par=guess, lower=c(log(0.6), 0.5, log(1.2), log(0.1)), 
                upper = c(log(2), 0.7, log(2.3), log(0)), method = "L-BFGS-B")


#this catches estimated parameter values from MLE 
mle_est_BC <- c(beta_r=exp(fit_BC$par[1]),p=expit(fit_BC$par[2]), beta_m=exp(fit_BC$par[3]),theta=exp(fit_BC$par[4]))

signif(mle_est_BC,3)


#estimated parameter values will now be used for for prediction 
coef(pomp_obj ) <- c(c(S_0=init[[1]],Er_0=init[[2]],Em_0=init[[3]],Ir_0=init[[4]],
                       Im_0=init[[5]],R_0=init[[6]],V_0=init[[7]],Erv_0=init[[8]], 
                       Emv_0=init[[9]],Irv_0=init[[10]],Imv_0=init[[11]],Rv_0=init[[12]],
                       W_0=init[[13]],Erw_0=init[[14]],Emw_0=init[[15]],Irw_0=init[[16]],
                       Imw_0=init[[17]],Rw_0=init[[18]] ,c(parameters,mle_est_BC)))

estim_parameters <- c(parameters,mle_est_BC,stngcy=0, eff_t=ymd("2022-02-01"),beff=0.7)

#an attempt to model changes in ascertainment probability 
test_prop <- 1 - 0.65/(1+ exp(-1.25*(1:nrow(dat_omic)-35)))

# model predictions (incidence, sigma*(Er +Em + Erv + Emv + Erw + Emw)) 

out_state <-  trajectory(pomp_obj )

model.pred <- test_prop*parameters[[1]]*(out_state["Er",,]+ 
                                           out_state["Erv",,] + out_state["Erw",,]+
                                           out_state["Em",,]+ out_state["Emv",,] + 
                                           out_state["Emw",,])


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

#create dates 
dat_sim = dat_sim %>% mutate(date=seq.Date(ymd(intro_date),
                                           ymd(intro_date)-1+nrow(dat_omic), 1))



######## forecast cases ############### 
check <- case_projection(out_state=out_state,
                         forecast_days=30, 
                         lag = 0,
                         parameters_estim = estim_parameters, 
                         simu_size = 100000,
                         pomp_obj = pomp_obj , 
                         dat_sim = dat_sim
)


ggplot() + geom_line(data=dat_sim,aes(x=date,y=`50%`),color='blue',size=1.2,alpha=0.5) +
  geom_ribbon(data=dat_sim,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='blue',alpha=0.1)+
  geom_point(data=dat_sim,aes(x=date, y=value),color='black') +
  geom_line(data=check,aes(x=date, y=`50%`),color='blue',size=1.2,alpha=0.5) +
  geom_ribbon(data=check,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='blue',alpha=0.1)+
  #geom_line(aes(y=typical),color='blue') +
  labs(y="cases",x="date") + ylim(c(0,50000))
