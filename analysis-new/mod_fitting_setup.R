library(ggplot2)
library(lubridate)
#install.packages("pomp")
library(pomp)
library(dplyr)
library(plyr)

source("analysis/functions.R")



#create a pomp object for fitting

pomp(
  data=dat_omic,
  times="day",t0=0,
  skeleton=vectorfield(
    Csnippet("
  DS  =  mu*N - (c*beta_r*(Ir + Irv + Irw)+c*beta_m*(Im + Imv + Imw))*S/N  + w1*R -(mu + nu*ve)*S;
  DEr  = c*beta_r*(Ir + Irv + Irw)*S/N + wf* epsilon_r*c*beta_r*(Ir + Irv + Irw)*R/N - (sigma+mu)*Er; 
  DEm  = c*beta_m*(Im + Imv + Imw)*S/N + wf* epsilon_m*c*beta_m*(Im + Imv + Imw)*R/N - (sigma+mu)*Em;
  DIr  = sigma*Er - (gamma + mu)*Ir;
  DIm  = sigma*Em - (gamma + mu)*Im;
  DR  = gamma*(Ir + Im) - wf*(epsilon_r*c*beta_r*(Ir + Irv + Irw) + epsilon_m*c*beta_m*(Im + Imv + Imw))*R/N - (mu + w1)*R;
  DV  = nu*ve*S + w2*Rv+w3*W - (epsilon_r*c*beta_r*(Ir + Irv + Irw) + epsilon_m*c*beta_m*(Im + Imv + Imw))*V/N - (mu + b*ve)*V; 
  DErv  = epsilon_r*c*beta_r*(Ir + Irv + Irw)*V/N + wf*epsilon_r*c*beta_r*(Ir + Irv + Irw)*Rv/N - (sigma+mu)*Erv; 
  DEmv  = epsilon_m*c*beta_m*(Im + Imv + Imw)*V/N + wf* epsilon_m*c*beta_m*(Im + Imv + Imw)*Rv/N - (sigma+mu)*Emv; 
  DIrv  =  sigma*Erv - (gamma + mu)*Irv;
  DImv  = sigma*Emv - (gamma + mu)*Imv;
  DRv  =  gamma*(Irv + Imv) -wf* (epsilon_r*c*beta_r*(Ir + Irv + Irw) + epsilon_m*c*beta_m*(Im + Imv + Imw))*Rv/N - (mu + w2)*Rv - b*ve*Rv;
  DW  =   b*ve*V + w2*Rw - (1-beff)*(epsilon_r*c*beta_r*(Ir + Irv + Irw) + epsilon_m*c*beta_m*(Im + Imv + Imw))*W/N -(mu+ w3)*W;
  DErw = (1-beff)*epsilon_r*c*beta_r*(Ir + Irv + Irw)*W/N + (1-beff)*wf*epsilon_r*c*beta_r*(Ir + Irv + Irw)*Rw/N - (sigma+mu)*Erw; 
  DEmw = (1-beff)*epsilon_m*c*beta_m*(Im + Imv + Imw)*W/N + (1-beff)*wf*epsilon_m*c*beta_m*(Im + Imv + Imw)*Rw/N - (sigma+mu)*Emw; 
  DIrw  = sigma*Erw - (gamma + mu)*Irw;
  DImw = sigma*Emw - (gamma + mu)*Imw;
  DRw  = gamma*(Irw + Imw) + b*ve*Rv- wf*(1-beff)*(epsilon_r*c*beta_r*(Ir + Irv + Irw) + epsilon_m*c*beta_m*(Im + Imv + Imw))*Rw/N - (mu + w2)*Rw;
 ")),
  rinit=Csnippet("
      S = S_0; V=V_0;       W= W_0;
      Er = Er_0;  Erv=Erv_0;  Erw = Erv_0;
      Em = Em_0;  Emv = Emv_0; Emw=Emw_0;
      Ir = Ir_0;  Irv = Irv_0; Irw=Irw_0;
      Im = Im_0;  Imv = Imv_0; Imw = Imw_0;
      R =  R_0;   Rv = Rv_0 ;  Rw = Rw_0;"),
  paramnames=c("mu","N","c","beta_r", "beta_m", "epsilon_r", "epsilon_m", "sigma", "wf", "gamma", "w1", "w2", "w3",
               "ve", "nu","b","beff","S_0","Er_0","Em_0", "Ir_0", "Im_0", "R_0", "V_0", "Erv_0",
               "Emv_0", "Irv_0", "Imv_0", "Rv_0", "W_0", "Erw_0", "Emw_0", "Irw_0","Imw_0", "Rw_0"),
  statenames=c("S","Er","Em", "Ir", "Im", "R", "V", "Erv",
               "Emv", "Irv", "Imv", "Rv", "W", "Erw", "Emw", "Irw","Imw", "Rw")) -> pomp_obj 



#define a likelihood function 
#theta is the over dispersion parameter 
#parameter p is sampling efficiency and ascertainment prob 

negbin.loglik <- function (params) {
  x <- trajectory(pomp_obj ,params=params)
  prediction <-  (x["Er",,]+ x["Erv",,] + x["Erw",,]
                   + x["Em",,]+ x["Emv",,] + x["Emw",,])*test_prop
 sum(dnbinom(x=obs(pomp_obj ),
              mu=params["p"]*prediction,size=1/params["theta"],
              log=TRUE))
}


# negbin.loglik_1 <- function (params) {
#   x <- trajectory(pomp_obj ,params=params)
#   prediction <-  (x["Er",,]+ x["Erv",,] + x["Erw",,]
#                   + x["Em",,]+ x["Emv",,] + x["Emw",,])
#   sum(dnbinom(x=obs(pomp_obj ),
#               mu=params["p"]*prediction,size=1/params["theta"],
#               log=TRUE))
# }

#use to generate neg log likelihood estimates  
f_loglik <- function (par) {
  params <- c(S_0=init[[1]],Er_0=init[[2]],Em_0=init[[3]],Ir_0=init[[4]],
              Im_0=init[[5]],R_0=init[[6]],V_0=init[[7]],Erv_0=init[[8]], 
              Emv_0=init[[9]],Irv_0=init[[10]],Imv_0=init[[11]],Rv_0=init[[12]],
              W_0=init[[13]],Erw_0=init[[14]],Emw_0=init[[15]],Irw_0=init[[16]],
              Imw_0=init[[17]],Rw_0=init[[18]] ,parameters,
              beta_r=exp(par[1]),p=expit(par[2]),beta_m=exp(par[3]),theta=exp(par[4]))
  -negbin.loglik(params)
}

# f_loglik_2 <- function (par) {
#   params <- c(S_0=last(out_state["S",,]),Er_0=last(out_state["Er",,]),Em_0=last(out_state["Em",,]),
#               Ir_0=last(out_state["Ir",,]),Im_0=last(out_state["Im",,]), R_0=last(out_state["R",,]),
#               V_0=last(out_state["V",,]),Erv_0=last(out_state["Erv",,]), Emv_0=last(out_state["Emv",,]), 
#               Irv_0=last(out_state["Irv",,]),Imv_0=last(out_state["Imv",,]), Rv_0=last(out_state["Rv",,]),
#               W_0=last(out_state["W",,]),Erw_0=last(out_state["Erw",,]), Emw_0=last(out_state["Emw",,]),
#               Irw_0=last(out_state["Irw",,]),Imw_0=last(out_state["Imw",,]), Rw_0=last(out_state["Rw",,])
#               ,parameters, beta_r=exp(par[1]),p=expit(par[2]),beta_m=exp(par[3]),theta=exp(par[4]))
#   -negbin.loglik_1(params)
# }

case_projection <- function(out_state=out_state,forecast_days=forecast_days, lag = lag ,
                            parameters_estim=estim_parameters, simu_size=simu_size,
                            pomp_obj=pomp_obj , dat_sim=dat_sim,test_prop=test_prop){
  init_current = c(S=last(out_state_2["S",,]),Er=last(out_state_2["Er",,]),Em=last(out_state_2["Em",,]),
                   Ir=last(out_state_2["Ir",,]),Im=last(out_state_2["Im",,]), R=last(out_state_2["R",,]),
                   V=last(out_state_2["V",,]),Erv=last(out_state_2["Erv",,]), Emv=last(out_state_2["Emv",,]), 
                   Irv=last(out_state_2["Irv",,]),Imv=last(out_state_2["Imv",,]), Rv=last(out_state_2["Rv",,]),
                   W=last(out_state_2["W",,]),Erw=last(out_state_2["Erw",,]), Emw=last(out_state_2["Emw",,]),
                   Irw=last(out_state_2["Irw",,]),Imw=last(out_state_2["Imw",,]), Rw=last(out_state_2["Rw",,])) 
  forecast_days = seq(1, forecast_days, 1)
  output = as.data.frame(deSolve::ode(y=init_current,times=forecast_days,func=sveirs,
                                      parms=c(estim_parameters)))       
  
  incidence_forecast =  last(test_prop)*estim_parameters[[1]]*lag_func(output$Er +
                                                      output$Erv+ output$Erw +
                                                      output$Em + output$Emv +
                                                      output$Emw, k=lag)
  
  uncert_bound = raply(simu_size,rnbinom(n=length(incidence_forecast),
                                         mu=coef(pomp_obj ,"p")*incidence_forecast,
                                         size=1/coef(pomp_obj ,"theta")))
  quantiles_proj =  as.data.frame(aaply(uncert_bound ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975)))
  
  project_dat = quantiles_proj %>% mutate(date=seq.Date(ymd(last(dat_sim$date)),
                                                        ymd(last(dat_sim$date))-1+last(forecast_days), 1))
  return(project_dat)
  
}

 
# case_projection_rel <- function(out_state=out_state_2,forecast_days=forecast_days, lag = lag ,
#                             parameters_estim=estim_parameters, simu_size=simu_size,
#                             pomp_obj=pomp_obj , dat_sim=dat_sim, test_prop=test_prop){
#   init_current = c(S=last(out_state_2["S",,]),Er=last(out_state_2["Er",,]),Em=last(out_state_2["Em",,]),
#                    Ir=last(out_state_2["Ir",,]),Im=last(out_state_2["Im",,]), R=last(out_state_2["R",,]),
#                    V=last(out_state_2["V",,]),Erv=last(out_state_2["Erv",,]), Emv=last(out_state_2["Emv",,]), 
#                    Irv=last(out_state_2["Irv",,]),Imv=last(out_state_2["Imv",,]), Rv=last(out_state_2["Rv",,]),
#                    W=last(out_state_2["W",,]),Erw=last(out_state_2["Erw",,]), Emw=last(out_state_2["Emw",,]),
#                    Irw=last(out_state_2["Irw",,]),Imw=last(out_state_2["Imw",,]), Rw=last(out_state_2["Rw",,])) 
#   forecast_days = seq(1, forecast_days, 1)
#   output = as.data.frame(deSolve::ode(y=init_current,times=forecast_days,func=sveirs,
#                                       parms=c(estim_parameters)))       
#   
#   incidence_forecast =  estim_parameters[[1]]*lag_func(output$Er + output$Erv+ output$Erw +
#                                                                          output$Em + output$Emv +
#                                                                          output$Emw, k=lag)
#   
#   uncert_bound = raply(simu_size,rnbinom(n=length(incidence_forecast),
#                                          mu=coef(pomp_obj ,"p")*incidence_forecast,
#                                          size=1/coef(pomp_obj ,"theta")))
#   quantiles_proj =  as.data.frame(aaply(uncert_bound ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975)))
#   
#   project_dat = quantiles_proj %>% mutate(date=seq.Date(ymd(last(dat_sim$date)),
#                                                         ymd(last(dat_sim$date))-1+last(forecast_days), 1))
#   return(project_dat)
#   
# }
# 
# 
# 
# case_projection_true <- function(out_state=out_state_3,forecast_days=forecast_days, lag = lag ,
#                                 parameters_estim=estim_parameters, simu_size=simu_size,
#                                 pomp_obj=pomp_obj , dat_sim=dat_sim, test_prop=test_prop){
#   init_current = c(S=last(out_state_3["S",,]),Er=last(out_state_3["Er",,]),Em=last(out_state_3["Em",,]),
#                    Ir=last(out_state_3["Ir",,]),Im=last(out_state_3["Im",,]), R=last(out_state_3["R",,]),
#                    V=last(out_state_3["V",,]),Erv=last(out_state_3["Erv",,]), Emv=last(out_state_3["Emv",,]), 
#                    Irv=last(out_state_3["Irv",,]),Imv=last(out_state_3["Imv",,]), Rv=last(out_state_3["Rv",,]),
#                    W=last(out_state_3["W",,]),Erw=last(out_state_3["Erw",,]), Emw=last(out_state_3["Emw",,]),
#                    Irw=last(out_state_3["Irw",,]),Imw=last(out_state_3["Imw",,]), Rw=last(out_state_3["Rw",,])) 
#   forecast_days = seq(1, forecast_days, 1)
#   output = as.data.frame(deSolve::ode(y=init_current,times=forecast_days,func=sveirs,
#                                       parms=c(estim_parameters)))       
#   
#   incidence_forecast =  estim_parameters[[1]]*lag_func(output$Er + output$Erv+ output$Erw +
#                                                          output$Em + output$Emv +
#                                                          output$Emw, k=lag)
#   
#   uncert_bound = raply(simu_size,rnbinom(n=length(incidence_forecast),
#                                          mu=incidence_forecast,
#                                          size=1/coef(pomp_obj ,"theta")))
#   quantiles_proj =  as.data.frame(aaply(uncert_bound ,2,quantile,na.rm=TRUE,probs=c(0.025,0.5,0.975)))
#   
#   project_dat = quantiles_proj %>% mutate(date=seq.Date(ymd(last(dat_sim$date)),
#                                                         ymd(last(dat_sim$date))-1+last(forecast_days), 1))
#   return(project_dat)
#   
# }
# 
