library(ggplot2)
library(lubridate)
#install.packages("pomp")
library(pomp)
library(dplyr)
library(plyr)

source("analysis/functions.R")

#importing data 
dat = readRDS("data/BC-dat.rds")

#including Omicron wave only 
intro_date <-  ymd("2021-11-28")
dat_omic <- filter(dat, date >= intro_date) %>% select(c("day", "value"))
dat_omic$day <- 1:nrow(dat_omic)

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
  DEmv  = epsilon_m*c*beta_m*(Im + Imv + Imw)*V/N +wf* epsilon_m*c*beta_m*(Im + Imv + Imw)*Rv/N - (sigma+mu)*Emv; 
  DIrv  =  sigma*Erv - (gamma + mu)*Irv;
  DImv  = sigma*Emv - (gamma + mu)*Imv;
  DRv  =  gamma*(Irv + Imv) -wf* (epsilon_r*c*beta_r*(Ir + Irv + Irw) + epsilon_m*c*beta_m*(Im + Imv + Imw))*Rv/N - (mu + w2)*Rv;
  DW  =   b*ve*V + w2*Rw - (epsilon_r*c*beta_r*(Ir + Irv + Irw) + epsilon_m*c*beta_m*(Im + Imv + Imw))*W/N -(mu+ w3)*W;
  DErw = epsilon_r*c*beta_r*(Ir + Irv + Irw)*W/N + wf*epsilon_r*c*beta_r*(Ir + Irv + Irw)*Rw/N - (sigma+mu)*Erw; 
  DEmw = epsilon_m*c*beta_m*(Im + Imv + Imw)*W/N + wf*epsilon_m*c*beta_m*(Im + Imv + Imw)*Rw/N - (sigma+mu)*Emw; 
  DIrw  = sigma*Erw - (gamma + mu)*Irw;
  DImw = sigma*Emw - (gamma + mu)*Imw;
  DRw  = gamma*(Irw + Imw) - wf*(epsilon_r*c*beta_r*(Ir + Irv + Irw) + epsilon_m*c*beta_m*(Im + Imv + Imw))*Rw/N - (mu + w2)*Rw;
 ")),
  rinit=Csnippet("
      S = S_0; V=V_0;       W= W_0;
      Er = Er_0;  Erv=Erv_0;  Erw = Erv_0;
      Em = Em_0;  Emv = Emv_0; Emw=Emw_0;
      Ir = Ir_0;  Irv = Irv_0; Irw=Irw_0;
      Im = Im_0;  Imv = Imv_0; Imw = Imw_0;
      R =  R_0;   Rv = Rv_0 ;  Rw = Rw_0;"),
  paramnames=c("mu","N","c","beta_r", "beta_m", "epsilon_r", "epsilon_m", "sigma", "wf", "gamma", "w1", "w2", "w3",
               "ve", "nu","b","S_0","Er_0","Em_0", "Ir_0", "Im_0", "R_0", "V_0", "Erv_0",
               "Emv_0", "Irv_0", "Imv_0", "Rv_0", "W_0", "Erw_0", "Emw_0", "Irw_0","Imw_0", "Rw_0"),
  statenames=c("S","Er","Em", "Ir", "Im", "R", "V", "Erv",
               "Emv", "Irv", "Imv", "Rv", "W", "Erw", "Emw", "Irw","Imw", "Rw")) -> BC_obj


#declearing fixed parameters 

parameters <-         c(sigma=1/3, # incubation period (3 days) (to fixed)
                        gamma=1/(5), #recovery rate (fixed)
                        nu =0.007, #vax rate: 0.7% per day (fixed)
                        mu=1/(82*365), # 1/life expectancy (fixed)
                        w1= 1/(3*365),# waning rate from R to S (fixed)
                        w2= 1/(3*365), # waning rate from Rv to V (fixed)
                        w3= 1/(3*365),# waning rate Rw to W (fixed)
                        ve=1, # I think this should be 1. it is not really efficacy  ( fixed)
                        #beta_r=0.72, #transmission rate (to estimate) (0.35)
                        #    beta_m=0.8*2.2, #transmission rate (to estimate)(*1.9)
                        epsilon_r = (1-0.8), # % this should be 1-ve 
                        epsilon_m = (1-0.6), # % escape capacity #(fixed)
                        b= 0.006, # booster rate  (fixed)
                        wf=0.2, # protection for newly recoverd #0.2
                        N=5e6,
                        c=1
)


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


#define a likelihood function 
#theta is the over dispersion parameter 
#parameter p is sampling efficiency and ascertainment prob 

negbin.loglik <- function (params) {
  x <- trajectory(BC_obj,params=params)
  prediction <-  (x["Er",,]+ x["Erv",,] + x["Erw",,]
                   + x["Em",,]+ x["Emv",,] + x["Emw",,])
 sum(dnbinom(x=obs(BC_obj),
              mu=params["p"]*prediction,size=1/params["theta"],
              log=TRUE))
}

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




guess <- c(log(1.03), logit(0.65), log(1.34), log(0.01)) #c(log(10),log(15),log(1))

#the parameters are constrained  accordingly (lower and upper)

fit_BC <- optim(fn=f_loglik,par=guess, lower=c(log(0.6), 0.5, log(1.2), log(0.1)), 
              upper = c(log(2), 0.7, log(2.3), log(0)), method = "L-BFGS-B")


#this catches estimated parameter values from MLE 
mle_est_BC <- c(beta_r=exp(fit_BC$par[1]),p=expit(fit_BC$par[2]), beta_m=exp(fit_BC$par[3]),theta=exp(fit_BC$par[4]))

signif(mle_est_BC,3)


#estimated parameter values will now be used for for prediction 
coef(BC_obj) <- c(c(S_0=init[[1]],Er_0=init[[2]],Em_0=init[[3]],Ir_0=init[[4]],
                    Im_0=init[[5]],R_0=init[[6]],V_0=init[[7]],Erv_0=init[[8]], 
                    Emv_0=init[[9]],Irv_0=init[[10]],Imv_0=init[[11]],Rv_0=init[[12]],
                    W_0=init[[13]],Erw_0=init[[14]],Emw_0=init[[15]],Irw_0=init[[16]],
                    Imw_0=init[[17]],Rw_0=init[[18]] ,c(parameters,mle_est_BC)))


#an attempt to model changes in ascertainment probability 
test_prop <- 1 - 0.5/(1+ exp(-1.25*(1:nrow(dat_omic)-32)))

# model predictions (incidence, sigma*(Er +Em + Erv + Emv + Erw + Emw)) 

model.pred <- test_prop*parameters[[1]]*(trajectory(BC_obj)["Er",,]+ 
              trajectory(BC_obj)["Erv",,] +trajectory(BC_obj)["Erw",,]+
              trajectory(BC_obj)["Em",,]+trajectory(BC_obj)["Emv",,] + 
              trajectory(BC_obj)["Emw",,])


#this generates multiple  model realizations 

raply(10000,rnbinom(n=length(model.pred),
                     mu=coef(BC_obj,"p")*model.pred,
                     size=1/coef(BC_obj,"theta"))) -> simdat

#this generates appriopriate quantiles

aaply(simdat,2,quantile,probs=c(0.025,0.5,0.975)) -> quantiles


#sample a typical model realization (sanity check)
typ <- sample(nrow(simdat),1) 


dat_sim <- cbind(as.data.frame(BC_obj),
                 quantiles,
                 typical=simdat[typ,]) 

#create dates 
dat_sim = dat_sim %>% mutate(date=seq.Date(ymd(intro_date),
                      ymd(intro_date)-1+nrow(dat_omic), 1))

ggplot(data=dat_sim,
       mapping=aes(x=date))+
  geom_line(aes(y=`50%`),color='blue',size=1.2,alpha=0.5)+
  geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`),fill='blue',alpha=0.1)+
  geom_point(aes(y=value),color='black')+
 geom_line(aes(y=typical),color='blue') +
  labs(y="cases",x="date") + ylim(c(0,7000))
 



