library(ggplot2)
library(lubridate)
#install.packages("pomp")
library(pomp)
library(dplyr)
library(plyr)

dat_omic <- filter(dat, date >= ymd("2021-12-05")) %>% select(c("day", "value"))
dat_omic$day <- 1:nrow(dat_omic)

ascFrac <- 0.6 

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

negbin.loglik <- function (params) {
  x <- trajectory(BC_obj,params=params)
  prediction <-  (x["Er",,]+ x["Erv",,] + x["Erw",,]
                   + x["Em",,]+ x["Emv",,] + x["Emw",,])
 sum(dnbinom(x=obs(BC_obj),
              mu=ascFrac*prediction,size=1/params["theta"],
              log=TRUE))
}


parameters <-         c(sigma=1/3, # incubation period (3 days) (to fixed)
                        gamma=1/(5), #recovery rate (fixed)
                        nu =0.007, #vax rate: 0.7% per day (fixed)
                        mu=1/(82*365), # 1/life expectancy (fixed)
                        w1= 1/(3*365),# waning rate from R to S (fixed)
                        w2= 1/(3*365), # waning rate from Rv to V (fixed)
                        w3= 1/(3*365),# waning rate Rw to W (fixed)
                        ve=1, # I think this should be 1. it is not really efficacy  ( fixed)
                        #beta_r=0.72, #transmission rate (to estimate) (0.35)
                        beta_m=0.8*2, #transmission rate (to estimate)(*1.9)
                        epsilon_r = (1-0.8), # % this should be 1-ve 
                        epsilon_m = (1-0.6), # % escape capacity #(fixed)
                        b= 0.006, # booster rate  (fixed)
                        wf=0.2, # protection for newly recoverd #0.2
                        N=5e6,
                        c=1
)

f7 <- function (par) {
  params <- c(S_0=make_init()[[1]],Er_0=make_init()[[2]],Em_0=make_init()[[3]],Ir_0=make_init()[[4]],
              Im_0=make_init()[[5]],R_0=make_init()[[6]],V_0=make_init()[[7]],Erv_0=make_init()[[8]], 
              Emv_0=make_init()[[9]],Irv_0=make_init()[[10]],Imv_0=make_init()[[11]],Rv_0=make_init()[[12]],
              W_0=make_init()[[13]],Erw_0=make_init()[[14]],Emw_0=make_init()[[15]],Irw_0=make_init()[[16]],
              Imw_0=make_init()[[17]],Rw_0=make_init()[[18]] ,parameters,
              beta_r=exp(par[1]),theta=exp(par[2]))
  -negbin.loglik(params)
}

guess <- c(log(10),log(1))

fit7 <- optim(fn=f7,par=guess); fit7
mle3 <- c(beta_r=exp(fit7$par[1]),theta=exp(fit7$par[2]))
signif(mle3,3)


coef(BC_obj) <- c(c(S_0=make_init()[[1]],Er_0=make_init()[[2]],Em_0=make_init()[[3]],Ir_0=make_init()[[4]],
                    Im_0=make_init()[[5]],R_0=make_init()[[6]],V_0=make_init()[[7]],Erv_0=make_init()[[8]], 
                    Emv_0=make_init()[[9]],Irv_0=make_init()[[10]],Imv_0=make_init()[[11]],Rv_0=make_init()[[12]],
                    W_0=make_init()[[13]],Erw_0=make_init()[[14]],Emw_0=make_init()[[15]],Irw_0=make_init()[[16]],
                    Imw_0=make_init()[[17]],Rw_0=make_init()[[18]] ,parameters,mle3))

model.pred <- parameters[[1]]*(trajectory(BC_obj)["Er",,]+ trajectory(BC_obj)["Erv",,] +trajectory(BC_obj)["Erw",,]+
                                 trajectory(BC_obj)["Em",,]+trajectory(BC_obj)["Emv",,] + 
                             trajectory(BC_obj)["Emw",,])


raply(100000,rnbinom(n=length(model.pred),
                     mu=ascFrac*model.pred,
                     size=1/coef(BC_obj,"theta"))) -> simdat

aaply(simdat,2,quantile,probs=c(0.025,0.5,0.975)) -> quantiles


typ <- sample(nrow(simdat),1)

ggplot(data=cbind(as.data.frame(BC_obj),
                  quantiles,
                  typical=simdat[typ,]),
       mapping=aes(x=day))+
  geom_line(aes(y=`50%`),color='red')+
  geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`),fill='red',alpha=0.2)+
  geom_point(aes(y=value),color='black')+
 # geom_line(aes(y=typical),color='blue')
  labs(y="cases",x="days") #+ ylim(c(0,4000))
 

