require(deSolve)
require(ggplot2)
require(reshape2)
library(lubridate)
library(dplyr)
library(data.table)
set.seed(3242)

sveirs <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    #wf=0.2 # NOTE - this was to test the impact of recovered people being more immune than 
    # vaccinated people. i think it probably makes sense - after all they *just* recovered .
    #c <- 1# effectiveness of NPIs, set as 1, change later to c(t)
    c <- 1 - stngcy/(1+ exp(-1.25*(time-eff_t)))   #intervention 
    N <- S+Er+Em+Ir+Im+R+V+Erv+Emv+Irv+Imv+Rv+W+Erw+Emw+Irw+Imw+Rw #total population 
    lambda_r <- c*beta_r*(Ir + Irv + Irw)
    lambda_m <- c*beta_m*(Im + Imv + Imw) #force of infection mutant strain
    
    
    dS <-  mu*N - (lambda_r+lambda_m)*S/N  + w1*R -(mu + nu*ve)*S
    dEr <- lambda_r*S/N + wf* epsilon_r*lambda_r*R/N - (sigma+mu)*Er 
    dEm <- lambda_m*S/N + wf* epsilon_m*lambda_m*R/N - (sigma+mu)*Em
    dIr <- sigma*Er - (gamma + mu)*Ir
    dIm <- sigma*Em - (gamma + mu)*Im
    dR <-  gamma*(Ir + Im) - wf*(epsilon_r*lambda_r + epsilon_m*lambda_m)*R/N - (mu + w1)*R
    dV <-  nu*ve*S + w2*Rv+w3*W - (epsilon_r*lambda_r + epsilon_m*lambda_m)*V/N - (mu + b*ve)*V 
    dErv <- epsilon_r*lambda_r*V/N + wf*epsilon_r*lambda_r*Rv/N - (sigma+mu)*Erv 
    dEmv <- epsilon_m*lambda_m*V/N +wf* epsilon_m*lambda_m*Rv/N - (sigma+mu)*Emv 
    dIrv <- sigma*Erv - (gamma + mu)*Irv
    dImv <- sigma*Emv - (gamma + mu)*Imv
    dRv <-  gamma*(Irv + Imv) -wf* (epsilon_r*lambda_r + epsilon_m*lambda_m)*Rv/N - (mu + w2 + b*ve)*Rv
    dW <-   b*ve*V + w2*Rw - (1-beff)*(epsilon_r*lambda_r + epsilon_m*lambda_m)*W/N -(mu+ w3)*W
    dErw <-(1-beff)*epsilon_r*lambda_r*W/N + (1-beff)*wf*epsilon_r*lambda_r*Rw/N - (sigma+mu)*Erw 
    dEmw <- (1-beff)*epsilon_m*lambda_m*W/N + (1-beff)*wf*epsilon_m*lambda_m*Rw/N - (sigma+mu)*Emw 
    dIrw <- sigma*Erw - (gamma + mu)*Irw
    dImw <- sigma*Emw - (gamma + mu)*Imw
    dRw <-  b*ve*Rv + gamma*(Irw + Imw) - wf*(1-beff)*(epsilon_r*lambda_r + epsilon_m*lambda_m)*Rw/N - (mu + w2)*Rw
    return(list(c(dS,dEr,dEm,dIr,dIm,dR,dV,dErv,dEmv,dIrv,dImv,dRv,dW,dErw,dEmw,dIrw,dImw,dRw)))
  })
}


# this function takes in some intuitive parameters and attempts to create a starting
# point for the ODE that sort of reflects them. 
make_init = function( N=N_pop, vaxlevel = vaxlevel_in,
                      port_wane = port_wane_in, 
                      past_infection = past_infection_in, incres = incres_in, incmut = incmut_in, 
                      pars=as.list(parameters)) {
  ff=6/7 # fudge factor . hard to get incidence right since it depends on other pars too (2/3)
  Vtot = vaxlevel*N*(1-port_wane) # allocate to V, Ev, Iv
  Wtot = vaxlevel*N*port_wane # allocate to W, Ew, Iw 
  # some have had covid. but they might also have been vaccinated. 
  # set up Rs 
  R0 = N*past_infection*(1-vaxlevel)  # past infections, but not vaccinated 
  Rv0 = N*past_infection*(vaxlevel) * (1-port_wane) # recovered, vaxd and not waned 
  Rw0 = N*past_infection*(vaxlevel) * port_wane # rec, vaxd, waned 
  # set up Es : resident 
  Ertot = ff*incres/pars$sigma # total Er 
  Ervw = vaxlevel*pars$ve*Ertot # vaccinated 
  Erw0 = Ervw * port_wane # waned
  Erv0 = Ervw *(1-port_wane) # not waned.  these two add to Ervwboth
  Er0 = (1-vaxlevel*pars$epsilon_r)*Ertot # unvaccinated 
  # set up Es : mutant  
  Emtot = ff*incmut/pars$sigma # total Er 
  Emvw = vaxlevel*pars$ve*Emtot # vaccinated 
  Emw0 = Emvw * port_wane # waned
  Emv0 = Emvw *(1-port_wane) # not waned.  these two add to Ervwboth
  Em0 = (1-vaxlevel*pars$epsilon_m)*Emtot # unvaccinated 
  # set up Is. for now just make them 2x Es since they last twice as long 
  Ir0=ff*2*Er0; Irv0 = ff*2*Erv0; Irw0 = ff*2*Erw0
  Im0 = ff*2*Em0; Imv0=ff*2*Emv0; Imw0= ff*2*Emw0
  
  # set up V0, W0
  # first line plus Rv0 adds to the V total but we have to leave room for the E, I
  V0 = N*(1-past_infection)*vaxlevel*(1-port_wane) - 
    (Erv0+Irv0+Emv0+Imv0)*(1-port_wane)
  W0 = N*(1-past_infection)*vaxlevel*(port_wane) - 
    ( Erw0+Irw0+Emw0+Imw0)*(port_wane) 
  # set up S 
  initmost = c(Er=Er0, Em =Em0, Ir=Ir0, Im=Im0, R=R0, 
               V=V0, Erv=Erv0, Emv=Emv0, Irv = Irv0,  Imv=Imv0, Rv=Rv0, 
               W=W0, Erw=Erw0, Emw=Emw0, Irw = Irw0,  Imw=Imw0, Rw=Rw0)
  state = c(S=N-sum(initmost), initmost)
  return(state)
}
# ----

lag_func <- function(x, k = 1, pad = NA){
  if(k == 0)
    return(x)
  nas <- rep(pad, min(length(x), abs(k)))
  if(k < 0)
    c(tail(x, k), nas) else c(nas, head(x, -k))
}

# this just pulls out some incidence values 
get_total_incidence = function(output, parameters, lag = 0) {
  with(as.list( parameters), {
    incid =  output %>% mutate(inc_res = ascFrac*sigma*lag_func(Er+Erv+Erw, k=lag), 
                               inc_mut = ascFrac*sigma*lag_func(Em +Emv +Emw, k=lag), 
                               inc_tot = ascFrac*sigma*lag_func(Er+Erv+Erw+Em +Emv +Emw, k=lag), 
                               inc_vax = ascFrac*sigma*lag_func(Erv+Erw + Emv +Emw, k=lag), 
                               inc_nonvax = ascFrac*sigma*lag_func(Er+Em), k=lag) %>% 
      select(time, inc_res, inc_mut, inc_tot, inc_vax, inc_nonvax)
    return(incid)})
}


get_vax = function(output) {
  vax = output %>% mutate(vaxtot = V+ Erv + Emv+ 
                            Irv+ Imv+ Rv+ W +Erw+Emw+ Irw+ Imw+ Rw,
                          vaxtwodose = V+ Erv + Emv+ 
                            Irv+ Imv+ Rv, boosted = W +Erw+Emw+ Irw+ Imw+ Rw) %>% 
    select(time, vaxtot, vaxtwodose, boosted) 
  return(vax) 
}

get_growth_rate = function(output, startoffset = 7, duration = 20) {
  tots = output %>% mutate( res = Er +Erv + Erw, 
                            mut = Em + Emv + Emw) %>% select(time, res, mut) %>% 
    filter(time > min(time) + startoffset & time < min(time)+startoffset+duration)
  return(list( resrate =  lm(log(tots$res) ~ tots$time)$coefficients[2], 
               mutrate = lm(log(tots$mut) ~ tots$time)$coefficients[2]))
}

get_population_immunity = function(output,N) {
  vax_induced = output %>% mutate(vax = V + W)
  inf_induced = output %>% mutate(inf = R + Rv + Rw)
  return(list(vax_induced=vax_induced$vax/N,inf_induced=inf_induced$inf/N))
}


get_doubling_time = function(growth_rate){ 
  resdoubling = log(2)/growth_rate$resrate
  mutdoubling = log(2)/growth_rate$mutrate
  return(list(resdoubling,mutdoubling))
}

get_selection_coef = function(growth_rate){ 
  sele_coef = growth_rate$mut - growth_rate$res  
  return(sele_coef)
}