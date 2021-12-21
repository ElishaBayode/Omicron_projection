
library("deSolve")
require(ggplot2)
require(reshape2)

sveirs <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    #wf=0.2 # NOTE - this was to test the impact of recovered people being more immune than 
    # vaccinated people. i think it probably makes sense - after all they *just* recovered .
    #c <- 1# effectiveness of NPIs, set as 1, change later to c(t)
    c <- 1 - stngcy/(2+ exp(-0.5*(time-eff_t)))   #intervention 
    N <- S+Er+Em+Ir+Im+R+V+Erv+Emv+Irv+Imv+Rv+W+Erw+Emw+Irw+Imw+Rw #total population 
    lambda_r <- c*beta_r*(Ir + Irv + Irw) #force of infection resident strain
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
    dRv <-  gamma*(Irv + Imv) -wf* (epsilon_r*lambda_r + epsilon_m*lambda_m)*Rv/N - (mu + w2)*Rv
    dW <-   b*ve*V + w2*Rw - (epsilon_r*lambda_r + epsilon_m*lambda_m)*W/N -(mu+ w3)*W
    dErw <- epsilon_r*lambda_r*W/N + wf*epsilon_r*lambda_r*Rw/N - (sigma+mu)*Erw 
    dEmw <- epsilon_m*lambda_m*W/N + wf*epsilon_m*lambda_m*Rw/N - (sigma+mu)*Emw 
    dIrw <- sigma*Erw - (gamma + mu)*Irw
    dImw <- sigma*Emw - (gamma + mu)*Imw
    dRw <-  gamma*(Irw + Imw) - wf*(epsilon_r*lambda_r + epsilon_m*lambda_m)*Rw/N - (mu + w2)*Rw
    return(list(c(dS,dEr,dEm,dIr,dIm,dR,dV,dErv,dEmv,dIrv,dImv,dRv,dW,dErw,dEmw,dIrw,dImw,dRw)))
  })
}

# this just pulls out some incidence values 
get_total_incidence = function(output, parameters) {
  with(as.list( parameters), {
    incid =  output %>% mutate(inc_res = ascFrac*sigma*(Er+Erv+Erw), 
                               inc_mut = ascFrac*sigma*(Em +Emv +Emw), 
                               inc_tot = ascFrac*sigma*(Er+Erv+Erw+Em +Emv +Emw), 
                               inc_vax = sigma*(Erv+Erw + Emv +Emw), 
                               inc_nonvax = sigma*(Er+Em)) %>% 
      select(time, inc_res, inc_mut, inc_tot, inc_vax, inc_nonvax)
    return(incid)})
}


get_vax = function(output) {
  vax = output %>% mutate(vaxtot = V+ Erv + Emv+ 
                            Irv+ Imv+ Rv+ W +Erw+Emw+ Irw+ Imw+ Rw,
                          vaxrecent = V+ Erv + Emv+ 
                            Irv+ Imv+ Rv, waned = W +Erw+Emw+ Irw+ Imw+ Rw) %>% 
    select(time, vaxtot, vaxrecent, waned) 
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
