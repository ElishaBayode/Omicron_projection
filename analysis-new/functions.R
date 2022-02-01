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
    c <- (1 - stngcy/(1+ exp(-1.25*(time-eff_t))))   #intervention 
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
    dW <-   b*ve*V + w2*Rw - (1-beff)*(lambda_r + lambda_m)*W/N -(mu+ w3)*W
    dErw <-(1-beff)*lambda_r*W/N + (1-beff)*wf*lambda_r*Rw/N - (sigma+mu)*Erw 
    dEmw <- (1-beff)*lambda_m*W/N + (1-beff)*wf*lambda_m*Rw/N - (sigma+mu)*Emw 
    dIrw <- sigma*Erw - (gamma + mu)*Irw
    dImw <- sigma*Emw - (gamma + mu)*Imw
    dRw <-  b*ve*Rv + gamma*(Irw + Imw) - wf*(1-beff)*(lambda_r + lambda_m)*Rw/N - (mu + w2)*Rw
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
# ---- show the simulation in the simplest plot, with the data 







# ----

lag_func <- function(x, k = 1, pad = NA){
  if(k == 0)
    return(x)
  nas <- rep(pad, min(length(x), abs(k)))
  if(k < 0)
    c(tail(x, k), nas) else c(nas, head(x, -k))
}

# ------


# this just pulls out some incidence values 
get_total_incidence = function(output, parameters, lag = 0 ) {
  ascFrac = parameters["p"]
  with(as.list( parameters), {
    incid =  output %>% mutate(inc_res = ascFrac*sigma*lag_func(Er+Erv+Erw, k=lag), 
                               inc_mut = ascFrac*sigma*lag_func(Em +Emv +Emw, k=lag), 
                               inc_tot = ascFrac*sigma*lag_func(Er+Erv+Erw+Em +Emv +Emw, k=lag), 
                               inc_vax = ascFrac*sigma*lag_func(Erv+Erw + Emv +Emw, k=lag), 
                               inc_nonvax = ascFrac*sigma*lag_func(Er+Em), k=lag) %>% 
      select(time, inc_res, inc_mut, 
             inc_tot, inc_vax, inc_nonvax)
    return(incid)})
}



get_true_incidence_plot = function(times, start_date, parameters_base,init) {
  with(as.list( parameters), {
    ascFrac=1 # for TRUE incidence, don't reduce by ascertainment, and don't lag
    lag = 0
    output   = as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                               parms=parameters_base))
    
    incid =  output %>% mutate(inc_res = ascFrac*sigma*lag_func(Er+Erv+Erw, k=lag), 
                               inc_mut = ascFrac*sigma*lag_func(Em +Emv +Emw, k=lag), 
                               inc_tot = ascFrac*sigma*lag_func(Er+Erv+Erw+Em +Emv +Emw, k=lag), 
                               inc_vax = ascFrac*sigma*lag_func(Erv+Erw + Emv +Emw, k=lag), 
                               inc_nonvax = ascFrac*sigma*lag_func(Er+Em), k=lag) %>% 
      select(time, inc_res, inc_mut, inc_tot, inc_vax, inc_nonvax) %>% 
      mutate(date=seq.Date(ymd(start_date),
      ymd(start_date)-1+length(time), 1))
  
      
     cols = c("total" = "orange", "incidence_vax" = "blue","incidence_nonvax" = "darkgreen" )
     plot = ggplot() + geom_line(data=incid,aes(x=date,y=inc_tot, colour="total"),size=1.2,alpha=0.4) +
       geom_line(data=incid,aes(x=date,y=inc_vax, colour = "incidence_vax" ),size=1.2,alpha=0.4) +
            geom_line(data=incid,aes(x=date,y=inc_nonvax, colour="incidence_nonvax"),size=1.2,alpha=0.4) +
           
            labs(y="True incidence",x="Date") + 
             scale_x_date(date_breaks = "10 days", date_labels = "%b-%d-%y") +theme_light() +
             theme(legend.position = "bottom") +
             scale_color_manual(values = cols) + labs(color = " ") 
     return(plot)})

 }


get_true_incidence_prop_plot = function(times, start_date, parameters_base,init) {
  with(as.list( parameters), {
    ascFrac=1 # for TRUE incidence, don't reduce by ascertainment, and don't lag
    lag = 0
    output   = as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                                          parms=parameters_base))
    
    incid_prop =  output %>% mutate(inc_vax_prop = ascFrac*sigma*lag_func((Erv+ Emv)/(V+Erv + Emv+Irv+Imv +Rv), k=lag), 
                  inc_boost_prop = ascFrac*sigma*lag_func((Erw+ Emw )/(W+Erw+Emw+Irw+Imw+Rw), k=lag),
                  inc_totvax_prop = ascFrac*sigma*lag_func((Erw+ Emw +Erv+ Emv)/(W+Erw+Emw+Irw+Imw+Rw+V+Erv + 
                                                                                Emv+Irv+Imv +Rv), k=lag),
                               inc_nonvax_prop = ascFrac*sigma*lag_func((Er+Em)/(S+R+Er+Em), k=lag)) %>% 
      select(time, inc_vax_prop, inc_boost_prop, inc_totvax_prop, inc_nonvax_prop) %>% 
      mutate(date=seq.Date(ymd(start_date),
                           ymd(start_date)-1+length(time), 1))
    
    
    cols = c("incidence_totvax_prop" = "orange", "incidence_boost_prop" = "blue","incidence_vax_prop" = "darkgreen", 
             "incidence_nonvax_prop"="purple" )
    plot = ggplot() + geom_line(data=incid_prop,aes(x=date,y=inc_totvax_prop, colour="incidence_totvax_prop"),size=1.2,alpha=0.4) +
      geom_line(data=incid_prop,aes(x=date,y=inc_boost_prop, colour =  "incidence_boost_prop"),size=1.2,alpha=0.4) +
      geom_line(data=incid_prop,aes(x=date,y=inc_vax_prop, colour="incidence_vax_prop"),size=1.2,alpha=0.4) +
      geom_line(data=incid_prop,aes(x=date,y=inc_nonvax_prop, colour="incidence_nonvax_prop"),size=1.2,alpha=0.4) +
      labs(y="True incidence proportion",x="Date") + 
      scale_x_date(date_breaks = "10 days", date_labels = "%b-%d-%y") +theme_light() +
      theme(legend.position = "bottom") +
      scale_color_manual(values = cols) + labs(color = " ") 
    return(plot)})
  
}






get_vax = function(output) {
  vax = output %>% mutate(unvax = S+Er + Em + Ir + Im + R, 
                          vaxtot = V+ Erv + Emv+ 
                            Irv+ Imv+ Rv+ W +Erw+Emw+ Irw+ Imw+ Rw,
                          vaxtwodose = V+ Erv + Emv+ 
                            Irv+ Imv+ Rv, boosted = W +Erw+Emw+ Irw+ Imw+ Rw) %>% 
    select(time,unvax, vaxtot, vaxtwodose, boosted) 
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

# ---- functions related to test_prop 

# mydat is a data frame with Reported_Date, under70 ("Yes" or "No") and the cases in that group on that day (totcases)
make_case_splines = function(mydat) {
  over70 = filter(mydat, under70=="No")
  logover70 = over70 %>% mutate(lcases = log(totcases)) %>% select(Reported_Date, lcases)
  ospline = smooth.spline(logover70$Reported_Date, logover70$lcases, df=15)
  pred = data.frame( predlcases = ospline$y,
                     cases = over70$totcases, lcases = log(over70$totcases),
                     Reported_Date = logover70$Reported_Date)
  # (2) under 70 
  under70 = filter(mydat, under70=="Yes")
  logunder70 = under70 %>% mutate(lcases = log(totcases)) %>% select(Reported_Date, lcases)
  uspline = smooth.spline(logunder70$Reported_Date, logunder70$lcases, df=15)
  upred = data.frame( predlcases = uspline$y, 
                      cases = under70$totcases, lcases = log(under70$totcases),
                      Reported_Date = logunder70$Reported_Date)
  return(list(pred=pred,upred=upred))
}
# splinetest = make_case_splines(mydat)
# the function returns two data fraems, pred and upred. 
# both have the date, the predicted log cases( from the spline), the cases from the data, and teh log cases from the data 

# do a sanity check for the spline to make sure that it looks ok. The 'df' parameter might need to change for example
# ggplot(data = splinetest$pred, aes(x=Reported_Date, y=lcases))+geom_point(color="blue", alpha=0.5) +
#  geom_line(data =splinetest$pred, aes(x=Reported_Date, y=predlcases))


# next we need a function to use this, a date, and an 'offset' to make test_prop. 


get_testprop = function(changedate, mysplines, halftime, steepness) { 
  # need to get start and end value for the new offset. 
  # the offset starts at the value of the estimated offset on the change date: 
  offstart = filter(mysplines$upred, Reported_Date == changedate)$predlcases -
    filter(mysplines$pred, Reported_Date == changedate)$predlcases 
  # and it will move, over a time about equal to 2*halftime, to its historical value: 
  offend = mean(filter(mysplines$upred, Reported_Date <= changedate)$predlcases -
                  filter(mysplines$pred, Reported_Date <= changedate)$predlcases)
  
  
  offconst = mysplines$upred$predlcases-mysplines$pred$predlcases
  L1 = changedate - min(mysplines$pred$Reported_Date)
  L2 = nrow(mysplines$pred)-L1 
  offconst[(L1+1):nrow(mysplines$pred)] = 
    getoffset(startvalue=offstart,   endvalue = offend, halftime= halftime,
              steepness = steepness, ndays = L2) 
  testpropdf = data.frame(date = mysplines$pred$Reported_Date, 
                          totalmodel = exp(mysplines$pred$predlcases)+ exp(mysplines$upred$predlcases), 
                          totaladjusted = exp(mysplines$pred$predlcases+offconst), 
                          test_prop = pmin(1,(exp(mysplines$pred$predlcases)+
                                                exp(mysplines$upred$predlcases))/ exp(mysplines$pred$predlcases+offconst)))
  return(testpropdf) 
}




#function to get incidence by vaccination status 


# mytest = get_testprop(changedate = ymd("2021-12-21"), 
# mysplines = splinetest, 
# halftime = 20, steepness = 0.1)
# glimpse(mytest)
# ggplot(mytest, aes(x=date, y=test_prop))+geom_line()

getoffset = function(startvalue = 3.73, endvalue=2.68, halftime=15, steepness=0.25, ndays=60) {
  return( startvalue - (startvalue-endvalue)/(1+exp(-steepness*(1:ndays-halftime))))  
}


# ---- two comparison functions -----
compare_two_incid = function(pars1, pars2,name1 = "first", name2="second", 
                             returnplot = T,numsamples = 1e3, numdays = 300, dispar=0.1,mode = "static") { 
  if (mode == "static") { 
    myfunction = sveirs 
  } else {
    myfunction = sveirs_evol
  }
  # run the first 
  out1 <- as.data.frame(deSolve::ode(y=init_BC,time=1:numdays,func= myfunction,
                                     parms=pars1)) 
  inc1 =  pars1[["p"]]*pars1["sigma"]*(out1$Er + out1$Erv + out1$Erw +
                            out1$Em + out1$Emv +
                            out1$Emw)
  unc1 = raply(numsamples,rnbinom(n=length(inc1),
                                  mu=inc1,
                                  size=1/dispar))
  
  proj1 =  as.data.frame(aaply(unc1,2,quantile,
                               na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
    mutate(date=seq.Date(ymd(intro_date),
                         ymd(intro_date)-1+numdays, 1))
  proj1$name = name1
  # run the second 
  out2 <- as.data.frame(deSolve::ode(y=init_BC,time=1:numdays,func= myfunction,
                                     parms=pars2)) 
  inc2 =  pars2[["p"]]*pars2["sigma"]*(out2$Er + out2$Erv + out2$Erw +
                            out2$Em + out2$Emv +
                            out2$Emw)
  unc2 = raply(numsamples,rnbinom(n=length(inc2),
                                  mu=inc2,
                                  size=1/dispar))
  
  proj2 =  as.data.frame(aaply(unc2,2,quantile,
                               na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
    mutate(date=seq.Date(ymd(intro_date),
                         ymd(intro_date)-1+numdays, 1))
  proj2$name = name2
  proj = rbind(proj1, proj2)
  if (returnplot) { 
    gg =   ggplot(data = proj) + geom_line(aes(x=date,y=`50%`, color=name),
                                           size=1.5,alpha=0.6) +
      geom_ribbon(aes(x=date,ymin=`2.5%`,ymax=`97.5%`, fill=name),
                  alpha=0.2, size = 1.5)+
      ylab("Incident infections") + theme_light() + theme(legend.position = "bottom",
                                                          axis.title.x = element_blank(),
                                                          legend.title = element_blank())
    
  } else {
    return(proj) 
  }
}

simple_prev_plot = function(pars1, numdays = 90){
  out1 <- as.data.frame(deSolve::ode(y=init,time=1:numdays,func= sveirs,
                                     parms=pars1)) 
  prev1 =  (out1$Ir + out1$Irv + out1$Irw + out1$Im + out1$Imv +out1$Imw)
  prev= data.frame(date = seq.Date(ymd(intro_date),ymd(intro_date)-1+numdays, 1),
                   preval = prev1)
return(  ggplot(prev, aes(x=date, y=preval))+geom_line() + 
    theme( axis.title.x = element_blank()))
}

compare_two_preval = function(pars1, pars2,name1 = "first", name2="second", 
                              returnplot = T,numsamples = 1e3, numdays = 300, dispar=0.1, 
                              mode = "static") { 
  if (mode == "static") { 
    myfunction = sveirs 
  } else {
    myfunction = sveirs_evol
  }
  # run the first 
  out1 <- as.data.frame(deSolve::ode(y=init_BC,time=1:numdays,func= myfunction,
                                     parms=pars1)) 
  prev1 =  (out1$Ir + out1$Irv + out1$Irw +
              out1$Im + out1$Imv +
              out1$Imw)
  unc1 = raply(numsamples,rnbinom(n=length(prev1),
                                  mu=prev1,
                                  size=1/dispar))
  
  proj1 =  as.data.frame(aaply(unc1,2,quantile,
                               na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
    mutate(date=seq.Date(ymd(intro_date),
                         ymd(intro_date)-1+numdays, 1))
  proj1$name = name1
  # run the second 
  out2 <- as.data.frame(deSolve::ode(y=init_BC,time=1:numdays,func= myfunction,
                                     parms=pars2)) 
  prev2 =  (out2$Ir + out2$Irv + out2$Irw +
              out2$Im + out2$Imv +
              out2$Imw)
  unc2 = raply(numsamples,rnbinom(n=length(prev2),
                                  mu=prev2,
                                  size=1/dispar))
  
  proj2 =  as.data.frame(aaply(unc2,2,quantile,
                               na.rm=TRUE,probs=c(0.025,0.5,0.975))) %>% 
    mutate(date=seq.Date(ymd(intro_date),
                         ymd(intro_date)-1+numdays, 1))
  proj2$name = name2
  proj = rbind(proj1, proj2)
  if (returnplot) { 
    gg =   ggplot(data = proj) + geom_line(aes(x=date,y=`50%`, color=name),
                                           size=1.5,alpha=0.6) +
      geom_ribbon(aes(x=date,ymin=`2.5%`,ymax=`97.5%`, fill=name),
                  alpha=0.2, size = 1.5)+
      ylab("Prevalence") + theme_light() + theme(legend.position = "bottom",
                                                 axis.title.x = element_blank(),
                                                 legend.title = element_blank())
    
  } else {
    return(proj) 
  }
}


sveirs_evol <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    #wf=0.2 # NOTE - this was to test the impact of recovered people being more immune than 
    # vaccinated people. i think it probably makes sense - after all they *just* recovered .
    #c <- 1# effectiveness of NPIs, set as 1, change later to c(t)
    c <- (1 - stngcy/(1+ exp(-1.25*(time-eff_t))))   #intervention 
    eff_evolve = (1 + deltaeff/(1+ exp(-0.05*(time-150))))
    epsilon_m = epsilon_m*eff_evolve # note: higher e_m is higher force of infection 
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
    dW <-   b*ve*V + w2*Rw - (1-beff)*(lambda_r + lambda_m)*W/N -(mu+ w3)*W
    dErw <-(1-beff)*lambda_r*W/N + (1-beff)*wf*lambda_r*Rw/N - (sigma+mu)*Erw 
    dEmw <- (1-beff)*lambda_m*W/N + (1-beff)*wf*lambda_m*Rw/N - (sigma+mu)*Emw 
    dIrw <- sigma*Erw - (gamma + mu)*Irw
    dImw <- sigma*Emw - (gamma + mu)*Imw
    dRw <-  b*ve*Rv + gamma*(Irw + Imw) - wf*(1-beff)*(lambda_r + lambda_m)*Rw/N - (mu + w2)*Rw
    return(list(c(dS,dEr,dEm,dIr,dIm,dR,dV,dErv,dEmv,dIrv,dImv,dRv,dW,dErw,dEmw,dIrw,dImw,dRw)))
  })
}
