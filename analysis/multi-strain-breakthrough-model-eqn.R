# original init 
init <- c(S=N-132,Er=100, Em= 5,Ir=20,
          Im=0,R=0,V=0,Erv=3,Emv=5,Irv=2,Imv=0, 
          Rv=0,W=0,Erw=1,Emw=1, Irw=0, Imw=0,Rw=0)

library("deSolve")
require(ggplot2)
require(reshape2)

sveirs <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    wf=0.2 # NOTE - this was to test the impact of recovered people being more immune than 
    # vaccinated people. i think it probably makes sense - after all they *just* recovered .
    c <- 1# effectiveness of NPIs, set as 1, change later to c(t)
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




# ----

# this function pulls out vaccinated and waned so we can easily compare to reported vaccination levels 
get_vax = function(output) {
  vax = output %>% mutate(vaxtot = V+ Erv + Emv+ 
                            Irv+ Imv+ Rv+ W +Erw+Emw+ Irw+ Imw+ Rw,
                          vaxrecent = V+ Erv + Emv+ 
                            Irv+ Imv+ Rv, waned = W +Erw+Emw+ Irw+ Imw+ Rw) %>% 
    select(time, vaxtot, vaxrecent, waned) 
  return(vax) 
}

# pulls out the growth rates in incidence in a time window in a simulation, 
# intended for setting up the model to have observed growth rates 
get_growth_rate = function(output, startoffset = 7, duration = 20) {
  tots = output %>% mutate( res = Er +Erv + Erw, 
                            mut = Em + Emv + Emw) %>% select(time, res, mut) %>% 
    filter(time > min(time) + startoffset & time < min(time)+startoffset+duration)
  return(list( resrate =  lm( log(tots$res) ~ tots$time)$coefficients[2], 
               mutrate = lm(log(tots$mut)~ tots$time)$coefficients[2]))
}

# this function takes in some intuitive parameters and attempts to create a starting
# point for the ODE that sort of reflects them. 
make_init = function( N=5.07e6, vaxlevel = 0.8,
                      port_wane = 0.1, 
                      past_infection = 0.1, incres = 400, incmut = 10, 
                      pars=as.list(parameters)) {
  ff=2/3 # fudge factor . hard to get incidence right since it depends on other pars too (2/3)
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


N=5.07e6
times <- seq(0,60,1)
ascFrac <- 0.6

# ---- pars ---- 
parameters <-         c(sigma=1/3, # incubation period (3 days) (to fixed)
                        gamma=1/(6), #recovery rate (fixed)
                        nu =0.007, #vax rate: 0.7% per day (fixed)
                        mu=1/(82*365), # 1/life expectancy (fixed)
                        w1= 1/(3*365),# waning rate from R to S (fixed)
                        w2= 1/(3*365), # waning rate from Rv to V (fixed)
                        w3= 1/(3*365),# waning rate Rw to W (fixed)
                        ve=1, # I think this should be 1. it is not really efficacy  ( fixed)
                        beta_r=0.5, #transmission rate (to estimate) (0.35)
                        beta_m=0.5*1.7, #transmission rate (to estimate)(*1.9)
                        epsilon_r = (1-0.8), # % this should be 1-ve 
                        epsilon_m = (1-0.6), # % escape capacity #(fixed)
                        b= 0.006 # booster rate  (fixed)
)



#odesolver output
output <- as.data.frame(ode(y = make_init(), times = times, func = sveirs , parms = parameters))
output <- output %>% mutate(total_pop =S+Er+Em+Ir+Im+R+V+Erv+Emv+Irv+Imv+Rv+W+Erw+Emw+Irw+Imw+Rw)


#helper fuctions 
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

growth_rate <- get_growth_rate(output)

get_doubling_time = function(growth_rate){ 
  resdoubling <- log(2)/growth_rate$resrate
  mutdoubling <- log(2)/growth_rate$mutrate
  return(list(resdoubling,mutdoubling))
}



incid = get_total_incidence(output,parameters)
incid = incid %>% select(time, inc_res, inc_mut, inc_tot)
incid = incid %>% mutate(date=seq.Date(ymd("2021-12-20"),ymd("2021-12-20")-1+length(times), 1))
incid = incid %>% mutate(rcases = MASS::rnegbin(length(incid$inc_tot), incid$inc_tot, 
                              theta=mean(fit$post$phi)))

incid_long <- melt(incid,  id.vars = 'time', variable.name = 'series')
ggplot(incid_long, aes(time,value)) + geom_line(aes(colour = series))

vax = get_vax(output) 
ggplot(vax, aes(x=time, y=vaxtot/N)) + geom_line()


#mutate(rcases = MASS::rnegbin(length(cases), cases, 
#               theta=mean(fit$post$phi)))

fit = readRDS("/Users/elishaare/Desktop/PHAC_forecasts_14DEC2021/COVID-PHAC-forecasts/data-generated/BC-fit.rds")
dat = readRDS("/Users/elishaare/Desktop/PHAC_forecasts_14DEC2021/COVID-PHAC-forecasts/data-generated/BC-dat.rds")


first_day <- min(dat$date)
lut <- dplyr::tibble(
  day = seq_len(300),
  date = seq(first_day, first_day + length(day) - 1, by = "1 day")
)



proj2 <- project_seir(fit, iter = seq_len(50), stan_model = stan_mod, forecast_days = 80)
proj2$day <- proj2$day + min(dat$day) - 1

modelproj = proj2 %>% tidy_seir()  
modelproj$date = modelproj$day + lut$date[1] - 1 



cols <- c("Omicron" = "#009E73", "Total"="#D55E00")
gg <- plot_projection(modelproj, dat, date_column = "date") +
  theme(axis.title.x = element_blank()) +
  geom_line(data=incid, aes(x=date, y=inc_tot, col ="Total"),  size=1.2) +
  geom_point(data=incid, aes(x=date, y=rcases), color="grey15", size=1.2, alpha=0.3,) +
  geom_line(data=incid, aes(x=date, y=inc_mut, col ="Omicron"), size=1.2) +
  #geom_ribbon(data = btout, aes(x=date, ymax = y_rep_0.95, ymin=y_rep_0.05), 
  #            alpha=0.3, fill="red")+
  coord_cartesian(ylim = c(0, 5000), expand = FALSE) + 
  scale_x_date(date_breaks = "months", date_labels = "%b") +theme_light() +
  scale_color_manual(values = cols) +  theme(axis.text=element_text(size=15),
                                             plot.title = element_text(size=15, face="bold"),
                                             legend.position = "bottom", legend.title = element_text(size=15),
                                             legend.text = element_text(size=15),
                                             axis.title=element_text(size=15,face="bold")) +
  labs(color = " ")





gg#proj_plot <- fan_plot(fit = fit, pred = modelproj, obs = dat)






#plot
# output_long <- melt(output,  id.vars = 'time', variable.name = 'series')
# ggplot(output_long, aes(time,value)) + geom_line(aes(colour = series))

# needs: 
# 1. accessors to get total cases, total wt and total omicron
# 2. input functions: create the starting state out of the info for how many people
# have had covid and how many we think have been vaccinated 
# 3. sample the solutions and make projections with ribbons, for total daily incidence 
# with each type and total 
