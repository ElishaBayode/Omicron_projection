# original init 
source("analysis/functions.R")
#init <- c(S=N-132,Er=100, Em= 5,Ir=20,
#          Im=0,R=0,V=0,Erv=3,Emv=5,Irv=2,Imv=0, 
#          Rv=0,W=0,Erw=1,Emw=1, Irw=0, Imw=0,Rw=0)


# ----
# this function takes in some intuitive parameters and attempts to create a starting
# point for the ODE that sort of reflects them. 
make_init = function( N=14.57e6, vaxlevel = 0.75,
                      port_wane = 0.1, 
                      past_infection = 0.15, incres = 1700, incmut = 40, 
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


N=14.57e6
times <- seq(0,60,1)
ascFrac <- 0.65
intro_date <- ymd("2021-12-08")
eff_date <- ymd("2022-01-01")
# ---- pars ---- 

# ---- pars ---- 
parameters <-         c(sigma=1/3, # incubation period (3 days) (to fixed)
                        gamma=1/(4), #recovery rate (fixed)
                        nu =0.007, #vax rate: 0.7% per day (fixed)
                        mu=1/(82*365), # 1/life expectancy (fixed)
                        w1= 1/(3*365),# waning rate from R to S (fixed)
                        w2= 1/(3*365), # waning rate from Rv to V (fixed)
                        w3= 1/(3*365),# waning rate Rw to W (fixed)
                        ve=1, # I think this should be 1. it is not really efficacy  ( fixed)
                        #beta_r=0.69, #transmission rate (to estimate) (0.35)
                        beta_m=0.97*2.2, #transmission rate (to estimate)(*1.9)
                        epsilon_r = (1-0.8), # % this should be 1-ve 
                        epsilon_m = (1-0.67), # % escape capacity #(fixed)
                        b= 0.006, # booster rate  (fixed)
                        wf=0.2,
                        stngcy= 2*0,#0.78, #(2*%(reduction)) strength of intervention (reduction in beta's)
                        eff_t = as.numeric(eff_date - intro_date) # time to
)



parameters_int <-         c(sigma=1/3, # incubation period (3 days) (to fixed)
                            gamma=1/(4), #recovery rate (fixed)
                            nu =0.007, #vax rate: 0.7% per day (fixed)
                            mu=1/(82*365), # 1/life expectancy (fixed)
                            w1= 1/(3*365),# waning rate from R to S (fixed)
                            w2= 1/(3*365), # waning rate from Rv to V (fixed)
                            w3= 1/(3*365),# waning rate Rw to W (fixed)
                            ve=1, # I think this should be 1. it is not really efficacy  ( fixed)
                            #beta_r=0.72, #transmission rate (to estimate) (0.35)
                            beta_m=0.97*2.2, #transmission rate (to estimate)(*1.9)
                            epsilon_r = (1-0.8), # % this should be 1-ve 
                            epsilon_m = (1-0.6), # % escape capacity #(fixed)
                            b= 0.006, # booster rate  (fixed)
                            wf=0.2, # protection for newly recoverd #0.2
                            stngcy= 2*0.75,#0.78, #(2*%(reduction)) strength of intervention (reduction in beta's)
                            eff_t = as.numeric(eff_date - intro_date) # time to 50% intervention effectiveness
                            
)




#sig = 1:500
#f_sig = 1 - 0.8*2/(2+ exp(-0.5*(sig-30)))
#plot(sig,f_sig) 

#odesolver output
m <- 0.97
s <- 0.3
location <- log(m^2 / sqrt(s^2 + m^2))
shape <- sqrt(log(1 + (s^2 / m^2)))
print(paste("location:", location))
print(paste("shape:", shape))
sim_size = 15
beta_r <-  rlnorm(n=sim_size, location, shape)

res <- vector(length(beta_r),mode="list")
for (k in seq_along(beta_r)){ #range of values for beta
  res[[k]] <- deSolve::ode(y=make_init(),times=times,func=sveirs,
                           parms=c(beta_r=beta_r[k], parameters))
}


res_int <- vector(length(beta_r),mode="list")
for (k in seq_along(beta_r)){ #range of values for beta
  res_int[[k]] <- deSolve::ode(y=make_init(),times=times,func=sveirs,
                               parms=c(beta_r=beta_r[k], parameters_int))
}


#names(res) <- beta_r
dd <- dplyr::bind_rows(lapply(res,as.data.frame),.id="beta_r")
dd <- select(dd, -beta_r)

dd_int <- dplyr::bind_rows(lapply(res_int,as.data.frame),.id="beta_r")
dd_int <- select(dd_int, -beta_r)
#dd$beta_r <- as.numeric(dd$beta_r)


ag <- aggregate(. ~ time, dd, function(x) c(mean = mean(x), sd = sd(x)))

ag_int <- aggregate(. ~ time, dd_int, function(x) c(mean = mean(x), sd = sd(x)))

output_ON <- data.frame(time=ag$time ,S=ag$S[,1],Er=ag$Er[,1],Em=ag$Em[,1],
                        Ir=ag$Ir[,1],Im=ag$Im[,1],R=ag$R[,1],
                        V=ag$V[,1],Erv=ag$Erv[,1],Emv=ag$Emv[,1],
                        Irv=ag$Irv[,1],Imv=ag$Imv[,1],Rv=ag$Rv[,1],
                        W=ag$W[,1],Erw=ag$Erw[,1],Emw=ag$Emw[,1],
                        Irw=ag$Irw[,1],Imw=ag$Imw[,1],Rw=ag$Rw[,1])


output_ON_sd <- data.frame(time=ag$time ,S=ag$S[,2],Er=ag$Er[,2],Em=ag$Em[,2],
                           Ir=ag$Ir[,2],Im=ag$Im[,2],R=ag$R[,2],
                           V=ag$V[,2],Erv=ag$Erv[,2],Emv=ag$Emv[,2],
                           Irv=ag$Irv[,2],Imv=ag$Imv[,1],Rv=ag$Rv[,2],
                           W=ag$W[,2],Erw=ag$Erw[,2],Emw=ag$Emw[,2],
                           Irw=ag$Irw[,2],Imw=ag$Imw[,2],Rw=ag$Rw[,2])





output_ON_int <- data.frame(time=ag_int$time ,S=ag_int$S[,1],Er=ag_int$Er[,1],Em=ag_int$Em[,1],
                            Ir=ag_int$Ir[,1],Im=ag_int$Im[,1],R=ag_int$R[,1],
                            V=ag_int$V[,1],Erv=ag_int$Erv[,1],Emv=ag_int$Emv[,1],
                            Irv=ag_int$Irv[,1],Imv=ag_int$Imv[,1],Rv=ag_int$Rv[,1],
                            W=ag_int$W[,1],Erw=ag_int$Erw[,1],Emw=ag_int$Emw[,1],
                            Irw=ag_int$Irw[,1],Imw=ag_int$Imw[,1],Rw=ag_int$Rw[,1])


output_ON_sd_int <- data.frame(time=ag_int$time ,S=ag_int$S[,2],Er=ag_int$Er[,2],Em=ag_int$Em[,2],
                               Ir=ag_int$Ir[,2],Im=ag_int$Im[,2],R=ag_int$R[,2],
                               V=ag_int$V[,2],Erv=ag_int$Erv[,2],Emv=ag_int$Emv[,2],
                               Irv=ag_int$Irv[,2],Imv=ag_int$Imv[,1],Rv=ag_int$Rv[,2],
                               W=ag_int$W[,2],Erw=ag_int$Erw[,2],Emw=ag_int$Emw[,2],
                               Irw=ag_int$Irw[,2],Imw=ag_int$Imw[,2],Rw=ag_int$Rw[,2])






#odesolver output
#output_ON <- as.data.frame(ode(y = make_init(), times = times, func = sveirs , parms = parameters))
output_ON <- output_ON %>% mutate(total_pop =S+Er+Em+Ir+Im+R+V+Erv+Emv+Irv+Imv+Rv+W+Erw+Emw+Irw+Imw+Rw)


#helper fuctions 

fit = readRDS("/Users/elishaare/Desktop/PHAC_forecasts_14DEC2021/COVID-PHAC-forecasts/data-generated/ON-fit.rds")
dat = readRDS("/Users/elishaare/Desktop/PHAC_forecasts_28DEC2021/COVID-PHAC-forecasts/data-generated/ON-dat.rds")


incid = get_total_incidence(output=output_ON,parameters=c(mean(beta_r),parameters)) #set output to Province output 
incid = incid %>% select(time, inc_res, inc_mut, inc_tot)
incid = incid %>% mutate(date=seq.Date(ymd(intro_date),ymd(intro_date)-1+length(times), 1))
incid = incid %>% mutate(rcases = MASS::rnegbin(length(incid$inc_tot), incid$inc_tot, 
                                                theta=mean(fit$post$phi)))


incid_sd_ON = get_total_incidence(output=output_ON_sd,parameters=c(mean(beta_r),parameters)) #set output to Province output 
incid_sd_ON = incid %>% select(time, inc_res, inc_mut, inc_tot)
incid_sd_ON = incid_sd_ON %>% mutate(err=qt(0.975,df=sim_size-1)*inc_tot/sqrt(sim_size))
incid_sd_ON = incid_sd_ON %>% mutate(lower = inc_tot-err, upper = inc_tot+err)

incid_int = get_total_incidence(output=output_ON_int,parameters=c(mean(beta_r),parameters_int)) #set output to Province output 
incid_int = incid_int %>% select(time, inc_res, inc_mut, inc_tot)
incid_int = incid_int %>% mutate(date=seq.Date(ymd(intro_date),ymd(intro_date)-1+length(times), 1))
incid_int = incid_int %>% mutate(rcases = MASS::rnegbin(length(incid$inc_tot), incid$inc_tot, 
                                                        theta=mean(fit$post$phi)))


incid_sd_ON_int = get_total_incidence(output=output_ON_sd_int,parameters=c(mean(beta_r),parameters_int)) #set output to Province output 
incid_sd_ON_int = incid_int %>% select(time, inc_res, inc_mut, inc_tot)
incid_sd_ON_int = incid_sd_ON_int %>% mutate(err=qt(0.975,df=sim_size-1)*inc_tot/sqrt(sim_size))
incid_sd_ON_int = incid_sd_ON_int %>% mutate(lower = inc_tot-err, upper = inc_tot+err)



#ON_Omicron_proj <- select(incid, c("date","inc_mut" , "rcases","inc_tot"))
#names(ON_Omicron_proj) <- c("date", "Omicron", "rcases","total")
#ON_Omicron_proj$Omicron <- round(ON_Omicron_proj$Omicron,2)
#ON_Omicron_proj$total <- round(ON_Omicron_proj$total,2)
#ON_Omicron_proj$rcases <- round(ON_Omicron_proj$rcases,2)
#ON_Omicron_proj <- filter(ON_Omicron_proj, date >= ymd("2021-12-10") & date <= ymd("2022-01-10"))
#write.csv(ON_Omicron_proj,'ON_Omicron_proj.csv')


incid_long <- melt(incid,  id.vars = 'time', variONle.name = 'series')
ggplot(incid_long, aes(time,value)) + geom_line(aes(colour = series))


growth_rate_ON <- get_growth_rate(output_ON)
doubling_time_ON <- get_doubling_time( growth_rate_ON)
get_selection_coef(growth_rate_ON)
check <- get_population_immunity(output_ON,N)
plot(NA,NA, xlim=c(0,60), ylim=c(0,1))
lines(check$vax_induced)
lines(check$inf_induced)
lines(check$inf_induced+check$vax_induced)

#mutate(rcases = MASS::rnegbin(length(cases), cases, 
#               theta=mean(fit$post$phi)))
saveRDS(incid, file.path("data/ON-incid.rds"))
incid =readRDS(file.path("data/ON-incid.rds"))

saveRDS(incid_sd_ON, file.path("data/ON-incid_er.rds"))
incid_sd_ON =readRDS(file.path("data/ON-incid_er.rds"))


saveRDS(incid_int, file.path("data/ON-incid_int.rds"))
incid_int =readRDS(file.path("data/ON-incid_int.rds"))

saveRDS(incid_sd_ON_int, file.path("data/ON-incid_er_int.rds"))
incid_sd_ON_int =readRDS(file.path("data/ON-incid_er_int.rds"))


#set path to recent PHAC forecasts data
#fit = readRDS("/Users/elishaare/Desktop/PHAC_forecasts_14DEC2021/COVID-PHAC-forecasts/data-generated/ON-fit.rds")
#dat = readRDS("/Users/elishaare/Desktop/PHAC_forecasts_14DEC2021/COVID-PHAC-forecasts/data-generated/ON-dat.rds")


first_day <- min(dat$date)
lut <- dplyr::tibble(
  day = seq_len(300),
  date = seq(first_day, first_day + length(day) - 1, by = "1 day")
)



proj2 <- project_seir(fit, iter = seq_len(50), stan_model = stan_mod, forecast_days = 0)
proj2$day <- proj2$day + min(dat$day) - 1

modelproj = proj2 %>% tidy_seir()  
modelproj$date = modelproj$day + lut$date[1] - 1 



cols <- c("Current" = "#D55E00", "75% reduction in transmission"="#009E73")

#cols <- c( "Total"="#D55E00") #"Omicron" = "#009E73",
gg_ON <- plot_projection(modelproj, dat, date_column = "date") +
  theme(axis.title.x = element_blank()) +
  # geom_line(data=incid, aes(x=date, y=inc_mut, col ="Omicron"), size=1.4) +
  geom_ribbon(data=incid, aes(x=date, ymin = incid_sd_ON$lower, ymax = incid_sd_ON$upper),
              inherit.aes = FALSE, fill = "grey", alpha = 0.5) +
  
  geom_line(data=incid_int, aes(x=date, y=inc_tot, col ="75% reduction in transmission"), 
            size=1, alpha = 0.5) +
  geom_ribbon(data=incid_int, aes(x=date, ymin = incid_sd_ON_int$lower, ymax = incid_sd_ON_int$upper),
              inherit.aes = FALSE, fill = "blue", alpha = 0.1) +
  geom_line(data=incid, aes(x=date, y=inc_tot, col ="Current"), 
            size=1, alpha = 0.5) +
  
  
  # geom_point(data=incid, aes(x=date, y=rcases), color="grey15", size=1.4, alpha=0.5) +
  
  #geom_ribbon(data = btout, aes(x=date, ymax = y_rep_0.95, ymin=y_rep_0.05), 
  #            alpha=0.3, fill="red")+
  coord_cartesian(ylim = c(0, 40000), expand = FALSE) + 
  scale_x_date(date_breaks = "months", date_labels = "%b") +theme_light() +
  scale_color_manual(values = cols) +  theme(axis.text=element_text(size=15),
                                             plot.title = element_text(size=15, face="bold"),
                                             legend.position = "bottom", legend.title = element_text(size=15),
                                             legend.text = element_text(size=15),
                                             axis.title=element_text(size=15,face="bold")) +
  
  labs(color = " ",title="ON")



saveRDS(gg_ON, file.path("figs/ON-fig.rds"))
ggsave(file="figs/ON_proj.png", gg_ON, width = 10, height = 8)
