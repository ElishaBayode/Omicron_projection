# original init 
source("analysis/functions.R")
#init <- c(S=N-132,Er=100, Em= 5,Ir=20,
#          Im=0,R=0,V=0,Erv=3,Emv=5,Irv=2,Imv=0, 
#          Rv=0,W=0,Erw=1,Emw=1, Irw=0, Imw=0,Rw=0)


# ----
# this function takes in some intuitive parameters and attempts to create a starting
# point for the ODE that sort of reflects them. 
make_init = function( N=1.174e6, vaxlevel = 0.72,
                      port_wane = 0.1, 
                      past_infection = 0.1, incres = 50, incmut = 4, 
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


N=1.174e6
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
output_SK <- as.data.frame(ode(y = make_init(), times = times, func = sveirs , parms = parameters))
output_SK <- output_SK %>% mutate(total_pop =S+Er+Em+Ir+Im+R+V+Erv+Emv+Irv+Imv+Rv+W+Erw+Emw+Irw+Imw+Rw)


#helper fuctions 

intro_date <- ymd("2021-12-23") 

incid = get_total_incidence(output=output_SK,parameters) #set output to Province output 
incid = incid %>% select(time, inc_res, inc_mut, inc_tot)
incid = incid %>% mutate(date=seq.Date(ymd(intro_date),ymd(intro_date)-1+length(times), 1))
incid = incid %>% mutate(rcases = MASS::rnegbin(length(incid$inc_tot), incid$inc_tot, 
                                                theta=mean(fit$post$phi)))

incid_long <- melt(incid,  id.vars = 'time', variable.name = 'series')
ggplot(incid_long, aes(time,value)) + geom_line(aes(colour = series))


growth_rate_SK <- get_growth_rate(output_SK)
doubling_time_SK <- get_doubling_time( growth_rate_SK)
#mutate(rcases = MASS::rnegbin(length(cases), cases, 
#               theta=mean(fit$post$phi)))
saveRDS(incid, file.path("data/SK-incid.rds"))
incid =readRDS(file.path("data/SK-incid.rds"))

#set path to recent PHAC forecasts data
fit = readRDS("/Users/elishaare/Desktop/PHAC_forecasts_14DEC2021/COVID-PHAC-forecasts/data-generated/SK-fit.rds")
dat = readRDS("/Users/elishaare/Desktop/PHAC_forecasts_14DEC2021/COVID-PHAC-forecasts/data-generated/SK-dat.rds")


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
  geom_line(data=incid, aes(x=date, y=inc_mut, col ="Omicron"), size=1.4) +
  geom_line(data=incid, aes(x=date, y=inc_tot, col ="Total"),  size=1.4) +
  geom_point(data=incid, aes(x=date, y=rcases), color="grey15", size=1.4, alpha=0.3,) +
  
  #geom_ribbon(data = btout, aes(x=date, ymax = y_rep_0.95, ymin=y_rep_0.05), 
  #            alpha=0.3, fill="red")+
  coord_cartesian(ylim = c(0, 2000), expand = FALSE) + 
  scale_x_date(date_breaks = "months", date_labels = "%b") +theme_light() +
  scale_color_manual(values = cols) +  theme(axis.text=element_text(size=15),
                                             plot.title = element_text(size=15, face="bold"),
                                             legend.position = "bottom", legend.title = element_text(size=15),
                                             legend.text = element_text(size=15),
                                             axis.title=element_text(size=15,face="bold")) +
  labs(color = " ",title="SK")

ggsave(file="figs/SK_proj.png", gg, width = 10, height = 8)
