# original init 
source("analysis/functions.R")
#init <- c(S=N-132,Er=100, Em= 5,Ir=20,
#          Im=0,R=0,V=0,Erv=3,Emv=5,Irv=2,Imv=0, 
#          Rv=0,W=0,Erw=1,Emw=1, Irw=0, Imw=0,Rw=0)

# create bc dat object. 
 source("analysis/get-data-bc.R") # NOTE check the last line of data, may need to remove.

# ----- covidseir stuff ---- 


#set path to recent PHAC forecasts data
#fit = readRDS("/Users/elishaare/Desktop/PHAC_forecasts_14DEC2021/COVID-PHAC-forecasts/data-generated/BC-fit.rds")
#dat = readRDS("/Users/elishaare/Desktop/PHAC_forecasts_14DEC2021/COVID-PHAC-forecasts/data-generated/BC-dat.rds")


first_day <- ymd("2020-06-22")# min(dat$date)
lut <- dplyr::tibble(
  day = seq_len(300),
  date = seq(first_day, first_day + length(day) - 1, by = "1 day")
)



proj2 <- project_seir(fit, iter = seq_len(50), stan_model = stan_mod, forecast_days = 0)
proj2$day <- proj2$day + min(dat$day) - 1

modelproj = proj2 %>% tidy_seir()  
modelproj$date = modelproj$day + lut$date[1] - 1 







# ----
# NOTE  I moved make_init to the functions file


N=5.07e6
times <- seq(0,120,1)
ascFrac <- 0.6
intro_date <- ymd("2021-12-02") # i think this is actually the start date for the simulation

eff_date <- ymd("2021-12-22") # date when measures are half effective. I have changed the steepness of the measures function
# so that it is steeper, and this is close to the date when measures take place (less of a range in time) 

# ---- pars ---- 
parameters <-         c(sigma=1/3, # incubation period (3 days) (to fixed)
                        gamma=1/(4), #recovery rate (fixed)
                        nu =0.007, #vax rate: 0.7% per day (fixed)
                        mu=1/(82*365), # 1/life expectancy (fixed)
                        w1= 1/(3*365),# waning rate from R to S (fixed)
                        w2= 1/(3*365), # waning rate from Rv to V (fixed)
                        w3= 1/(3*365),# waning rate Rw to W (fixed)
                        ve=1, # I think this should be 1. it is not really efficacy  ( fixed)
                        #beta_r=0.72, #transmission rate (to estimate) (0.35)
                        beta_m=0.78*2.2, #transmission rate (to estimate)(*1.9)
                        epsilon_r = (1-0.8), # % this should be 1-ve 
                        epsilon_m = (1-0.6), # % escape capacity #(fixed)
                        b= 0.006, # booster rate  (fixed)
                        wf=0.2, # protection for newly recoverd #0.2
                        stngcy= 0,#0.78, #(2*%(reduction)) strength of intervention (reduction in beta's)
                        eff_t = as.numeric(eff_date - intro_date) # time to 50% intervention effectiveness
                         
)

parameters_int <-       c(sigma=1/3, # incubation period (3 days) (to fixed)
                        gamma=1/(4), #recovery rate (fixed)
                        nu =0.007, #vax rate: 0.7% per day (fixed)
                        mu=1/(82*365), # 1/life expectancy (fixed)
                        w1= 1/(3*365),# waning rate from R to S (fixed)
                        w2= 1/(3*365), # waning rate from Rv to V (fixed)
                        w3= 1/(3*365),# waning rate Rw to W (fixed)
                        ve=1, # I think this should be 1. it is not really efficacy  ( fixed)
                        #beta_r=0.72, #transmission rate (to estimate) (0.35)
                        beta_m=0.78*2.2, #transmission rate (to estimate)(*1.9)
                        epsilon_r = (1-0.8), # % this should be 1-ve 
                        epsilon_m = (1-0.6), # % escape capacity #(fixed)
                        b= 0.006, # booster rate  (fixed)
                        wf=0.2, # protection for newly recoverd #0.2
                        stngcy= 0.5,#0.78, #(*%(reduction)) strength of intervention (reduction in beta's)
                        eff_t = as.numeric(eff_date - intro_date) # time to 50% intervention effectiveness
)


# ---- sample betar and get odesolver output ----
m <- 0.78
s <- 0.25
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

# ---- reshape and collect output ---- 
dd <- dplyr::bind_rows(lapply(res,as.data.frame),.id="beta_r")
dd <- select(dd, -beta_r)

dd_int <- dplyr::bind_rows(lapply(res_int,as.data.frame),.id="beta_r")
dd_int <- select(dd_int, -beta_r)
#dd$beta_r <- as.numeric(dd$beta_r)


ag <- aggregate(. ~ time, dd, function(x) c(mean = mean(x), sd = sd(x)))
ag_int <- aggregate(. ~ time, dd_int, function(x) c(mean = mean(x), sd = sd(x)))

output_BC <- data.frame(time=ag$time ,S=ag$S[,1],Er=ag$Er[,1],Em=ag$Em[,1],
                        Ir=ag$Ir[,1],Im=ag$Im[,1],R=ag$R[,1],
                        V=ag$V[,1],Erv=ag$Erv[,1],Emv=ag$Emv[,1],
                        Irv=ag$Irv[,1],Imv=ag$Imv[,1],Rv=ag$Rv[,1],
                        W=ag$W[,1],Erw=ag$Erw[,1],Emw=ag$Emw[,1],
                        Irw=ag$Irw[,1],Imw=ag$Imw[,1],Rw=ag$Rw[,1])

output_BC_sd <- data.frame(time=ag$time ,S=ag$S[,2],Er=ag$Er[,2],Em=ag$Em[,2],
                           Ir=ag$Ir[,2],Im=ag$Im[,2],R=ag$R[,2],
                           V=ag$V[,2],Erv=ag$Erv[,2],Emv=ag$Emv[,2],
                           Irv=ag$Irv[,2],Imv=ag$Imv[,1],Rv=ag$Rv[,2],
                           W=ag$W[,2],Erw=ag$Erw[,2],Emw=ag$Emw[,2],
                           Irw=ag$Irw[,2],Imw=ag$Imw[,2],Rw=ag$Rw[,2])

output_BC_int <- data.frame(time=ag_int$time ,S=ag_int$S[,1],Er=ag_int$Er[,1],Em=ag_int$Em[,1],
                        Ir=ag_int$Ir[,1],Im=ag_int$Im[,1],R=ag_int$R[,1],
                        V=ag_int$V[,1],Erv=ag_int$Erv[,1],Emv=ag_int$Emv[,1],
                        Irv=ag_int$Irv[,1],Imv=ag_int$Imv[,1],Rv=ag_int$Rv[,1],
                        W=ag_int$W[,1],Erw=ag_int$Erw[,1],Emw=ag_int$Emw[,1],
                        Irw=ag_int$Irw[,1],Imw=ag_int$Imw[,1],Rw=ag_int$Rw[,1])

output_BC_sd_int <- data.frame(time=ag_int$time ,S=ag_int$S[,2],Er=ag_int$Er[,2],Em=ag_int$Em[,2],
                           Ir=ag_int$Ir[,2],Im=ag_int$Im[,2],R=ag_int$R[,2],
                           V=ag_int$V[,2],Erv=ag_int$Erv[,2],Emv=ag_int$Emv[,2],
                           Irv=ag_int$Irv[,2],Imv=ag_int$Imv[,1],Rv=ag_int$Rv[,2],
                           W=ag_int$W[,2],Erw=ag_int$Erw[,2],Emw=ag_int$Emw[,2],
                           Irw=ag_int$Irw[,2],Imw=ag_int$Imw[,2],Rw=ag_int$Rw[,2])

incid = get_total_incidence(output=output_BC,parameters=c(mean(beta_r),parameters)) #set output to Province output 
incid = incid %>% select(time, inc_res, inc_mut, inc_tot)
incid = incid %>% mutate(date=seq.Date(ymd(intro_date),ymd(intro_date)-1+length(times), 1))
incid = incid %>% mutate(rcases = MASS::rnegbin(length(incid$inc_tot), incid$inc_tot, 
                                                theta=mean(fit$post$phi)))

incid_sd_BC = get_total_incidence(output=output_BC_sd,parameters=c(mean(beta_r),parameters)) #set output to Province output 
incid_sd_BC = incid %>% select(time, inc_res, inc_mut, inc_tot)
incid_sd_BC = incid_sd_BC %>% mutate(err=qt(0.975,df=sim_size-1)*inc_tot/sqrt(sim_size))
incid_sd_BC = incid_sd_BC %>% mutate(lower = inc_tot-err, upper = inc_tot+err)

incid_int = get_total_incidence(output=output_BC_int,parameters=c(mean(beta_r),parameters_int)) #set output to Province output 
incid_int = incid_int %>% select(time, inc_res, inc_mut, inc_tot)
incid_int = incid_int %>% mutate(date=seq.Date(ymd(intro_date),ymd(intro_date)-1+length(times), 1))
incid_int = incid_int %>% mutate(rcases = MASS::rnegbin(length(incid$inc_tot), incid$inc_tot, 
                                                theta=mean(fit$post$phi)))

incid_sd_BC_int = get_total_incidence(output=output_BC_sd_int,parameters=c(mean(beta_r),parameters_int)) #set output to Province output 
incid_sd_BC_int = incid_int %>% select(time, inc_res, inc_mut, inc_tot)
incid_sd_BC_int = incid_sd_BC_int %>% mutate(err=qt(0.975,df=sim_size-1)*inc_tot/sqrt(sim_size))
incid_sd_BC_int = incid_sd_BC_int %>% mutate(lower = inc_tot-err, upper = inc_tot+err)

ggplot(data = incid, aes(x=date, y=inc_tot))+geom_line(color="red") + 
  geom_line(data = incid_int, aes(x=date, y=inc_tot), color="blue")

incid_long <- melt(incid,  id.vars = 'time', variable.name = 'series')
# ggplot(incid_long, aes(time,value)) + geom_line(aes(colour = series))


# ---- make ggplot with covidseir fit and new model fit ---- 
cols <- c("Current" = "#D55E00", "75% reduction in transmission"="#009E73")
mindate = ymd("2021-10-01")
maxdate = max(incid$date)
gg_BC <- plot_projection(modelproj, dat, date_column = "date") +
  theme(axis.title.x = element_blank()) +
  # geom_line(data=incid, aes(x=date, y=inc_mut, col ="Omicron"), size=1.4) +
  geom_ribbon(data=incid, aes(x=date, ymin = incid_sd_BC$lower, ymax = incid_sd_BC$upper),
             inherit.aes = FALSE, fill = "grey", alpha = 0.5) +
  
  geom_line(data=incid_int, aes(x=date, y=inc_tot, col ="Reduction in transmission"), 
            size=1, alpha = 0.5) +
  geom_ribbon(data=incid_int, aes(x=date, ymin = incid_sd_BC_int$lower, ymax = incid_sd_BC_int$upper),
              inherit.aes = FALSE, fill = "blue", alpha = 0.1) +
  geom_line(data=incid, aes(x=date, y=inc_tot, col ="Current"), 
            size=1, alpha = 0.5) +
  coord_cartesian(ylim = c(0, 15000), xlim=c(mindate, maxdate), expand = FALSE) + 
  scale_x_date(date_breaks = "months", date_labels = "%b") +theme_light() +
  scale_color_manual(values = cols) +  theme(axis.text=element_text(size=15),
                                             plot.title = element_text(size=15, face="bold"),
                                             legend.position = "bottom", legend.title = element_text(size=15),
                                             legend.text = element_text(size=15),
                                             axis.title=element_text(size=15,face="bold")) +
  labs(color = " ",title="BC")

gg_BC



ggsave(file="figs/BC_proj.png", gg_BC, width = 10, height = 8)
saveRDS(gg_BC, file.path("figs/BC-fig.rds"))
#saveRDS(g, file.path("figs/Show_freqBC.rds"))



gg#proj_plot <- fan_plot(fit = fit, pred = modelproj, obs = dat)

#check transm multiplier 
sig = 1:500
f_sig = 1 - 0.8*2/(2+ exp(-0.5*(sig-30)))
plot(sig,f_sig)  

