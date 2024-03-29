# got BC hosp from Jens' code
source("~/BCCovidSnippets/R/helpers.R")

d<-get_can_covid_tracker_data("bc") %>%
  mutate(Date=as.Date(date))

dd<-get_british_columbia_case_data() %>%
  count(Date=`Reported Date`,name="Cases")

jens = d %>%
  pivot_longer(c("total_hospitalizations","total_criticals","change_cases")) %>%
  select(Date,name,value) %>%
  mutate(value=as.integer(value)) %>%
  bind_rows(dd %>% rename(value=Cases) %>% mutate(name="Cases")) %>%
  filter(Date>=as.Date("2020-08-01")) %>%
  mutate(value=pmax(1,value)) %>%
  group_by(name) %>%
  arrange(desc(Date)) %>%
  filter(cumsum(value)>0) %>%
  mutate(stl=add_stl_trend_m(value)) %>%
  mutate(trend=stl$trend) %>%
  filter(name!="change_cases") %>%
  filter(Date>=as.Date("2020-09-01")) %>%
  mutate(name=recode(name,"total_criticals"="ICU census","total_hospitalizations"="Hospital census"))
glimpse(jens)

tothosp <- jens %>% filter( name == "Hospital census") %>% ungroup() %>% 
  rename( date = Date, hosp= value) %>% 
  select(date, hosp)
ggplot(tothosp, aes(x=date, y = hosp))+geom_point()


# case data 
bcpub <- get_british_columbia_case_data() #readr::read_csv("data-raw/BCCDC_COVID19_Dashboard_Case_Details.csv")
#bcpub <- readr::read_rds("bcpub-2020-09-29.rds")
stop_date <- ymd("2022-03-10")   # I don't know why we would stop on this date but I will 

bcpub$Reported_Date <- lubridate::ymd(bcpub$`Reported Date`)
bcpub$thiscase <- 1
dat <- group_by(bcpub, Reported_Date) %>%
  dplyr::summarise(cases = sum(thiscase)) %>%
  filter(Reported_Date >= ymd("2020-02-27")) %>% 
  dplyr::rename(date = Reported_Date) 
dat <- dat[order(dat$date), ]
dat$day <- seq(1, nrow(dat))
dat$value <- dat$cases

alldat <- group_by(bcpub, Reported_Date) %>%
  dplyr::summarise(cases = sum(thiscase)) %>%
  filter(Reported_Date >= ymd("2020-02-27")) %>% dplyr::rename(date = Reported_Date) 
alldat <- alldat[order(alldat$date), ] %>% select(date, cases)



cashos = merge( dat, tothosp, by = "date") %>% select(date, cases, hosp) %>% 
  mutate(scaledcases = lag(cases, n=18)*0.55) %>%
     pivot_longer(c(3,4), names_to = "type", values_to = "number")
ggplot(cashos, aes(x=date, y = number, color=type))+geom_line() +
  scale_x_date(date_breaks = "months", date_labels = "%b-%d") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "bottom")


# fit the model to the cases and adjust to get a projection for hosp 

# continuing as in the BC projection file 
# setup stuff 
forecasts_days <- 30 
intro_date <-  ymd("2021-11-30")
# come back to this later 
#import data 
# dat = readRDS("data/ON-dat.rds")
#include Omicron wave only
dat <- dat %>% filter(date >= intro_date &  date <= stop_date)
dat_omic <- dat
dat_omic <- filter(dat_omic, date >= intro_date) %>% dplyr::select(c("day", "value", "date"))
dat_omic$day <- 1:nrow(dat_omic)



########## Test_prop set up


# am going to do something very simple. test_prop from 1 to about 0.3 over dec 15-about 25
x=intro_date+1:nrow(dat_omic)
tptest = 1 -0.7/(1+exp(-0.2*as.numeric((x - ymd("2021-12-18")))))
plot(x,tptest)
test_prop= tptest


## parameter and model setup
sum(filter(cashos, type=="hosp", date< ymd("2021-11-30"))$cases) # 425334 or only 8.4% of BC 

N=5.07e6
N_pop=N
vaxlevel_in =  0.88 # 
simu_size = 10000 # number of times to resample the negative binom (for ribbons)
forecasts_days =30 # how long to forecast for 
times = 1:nrow(dat_omic)

#declaring  parameters 
eff_date <-   ymd("2021-12-30")  # intervention date # restrictions were updated on 24 December 
intv_date <-  ymd("2022-02-10") 
fur_intv_date <- ymd("2022-04-03")
parameters <-         c(sigma=1/3, # incubation period (days) 
                        gamma=1/(5), #recovery rate 
                        nu =0.007, #vax rate: 0.7% per day 
                        mu=1/(82*365), # 1/life expectancy 
                        w1= 1/(0.5*365),# waning rate from R to S 
                        w2= 1/(0.5*365), # waning rate from Rv to V 
                        w3= 1/(0.5*365),# waning rate Rw to W 
                        ve=1, # I think this should be 1. it is not really efficacy  
                        beta_r=0.6, #transmission rate 
                        beta_m=1, #transmission rate 
                        epsilon_r = (1-0.8), # % this should be 1-ve 
                        epsilon_m = 1-0.15, # % escape capacity 
                        b= 0.012, # booster rate
                        beff = 0.7, # booster efficacy
                        wf=0.1, # protection for newly recovered
                        N=5.07e6,
                        stngcy= 0.4, #(*%(reduction)) strength of intervention (reduction in beta's)
                        eff_t = as.numeric(eff_date - intro_date),
                        relx_level = 0,
                        fur_relx_level = 0,
                        rlx_t = as.numeric(intv_date - intro_date),
                        fur_rlx_t = as.numeric(fur_intv_date - intro_date),
                        p = 0.5, # asc fraction 
                        theta = 0.1 #negative binomial dispersion
)

# initial conditions
port_wane_in = 0.05 # portion boosted at start time .... on Nov 30
past_infection_in = 0.16  # should be related to p, ascertainment in 2021
incres_in = 335 #866#(made up factor) # resident strain (delta) incidence at start (938, reported cases) 
incmut_in = 20 # new (omicron) inc at stat (presumably very low)

init <- make_init()   #generate initial states



########## Run the model and compare the model to the data at least for IC checking
#          (i.e. without fitting anything)
times <- 1:(nrow(dat_omic) + forecasts_days)
outtest <- as.data.frame(deSolve::ode(y=init,time=times,func= sveirs,
                                      parms=parameters)) 
inctest = get_total_incidence(outtest, parameters = parameters) # has p but not test_prop! 
inctest$date = intro_date-1+1:nrow(inctest)
thisvec=extendtp(nrow(inctest), test_prop)
inctest$inc_reported = inctest$inc_tot*thisvec

# sanity check - should have delta in a decline of about 2% /day (-0.02) and around 
# a 0.2 to 0.25 difference between the two, so omicron at about 0.2 
get_growth_rate(outtest, startoffset = 2, duration = 10)

ggplot(data =inctest, aes(x=date, y = inc_reported))+geom_line() +
  geom_line(aes(x=date, y= inc_res), color = "blue") +
  geom_line(aes(x=date, y= inc_mut), color = "red") +
  geom_point(data = dat, aes(x=date, y=cases), alpha=0.5) +
  ylim(c(0,15000)) + xlim(c(ymd("2021-11-30"), ymd("2022-02-28")))


############## do the fit 

known_prop <- 0.5
date_known_prop <- "2021-12-12" # confirmed with plot p4 of https://www.publichealthontario.ca/-/media/Documents/nCoV/epi/covid-19-sars-cov2-whole-genome-sequencing-epi-summary.pdf?sc_lang=en

# 2. growth advantage of mutant strain
known_growth <- 0.2
period_known_growth <- c("2021-12-05", "2021-12-15")
penalties <- list(known_prop = known_prop, date_known_prop = date_known_prop, 
                  known_growth = known_growth, period_known_growth = period_known_growth)

# Determine weight of penalty. 
# Qu: how strong should penalty be on scale of 0-1? 0 = no penalty. 1 = relatively as impactful as the likelihood
pen.size <- 0.3
guess <- c( beta_m=1, stngcy=0.4,beta_r=0.6) 
pen.fit_BC <- optim(fn=func_penloglik,  par=guess, lower=c(0,0,0,0.001), upper = c(Inf,1,Inf,Inf), 
                    method = "L-BFGS-B", 
                    parameters = parameters, test_prop=test_prop, 
                    pen.size=pen.size, penalties = penalties, dat_omic=dat_omic, hessian=T)
pen.fit_BC

# check the fit 
func_loglik(pen.fit_BC$par, test_prop, dat_omic,parameters) 
parameters[names(guess)] <- pen.fit_BC$par


gg = simple_prev_plot(parameters, numdays = 190, mode = "both"); gg  # cc's simple prevalence plot 


# i want to plot the incidence against reported cases. 
# for this i need to get the incidence out, and adjust by test_prop, and plot 
# that against reported cases 

out_samp <- as.data.frame(deSolve::ode(y=init,time=times,func=sveirs,
                                       parms=parameters)) 
reportable = parameters[["p"]]*parameters[["sigma"]]*(out_samp$Er + out_samp$Erv + out_samp$Erw +
                                                        out_samp$Em + out_samp$Emv +
                                                        out_samp$Emw)
incidence = extendtp(n= length(reportable), test_prop)*reportable
incid= data.frame(date = intro_date+out_samp$time, 
                  inc = incidence, reportable = reportable)
tmp =  merge(dat, incid, by="date") %>% select(date, cases, inc, reportable) %>% 
  pivot_longer(c(2,3,4), names_to = "type", values_to = "number")
ggplot(tmp, aes(x=date, y=number, color=type))+geom_line()        



# this is pretty  good - sets the model up with reasonable immunity etc
# based on what i think I know.
# I think what I want to do here is just project forward
# with the cases shown, rather than fiddling with fitting more 


# save.image("~/Dropbox/COVID/covid-modelling-snippets/on-cases-hosps.Rdata")

# 
#
#pars_good1 = parameters; # had test_prop going down to 0.3. changing to 0.15
# pars_good2 = parameters; # had test_prop going down to 0.15 

tptest = 1 -0.85/(1+exp(-0.2*as.numeric((x - ymd("2021-12-18")))))
# now you have to run everything *again* and get parameters2 ! 
plot(x,tptest)
# test_prop2= tptest set 1 and 2 with 0.7 and 0.85 

intv_date <-  ymd("2022-03-07") 
fur_intv_date = ymd("2022-04-05")
projpar = pars_good1
projpar["relx_level"]=0.35
projpar["fur_relx_level"]=0.55
projpar["rlx_t"] = as.numeric(intv_date - intro_date)
projpar["fur_rlx_t"] = as.numeric(fur_intv_date - intro_date)
projpar2 = projpar
projpar2 = projpar + pars_good2-pars_good1 
# change - going to keep model 2, but make model 1 have a lower transmission 
projpar  = projpar2
projpar["fur_relx_level"]=0.3

end_proj = ymd("2022-07-01")
times = 1:(end_proj - intro_date)
tpproj1=extendtp(n=length(times), test_prop1)
tpproj2=extendtp(n=length(times), test_prop2)
out1 <- as.data.frame(deSolve::ode(y=init,time=times,func=sveirs,
                                       parms=projpar)) 
out2 <- as.data.frame(deSolve::ode(y=init,time=times,func=sveirs,
                                   parms=projpar2)) 

reportable1 = projpar[["p"]]*projpar[["sigma"]]*(out1$Er + out1$Erv + out1$Erw +
                                                        out1$Em + out1$Emv +
                                                        out1$Emw)
reportable2 = projpar[["p"]]*projpar[["sigma"]]*(out2$Er + out2$Erv + out2$Erw +
                                                   out2$Em + out2$Emv +
                                                   out2$Emw)

incidence1 = tpproj2*reportable1
incidence2 = tpproj2*reportable2

incid= data.frame(date = intro_date+out1$time, 
                  model1 = incidence1, reportable1 = reportable1,
                  model2=incidence2, reportable2=reportable2)
incid$modelhosp1 = lag(incid$model1, n=20)*0.55
incid$modelhosp2 = lag(incid$model2, n=20)*0.55

incid$cases = NA 

incid$cases =alldat$cases[match(incid$date, alldat$date)] # need the full data not filterd for date
incid$hosp = NA
incid$hosp = tothosp$hosp[match(incid$date, tothosp$date)]
tmp =  incid %>% 
  pivot_longer(c(2,4), names_to = "type", values_to = "number")
tmp = tmp %>% mutate(upper = make_ribbon(number, theta=0.1)$ymax,
                     lower = make_ribbon(number, theta=0.1)$ymin)
p1= ggplot(tmp, aes(x=date, y=number, color=type))+geom_line()   +
  geom_ribbon(aes(x=date, ymin=lower, ymax = upper, fill=type),
              alpha=0.3,color=NA) +
  scale_x_date(date_breaks = "months", date_labels = "%Y-%m-%d") +
  geom_point(data = filter(alldat, date>= intro_date), aes(x=date, y=cases), alpha=0.5, inherit.aes = F)+
  theme_light() +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank()) + ylab("Reported cases") +
   guides(fill=guide_legend(title="Testing"), color=guide_legend(title="Testing"))



htmp = incid %>% pivot_longer(c(6,7), names_to="type", values_to = "hosp_number")
htmp = htmp %>% mutate(upper = make_ribbon(hosp_number, theta=0.2)$ymax,
                     lower = make_ribbon(hosp_number, theta=0.2)$ymin)

p2 = ggplot(filter(htmp,date> ymd("2022-02-01")), aes(x=date, y=hosp_number, color=type))+geom_line()   +
  geom_ribbon(aes(x=date, ymin=lower, ymax = upper, fill=type),alpha=0.3,color=NA) + 
  scale_x_date(date_breaks = "months", date_labels = "%b-%d") +  # was %Y-%m-%d
  geom_point(data = filter(tothosp,date> ymd("2021-12-01")), 
             aes(x=date, y=hosp), alpha=0.5, inherit.aes = F)+
   ylab("Hospitalizations") +
  annotate("rect", xmin =  ymd("2021-12-01"), xmax = ymd("2022-02-01"), 
           ymin = 0, ymax = 3000,
           alpha = .3,fill = "grey")+theme_light()+ 
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        axis.title.x = element_blank())+
  guides(fill=guide_legend(title="Testing"), color=guide_legend(title="Testing"))

ggarrange(p1, p2, nrow = 2, common.legend = T, legend = "top")



#  scale_x_date(date_breaks = "2 weeks", date_labels = "%Y-%m-%d") 





make_ribbon = function( somedata, theta,lowp=0.1, highp=0.9) { 
  ymin = qnbinom(p=lowp, mu=somedata, size = 1/theta)
  ymax = qnbinom(p=highp, mu=somedata, size = 1/theta)
 return(data.frame(ymin = ymin, ymax = ymax))
}
tt= make_ribbon(incid$inc1, theta = 0.1)

# interpreting what i just did: how much have these things changed the infectivity
# one big question is how we could have had that stringency going on and for so long? 

with(as.list( projpar), {
c <- (1 - stngcy/(1+ exp(-1.25*(times-eff_t)))) 
rlx <- (1 + relx_level/(1+ exp(-1.25*(times-rlx_t)))) # relaxation 
further_rlx <- (1 + fur_relx_level/(1+ exp(-1.25*(times-fur_rlx_t))))
infectionfactor <- c*rlx*further_rlx
plot(infectionfactor)})

# so in the end, I want to make two great plots in two pairs (cases, hosps) 
# both with good fits 
# one with maybe lower wf or higher wf, maybe lower or higher p 
