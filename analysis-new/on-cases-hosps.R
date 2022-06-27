# got ON hosp from here, there is a download link
# https://data.ontario.ca/en/dataset/covid-19-cases-in-hospital-and-icu-by-ontario-health-region


onhosp = readr::read_csv("https://data.ontario.ca/dataset/8f3a449b-bde5-4631-ada6-8bd94dbc7d15/resource/e760480e-1f95-4634-a923-98161cfb02fa/download/region_hospital_icu_covid_data.csv")
#onhosp = readr::read_csv("data/region_hospital_icu_covid_data_ON.csv")

ggplot(onhosp, aes(x=date, y = hospitalizations, color=oh_region))+geom_point()

tothosp <- onhosp %>% group_by(date) %>% dplyr::summarise(hosp = sum(hospitalizations))

ggplot(tothosp, aes(x=date, y = hosp))+geom_point()



# case data 
ondat <- readr::read_csv("https://data.ontario.ca/dataset/f4112442-bdc8-45d2-be3c-12efae72fb27/resource/455fd63b-603d-4608-8216-7d8647f43350/download/conposcovidloc.csv")

oncases <- ondat %>%
  group_by(Case_Reported_Date) %>%
  dplyr::summarise(cases = n())
ggplot(data = oncases, aes(x = Case_Reported_Date, y = cases)) +
  geom_point(color = "blue") +
  #   geom_point(data = dat, aes(x = date, y = cases), color = "red") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%Y-%m-%d") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
oncases$date = oncases$Case_Reported_Date
# now just fiddle 

cashos = merge(oncases, tothosp, by = "date") %>% select(date, cases, hosp) %>% 
  mutate(scaledcases = lag(cases, n=9)*0.55) %>%
     pivot_longer(c(3,4), names_to = "type", values_to = "number")
ggplot(cashos, aes(x=date, y = number, color=type)) + geom_line() +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%Y-%m-%d") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))


# fit the model to the cases and adjust to get a projection for hosp 

# needs ON_data.R but am going to look at the spline thing 
# looks great for both pred and upred, because splines work fine
ggplot(pivot_longer(test_ONspline$upred,
                    cols = c(1,3),
                    values_to = "number",
                    names_to = "type"), 
       aes(x=Reported_Date, y=number,color=type))+geom_point() # comparison stupid bc log

# now continuing from ON_projection.R 

# setup stuff 
forecasts_days <- 60 
intro_date <-  ymd("2021-11-30")
stop_date <- ymd("2022-02-12")   # I don't know why we would stop on this date but I will 
# come back to this later 
#import data 
# dat = readRDS("data/ON-dat.rds")
#include Omicron wave only
dat <- dat %>% filter(date >= intro_date &  date <= stop_date)
dat_omic <- dat
dat_omic <- filter(dat_omic, date >= intro_date) %>% dplyr::select(c("day", "value", "date"))
dat_omic$day <- 1:nrow(dat_omic)




########## Test_prop set up
#-------EB: test_prop now has dates and this ensures it starts at the right date 
# note - the fit at line 30 with tp_approx looks terrible so I need to use the 
# mytest_ON one explicitly, instead 

# Elisha had: test_prop_proj_ON <- data.frame(tp_approx[2]) # should have date and tp and 
# should start at the start date and end at in this case April 13?) 
# I will use: 
test_prop_proj_ON <- mytest_ON %>% select(date, test_prop) %>% dplyr::rename(tp=test_prop)
test_prop <- test_prop_proj_ON$tp[1:length(dat_omic$day)]

# am going to do something even dumber. test_prop from 1 to about 0.3 over dec 15-about 25
x=intro_date+1:nrow(dat_omic)
tptest = 1 -0.85/(1+exp(-0.2*as.numeric((x - ymd("2021-12-18")))))
plot(x,tptest)
test_prop= tptest


## parameter and model setup
sum(filter(cashos, date< ymd("2021-11-30"))$cases) # 1,230,372 reported cases or 8.4% of ON

N=14.57e6
N_pop=N
vaxlevel_in =  0.75 # 
simu_size = 10000 # number of times to resample the negative binom (for ribbons)
forecasts_days =60 # how long to forecast for 
times = 1:nrow(dat_omic)

#declaring  parameters 
eff_date <-   ymd("2021-12-30")  # intervention date # restrictions were updated on 24 December 
intv_date <-  ymd("2022-02-10") 
fur_intv_date <- ymd("2022-05-10") 
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
                        N=14.57e6,
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
incres_in = 738 #866#(made up factor) # resident strain (delta) incidence at start (938, reported cases) 
incmut_in = 150 # new (omicron) inc at stat (presumably very low)

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
pen.fit_ON <- optim(fn=func_penloglik,  par=guess, lower=c(0,0,0,0.001), upper = c(Inf,1,Inf,Inf), 
                    method = "L-BFGS-B", 
                    parameters = parameters, test_prop=test_prop, 
                    pen.size=pen.size, penalties = penalties, dat_omic=dat_omic, hessian=T)
pen.fit_ON

# check the fit 
func_loglik(pen.fit_ON$par, test_prop, dat_omic,parameters) 
parameters[names(guess)] <- pen.fit_ON$par


gg = simple_prev_plot(parameters, numdays = 190, mode = "both"); gg  # cc's simple prevalence plot 


# i want to plot the incidence against reported cases. 
# for this i need to get the incidence out, and adjust by test_prop, and plot 
# that against reported cases 

out_samp <- as.data.frame(deSolve::ode(y=init,time=times,func=sveirs,
                                       parms=parameters)) 

reportable = parameters[["p"]]*parameters[["sigma"]]*(out_samp$Er + out_samp$Erv + out_samp$Erw +
                                                        out_samp$Em + out_samp$Emv +
                                                        out_samp$Emw)
incidence = test_prop*reportable #(length differs)

incid = data.frame(date = intro_date+out_samp$time, 
                  inc = incidence, reportable = reportable)
tmp =  merge(dat, incid, by="date") %>% select(date, cases, inc, reportable) %>% 
  pivot_longer(c(2,3,4), names_to = "type", values_to = "number")
ggplot(tmp, aes(x=date, y=number, color=type))+geom_line()        



# this is pretty  good - sets the model up with reasonable immunity etc
# based on what i think I know.
# I think what I want to do here is just project forward
# with the cases shown, rather than fiddling with fitting more 


# save.image("~/Dropbox/COVID/covid-modelling-snippets/on-cases-hosps.Rdata")

# pars_good1 = parameters; # had test_prop going down to 0.3. changing to 0.15
# pars_good2 = parameters; # had test_prop going down to 0.15 

tptest = 1 -0.7/(1+exp(-0.2*as.numeric((x - ymd("2021-12-18")))))
plot(x,tptest)

#test_prop1= tptest #set 1 and 2 with 0.7 and 0.85 

tptest = 1 -0.85/(1+exp(-0.2*as.numeric((x - ymd("2021-12-18")))))
plot(x,tptest)

test_prop2= tptest #set 1 and 2 with 0.7 and 0.85 

intv_date <-  ymd("2022-03-07") 
fur_intv_date = ymd("2022-04-05")
projpar = pars_good1
projpar["relx_level"]=0.7
projpar["fur_relx_level"]=0.4
projpar["rlx_t"] = as.numeric(intv_date - intro_date)
projpar["fur_rlx_t"] = as.numeric(fur_intv_date - intro_date)
projpar2 = projpar
projpar2 = projpar + pars_good2-pars_good1

end_proj = ymd("2022-07-01")
times = 1:(end_proj - intro_date)
tpproj1=extendtp(n=length(times), test_prop1)  #test_prop1 not defined 
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

incidence1 = tpproj1*reportable1
incidence2 = tpproj2*reportable2

incid= data.frame(date = intro_date+out2$time-1,#intro_date+out_samp$time, 
                  model1 = incidence1, reportable1 = reportable1,
                  model2=incidence2, reportable2=reportable2)
incid$modelhosp1 = lag(incid$model1, n=9)*0.55
incid$modelhosp2 = lag(incid$model2, n=9)*0.55

incid$cases = NA 
incid$cases = oncases$cases[match(incid$date, oncases$date)]
incid$hosp = NA
incid$hosp = tothosp$hosp[match(incid$date, tothosp$date)]
tmp =  incid %>% 
  pivot_longer(c(2,4), names_to = "type", values_to = "number")
tmp = tmp %>% mutate(upper = make_ribbon(number, theta=0.1)$ymax,
                     lower = make_ribbon(number, theta=0.1)$ymin)
p1 = ggplot(tmp, aes(x=date, y=number, color=type))+geom_line()   +
  geom_ribbon(aes(x=date, ymin=lower, ymax = upper, fill=type),alpha=0.3,color=NA) + 
  scale_x_date(date_breaks = "months", date_labels = "%Y-%m-%d") +
  geom_point(data = incid, aes(x=date, y=cases), alpha=0.5, inherit.aes = F)+
  theme_light() +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank()) + ylab("Reported cases") +
   guides(fill=guide_legend(title="Testing"), color=guide_legend(title="Testing"))



htmp = incid %>% pivot_longer(c(8,9), names_to="type", values_to = "hosp_number")
htmp = htmp %>% mutate(upper = make_ribbon(hosp_number, theta=0.03)$ymax,
                     lower = make_ribbon(hosp_number, theta=0.03)$ymin)

p2 = ggplot(filter(htmp,date> ymd("2022-01-20")), aes(x=date, y=hosp_number, color=type))+geom_line()   +
  geom_ribbon(aes(x=date, ymin=lower, ymax = upper, fill=type),alpha=0.3,color=NA) + 
  scale_x_date(date_breaks = "months", date_labels = "%b-%d") +  # was %Y-%m-%d
  geom_point(data = filter(tothosp,date> ymd("2021-12-01")), 
             aes(x=date, y=hosp), alpha=0.5, inherit.aes = F)+
   ylab("Hospitalizations") +
  annotate("rect", xmin =  ymd("2021-12-01"), xmax = ymd("2022-01-22"), 
           ymin = 0, ymax = 12500,
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
