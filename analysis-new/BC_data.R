library(CanCovidData)
library(lubridate)
library(tidyverse)
library(tidyr)
library(dplyr)
source("analysis-new/mod_fitting_setup.R")

-------#gets updated data 

bcpub <-get_british_columbia_case_data() #readr::read_csv("data-raw/BCCDC_COVID19_Dashboard_Case_Details.csv")
#bcpub <- readr::read_rds("bcpub-2020-09-29.rds")

bcpub$Reported_Date <- lubridate::ymd(bcpub$`Reported Date`)
bcpub$thiscase <- 1
dat <- group_by(bcpub, Reported_Date) %>%
  dplyr::summarise(cases = sum(thiscase)) %>%
  filter(Reported_Date >= ymd("2020-02-27"))
dat <- dat[order(dat$Reported_Date), ]
dat$day <- seq(1, nrow(dat))
dat$value <- dat$cases

dat$Reported_Date[1:2] <- c(ymd("2020-03-01"), ymd("2020-03-02"))

dat$day <- seq(1, nrow(dat))
dat$value <- dat$cases
dat$date <- dat$Reported_Date
#dat <- dat[-nrow(dat), ] # remove today - incomplete reporting before 3pm and other fiddles
tail(dat)

ggplot(dat, aes(Reported_Date, cases)) +
geom_point() +
scale_x_date(date_breaks = "2 weeks", date_labels = "%Y-%m-%d") +
theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid.major = element_line(color = "grey"))


dat_add <- filter(dat, date >= "2020-06-22")
dat_add$day <- 1:length(dat_add$cases)

dat <- filter(dat, date >= "2020-06-22")
dat$day <- 1:length(dat$cases)
saveRDS(dat, file.path("data/BC-dat.rds"))





###### make test_prop (testing correction) ################

dat = get_british_columbia_case_data()

# going to start this whole thin in fall 2021 
agedat <- group_by(dat,`Reported Date`, `Age group`) %>%
  dplyr::summarise(cases = n()) %>%
  filter(`Reported Date` >= ymd("2021-09-01"))

# get a time series for <70 and 70+ 
lowerages = c("<10" , "10-19","20-29", "30-39",  "40-49" , "50-59" , "60-69")
agedat$under70 = "No" 
agedat$under70[which(agedat$`Age group` %in% lowerages)] = "Yes"

mydat_BC =  group_by(agedat,`Reported Date`, under70) %>% 
  dplyr::summarise(totcases = sum(cases)) 

mydat_BC$Reported_Date <- mydat_BC$`Reported Date`


test_BCspline  <- make_case_splines(mydat_BC)

mytest_BC = get_testprop(changedate = ymd("2021-12-21"), 
                         mysplines = test_BCspline, 
                         halftime = 60, steepness = 0.05)
glimpse(mytest_BC)
ggplot(mytest_BC, aes(x=date, y=test_prop)) + geom_line()


-------#EB: fits test_prop and add dates   
tp_approx_fit = function(mytest_BC, forecasts_days=forecasts_days, howlow = howlow, 
   slope = slope, dat= dat_omic, midpoint=midpoint, intro_date= intro_date, stop_date=stop_date){
   
  dat <- dat %>% filter(date >= intro_date &  date <= stop_date)
  mytest_BC <- mytest_BC %>% filter(date >= intro_date &  date <= stop_date)
   
   fake_tp = data.frame(tp=(1 - (1-howlow)/(1 + exp(-slope*(1:(nrow(dat)+forecasts_days)-midpoint)))))
   fake_tp  = fake_tp %>% mutate(date=seq.Date(ymd(intro_date), ymd(intro_date)-1+(nrow(dat)+forecasts_days), 1))
   
   return(list(ggplot() + geom_line(data= fake_tp, aes(x=date, y=tp, colour="fit")) +
   geom_point(data= mytest_BC, aes(x=date, y=test_prop, colour="tp")) + 
     scale_x_date(date_breaks = "20 days", date_labels = "%b-%d-%y"),fake_tp))
  }
          



















