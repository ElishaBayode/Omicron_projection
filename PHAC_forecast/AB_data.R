library(CanCovidData)
library(lubridate)
library(dplyr)

#----------case data
abdat <- get_alberta_case_data()
abdat$thiscase <- 1
abdat <- abdat %>% group_by(`Date reported`) %>%
  dplyr::summarise(cases = sum(thiscase)) %>% dplyr::rename(date = `Date reported`)

ggplot(data = abdat, aes(x = date, y = cases)) +
  geom_point() +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%Y-%m-%d") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid.major = element_line(color = "grey"))
startDate <- lubridate::ymd("2020-03-06")
dat <- filter(abdat, date > startDate)

tail(dat)

dat$value <- dat$cases
#dat$adjust_cases <-dat$cases 
dat$day <- 1:nrow(dat) # may want to check the diff(date) again each time
diff(dat$date)
dat$day <- seq_len(nrow(dat)) # reset day again? ?? 
dat
#tail(dat,20)


dat <- filter(dat, date >= "2020-06-22")

#tail(dat)
dat$day <- 1:length(dat$cases)
saveRDS(dat, file.path("data/AB-dat.rds"))


#------------- get test_prop  

dat = get_alberta_case_data()

# going to start this whole thin in fall 2021 
agedat <- group_by(dat, `Date reported`, `Age group`) %>%
  dplyr::summarise(cases = n()) %>%
  filter(`Date reported` >= ymd("2021-09-01"))

# get a time series for <70 and 70+ 
#lowerages = c("<10" , "10-19","20-29", "30-39",  "40-49" , "50-59" , "60-69")
lowerages = c("Under 1 year", "1-4 years" , "5-9 years","10-19 years", 
              "20-29 years",  "30-39 years",  "40-49 years" ,"50-59 years", "60-69 years",
              "Unknown")

agedat$under70 = "No" 
agedat$under70[which(agedat$`Age group` %in% lowerages)] = "Yes"

mydat_AB =  group_by(agedat,  `Date reported`, under70) %>% 
  dplyr::summarise(totcases = sum(cases)) 
mydat_AB$`Reported_Date` <- mydat_AB$`Date reported`

test_ABspline  <- make_case_splines(mydat_AB)

mytest_AB = get_testprop(changedate = ymd("2021-12-21"), 
                         mysplines = test_ABspline, 
                         halftime = 30, steepness = 0.1)
glimpse(mytest_AB)
ggplot(mytest_AB, aes(x=date, y=test_prop))+geom_line()



