library(CanCovidData)
ondat <- readr::read_csv("https://data.ontario.ca/dataset/f4112442-bdc8-45d2-be3c-12efae72fb27/resource/455fd63b-603d-4608-8216-7d8647f43350/download/conposcovidloc.csv")

ondat <- ondat %>%
  group_by(Case_Reported_Date) %>%
  dplyr::summarise(cases = n())
ggplot(data = ondat, aes(x = Case_Reported_Date, y = cases)) +
  geom_point(color = "blue") +
  #   geom_point(data = dat, aes(x = date, y = cases), color = "red") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%Y-%m-%d") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
# they are NOT THE SAME> interesting. need to adjust a0peak on may29, and there is
# one at april 17
# neither are as bad as the Apr 1 one so for now I'm leaving it - as of Sept 15 I am going to fit from these data
tail(ondat) # last row NA?
#ondat$cases[nrow(ondat)] <- 1525 #delete this...

dat <- ondat#[-c(nrow(ondat), nrow(ondat)), ] # FIX THIS NEXT TIME
dat <- dplyr::filter(dat, Case_Reported_Date >= lubridate::ymd("2020-02-29"))
dat$date <- dat$Case_Reported_Date
dat$date[1:2] <- c(ymd("2020-03-01"), ymd("2020-03-02"))
dat$value <- dat$cases

dat$day <- 1:nrow(dat)
dat$date[1] + max(dat$day) - 1 # I think this should be today.
diff(dat$date)
# if not, check diff() of the dates
dat <- filter(dat, date >= "2020-07-22") 
dat$day <- 1:length(dat$cases)
saveRDS(dat, file.path("data/ON-dat.rds"))
dat = readRDS( file.path("data/ON-dat.rds"))


#________________test_prop __________________

ondat <- readr::read_csv("https://data.ontario.ca/dataset/f4112442-bdc8-45d2-be3c-12efae72fb27/resource/455fd63b-603d-4608-8216-7d8647f43350/download/conposcovidloc.csv")

agedat_on <- group_by(ondat, Test_Reported_Date,Age_Group) %>%
  dplyr::summarise(cases = n()) %>%
  filter(Test_Reported_Date >= ymd("2021-09-01"))

lowerages = c("<20","20s", "30s", "40s", "50s", "60s","UNKNOWN")
agedat_on$under70 = "No" 
agedat_on$under70[which(agedat_on$Age_Group %in% lowerages)] = "Yes"

mydat_ON =  group_by(agedat_on, Test_Reported_Date , under70) %>% 
  dplyr::summarise(totcases = sum(cases)) 
mydat_ON$`Reported_Date` <- mydat_ON$Test_Reported_Date 

test_ONspline  <- make_case_splines(mydat_ON)

mytest_ON = get_testprop(changedate = ymd("2021-12-21"), 
                         mysplines = test_ONspline, 
                         halftime = 35, steepness = 0.15)
glimpse(mytest_ON)
ggplot(mytest_ON, aes(x=date, y=test_prop))+geom_line()






