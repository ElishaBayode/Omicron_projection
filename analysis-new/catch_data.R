library(CanCovidData)
library(lubridate)
library(tidyverse)
library(tidyr)
library(dplyr)
location <- "CAN.csv"
linkRaw <- "https://raw.githubusercontent.com/ishaberry/Covid19Canada/master/timeseries_prov/cases_timeseries_prov.csv"
data <- readr::read_csv(linkRaw)
readr::write_csv(data, file.path("data", location))


##################AB ##########################
#download .csv manually and save in data
# https://www.alberta.ca/stats/covid-19-alberta-statistics.htm

#abdat <- readr::read_csv("data/covid-19-alberta-statistics-data.csv")
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


################### BC ###############################

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

#saveRDS(dat_add, file.path("data/BC-dat_add.rds")) # was PUB-dat but moving to this
#dat_add = readRDS( file.path("data/BC-dat_add.rds"))


dat <- filter(dat, date >= "2020-06-22")
dat$day <- 1:length(dat$cases)
saveRDS(dat, file.path("data/BC-dat.rds"))


#################### MB ############################



dat <- readr::read_csv(file.path("data/CAN.csv"))
dat$date <- lubridate::dmy(dat$date_report)
dat <- dplyr::filter(dat, province == "Manitoba")
# View(dat)
ggplot(dat, aes(date, cases)) +
  geom_point()
# Pick a reasonable starting date:
dat <- dplyr::filter(dat, date >= lubridate::ymd("2020-03-10"))
dat$day <- seq_len(nrow(dat))
ggplot(dat, aes(date, cases)) +
  geom_point() +
  scale_x_continuous(breaks = seq(ymd("2020-01-01"), ymd("2020-10-01"), by = 30))

tail(dat)
dat$value <- dat$cases

#dat <- dat[-c(nrow(dat), nrow(dat)), ]






realloc <- function(total, numdays) {
  fs <- runif(numdays, min = (1 / numdays) * 0.75, max = (1 / numdays) * 1.25)
  return(c(round(total * fs)))
}



# fix 26 Dec 1:
dat$adjust_cases <- dat$cases
ii_jul1 <- which(dat$date == ymd("2020-12-26"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)



dat$adjust_cases <- dat$cases
ii_jul1 <- which(dat$date == ymd("2020-12-25"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)



dat$adjust_cases <- dat$cases
ii_jul1 <- which(dat$date == ymd("2021-01-01"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)




ggplot(dat, aes(date, adjust_cases)) +
  geom_point()

# fix august 1 and any other long weeekends
isStartLongWeekend <- (dat$adjust_cases == 0 &
                         dat$adjust_cases[c(2:nrow(dat), 1)] == 0 &
                         dat$adjust_cases[c(3:nrow(dat), 1:2)] == 0
                       & dat$date > ymd("2020-12-24"))
# ii = which(dat$date==ymd("2020-08-01"))
for (ii in which(isStartLongWeekend)) {
  dat$adjust_cases[ii:(ii + 3)] <- realloc(total = dat$cases[ii + 3], numdays = 4)
}
# round(0.25*dat$adjust_cases[ii+3])


# fix other weekends
isStartWeekend <- (dat$adjust_cases == 0 &
                     dat$adjust_cases[c(2:nrow(dat), 1)] == 0
                   & dat$date > ymd("2020-12-24"))
for (ii in which(isStartWeekend)) {
  dat$adjust_cases[ii:(ii + 2)] <- realloc(total = dat$adjust_cases[ii + 2], numdays = 3)
}


isStartWeekend <- (dat$adjust_cases == 0 &
                     dat$adjust_cases[c(2:nrow(dat), 1)] == 0
                   & dat$date > ymd("2021-02-16"))
for (ii in which(isStartWeekend)) {
  dat$adjust_cases[ii:(ii + 2)] <- realloc(total = dat$adjust_cases[ii + 2], numdays = 3)
}



ii_jul1 <- which(dat$date == ymd("2021-02-15"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)



ii_jul1 <- which(dat$date == ymd("2021-04-04"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)

ii_jul1 <- which(dat$date == ymd("2021-04-02"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)

ii_jul1 <- which(dat$date == ymd("2021-07-01"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)

ii_jul1 <- which(dat$date == ymd("2021-09-30"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)

ii_jul1 <- which(dat$date == ymd("2021-11-11"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)



ggplot(dat, aes(date, adjust_cases)) +
  geom_point() +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%Y-%m-%d") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid.major = element_line(color = "grey"))
dat$value <- dat$adjust_cases


dat <- filter(dat, date >= "2020-06-22")

dat$day <- 1:length(dat$cases)
saveRDS(dat, file.path("data/MB-dat.rds"))
dat =readRDS(file.path("data/MB-dat.rds"))

tail(dat)
######################## ON ######################


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

dat <- ondat[-c(nrow(ondat), nrow(ondat)), ] # FIX THIS NEXT TIME
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

######################## QC ###############################


dat <- readr::read_csv(file.path("data/CAN.csv"))
dat$date <- lubridate::dmy(dat$date_report)
dat <- dplyr::filter(dat, date >= lubridate::ymd("2020-03-14"))
qc_dat <- dplyr::filter(dat, province == "Quebec")
qc_dat$day <- seq_len(nrow(qc_dat))
qc_dat$cases




ggplot(qc_dat, aes(date, cases)) +
  geom_point() +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%Y-%m-%d") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid.major = element_line(color = "grey"))

# Quebec public health officials announced Sunday that a computer error resulted
# in 1,317 missing positive COVID-19 cases between April 2-30.
# https://montrealgazette.com/news/latest-covid-19-statistics-for-quebec-include-1317-previously-unreported-cases-for-april/
outlier <- which(qc_dat$cases > 2000 & qc_dat$day < 65)
qc_dat$date[outlier]
qc_dat$cases[outlier]
qc_dat$cases[outlier] <- qc_dat$cases[outlier] - 1317

dates_with_missing <- qc_dat$date >= ymd("2020-04-02") & qc_dat$date <= ymd("2020-04-30")

qc_dat$cases_adjusted <- qc_dat$cases
qc_dat$cases_adjusted[dates_with_missing] <- qc_dat$cases[dates_with_missing] +
  1317 / length(dates_with_missing)
qc_dat$cases_adjusted <- round(qc_dat$cases_adjusted)
qc_dat$value <- qc_dat$cases_adjusted # for plotting function

stopifnot(identical(qc_dat$value[105:107], c(0, 0, 0)))
fill_in <- qc_dat$value[108] / 4
fill_in <- c(rep(ceiling(fill_in), 3), floor(fill_in))
stopifnot(sum(fill_in) == qc_dat$value[108])
qc_dat$value[105:108] <- fill_in



stopifnot(identical(qc_dat$value[118:120], c(0, 0, 0)))
fill_in <- qc_dat$value[121] / 4
fill_in <- c(rep(ceiling(fill_in), 3), floor(fill_in))
stopifnot(sum(fill_in) == qc_dat$value[121])
qc_dat$value[118:121] <- fill_in
plot(qc_dat$value)

ggplot(qc_dat, aes(date,value)) +
  geom_point() +
  geom_line(aes(y = value), colour = "blue") +
  geom_vline(xintercept = ymd("2020-03-23")) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%Y-%m-%d") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid.major = element_line(color = "grey"))






dat = qc_dat

dat <- dat[-c(nrow(dat)+ 1, nrow(dat)), ]
# dat <- filter(qc_dat, date > ymd("2020-03-14"))
# dat$day <- seq_len(nrow(dat)) # trying this 



realloc <- function(total, numdays) {
  fs <- runif(numdays, min = (1 / numdays) * 0.75, max = (1 / numdays) * 1.25)
  return(c(round(total * fs)))
}



# fix 26 Dec 1:
dat$adjust_cases <- dat$cases
ii_jul1 <- which(dat$date == ymd("2020-12-26"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)

#dat$adjust_cases <- dat$cases
ii_jul1 <- which(dat$date == ymd("2020-12-25"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)

ii_jul1 <- which(dat$date == ymd("2021-07-01"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)

ii_jul1 <- which(dat$date == ymd("2021-10-15"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)


ii_jul1 <- which(dat$date == ymd("2022-01-04"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)

# fix other weekends
isStartWeekend <- (dat$adjust_cases == 0 &
                     dat$adjust_cases[c(2:nrow(dat), 1)] == 0
                   & dat$date > ymd("2020-12-30"))
for (ii in which(isStartWeekend)) {
  dat$adjust_cases[ii:(ii + 2)] <- realloc(total = dat$adjust_cases[ii + 2], numdays = 3)
}

dat <- as.data.frame(dat)

length(dat$adjust_cases)

isStartWeekend <- (dat$adjust_cases == 0 &
                     dat$adjust_cases[c(2:nrow(dat), 1)] == 0
                   & dat$date > ymd("2020-07-10"))
for (ii in which(isStartWeekend)) {
  dat$adjust_cases[ii:(ii + 2)] <- realloc(total = dat$adjust_cases[ii + 2], numdays = 3)
}



ggplot(dat, aes(date, adjust_cases)) +
  geom_point() +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%Y-%m-%d") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid.major = element_line(color = "grey"))


dat$value <- abs(dat$adjust_cases)

tail(dat,10)


ggplot(dat, aes(date, adjust_cases)) +
  geom_point()


dat <- filter(dat, date >= "2020-08-01")

dat$day <- 1:length(dat$cases)
saveRDS(dat, file.path("data/QC-dat.rds"))
dat =readRDS(file.path("data/QC-dat.rds"))


############ SK ##########################




# Notes ---------------------------------------------------------------------



# Read and prepare data -----------------------------------------------------

#dat <- readr::read_csv("https://dashboard.saskatchewan.ca/export/cases/2097.csv")
 dat <- readr::read_csv("https://raw.githubusercontent.com/ishaberry/Covid19Canada/master/timeseries_prov/cases_timeseries_prov.csv")

dat <- readr::read_csv(file.path("data/CAN.csv"))
dat$date <- lubridate::dmy(dat$date_report)
dat <- dplyr::filter(dat, province == "Saskatchewan")
# View(dat)
ggplot(dat, aes(date, cases)) +
  geom_point() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b")
# Pick a reasonable starting date - for SK let's ignore the early outbreaks entirely and start on June 1
startDate <- lubridate::ymd("2020-06-10")
dat <- dplyr::filter(dat, date >= startDate) # had been 1st but ..
dat$value <- dat$cases
dat$day <- seq_len(nrow(dat))
ggplot(dat, aes(date, cases)) +
  geom_point()

# ##### public data (more current)###############
#dat <- readr::read_csv("https://dashboard.saskatchewan.ca/export/cases/2097.csv")
# dat <- group_by(dat, Date) %>%
#   summarise(cases = sum(`New Cases`)) %>%
#   filter(Date >= ymd("2020-06-10"))
# dat <- dat[order(dat$Date), ]
# dat$date <- dat$Date
# 
# ggplot(dat, aes(date, cases)) +
#   geom_point()
# 
# ###########################
realloc <- function(total, numdays) {
  fs <- runif(numdays, min = (1 / numdays) * 0.75, max = (1 / numdays) * 1.25)
  return(c(round(total * fs)))
}



# fix 26 Dec :


dat$adjust_cases <- dat$cases
ii_jul1 <- which(dat$date == ymd("2020-12-26"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)



#dat$adjust_cases <- dat$cases
ii_jul1 <- which(dat$date == ymd("2020-12-25"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)


ii_jul1 <- which(dat$date == ymd("2020-12-28"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)

ii_jul1 <- which(dat$date == ymd("2021-01-01"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)

ii_jul1 <- which(dat$date == ymd("2021-01-10"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)

ii_jul1 <- which(dat$date == ymd("2021-02-15"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)

ii_jul1 <- which(dat$date == ymd("2021-12-27"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)

ii_jul1 <- which(dat$date == ymd("2021-12-26"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)

ii_jul1 <- which(dat$date == ymd("2021-12-25"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)

ii_jul1 <- which(dat$date == ymd("2021-12-24"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)


###temporary: check next week if this is  needed
ii_jul1 <- which(dat$date == ymd("2022-01-03"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)

ii_jul1 <- which(dat$date == ymd("2022-01-02"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)

ii_jul1 <- which(dat$date == ymd("2022-01-01"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)



ggplot(dat, aes(date, adjust_cases)) +
  geom_point() +
  # x <- seq(0, 40, length.out = 200)
  # plot(x, dlnorm(x, log(12), 0.1), type = "l", xaxs = "i", yaxs = "i")
  scale_x_date(date_breaks = "2 weeks", date_labels = "%Y-%m-%d") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid.major = element_line(color = "grey"))


tail(dat,20)
#dat$value <- dat$cases
dat$value <- dat$adjust_cases
dat <- dat[-c(nrow(dat)-1, nrow(dat)), ]

dat <- filter(dat, date >= "2020-09-22")
#dat <- filter(dat, date <= "2021-12-23")
dat$day <- 1:length(dat$cases)
saveRDS(dat, file.path("data/SK-dat.rds"))
dat =readRDS(file.path("data/SK-dat.rds"))






###### AB_age structured data ############


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


################################  QC test Prop  ###########################

dat_qc <- readr::read_csv("data/DomesticSurveillanceData_2022-01-13_DISCOVER.csv")

dat_qc <- filter(dat_qc,  reporteddate >= ymd("2021-09-01")) 
# going to start this whole thin in fall 2021 
agedat_qc <- group_by(dat_qc, reporteddate,agegroup10) %>%
  dplyr::summarise(cases = n()) %>%
  filter(reporteddate >= ymd("2021-09-01"))
lowerages = c("0 to 19" ,"20 to 29", "30 to 39", "40 to 49", "50 to 59", "60 to 69", "unknown")
agedat_qc$under70 = "No" 
agedat_qc$under70[which(agedat_qc$agegroup10 %in% lowerages)] = "Yes"

mydat_QC =  group_by(agedat_qc, reporteddate , under70) %>% 
  dplyr::summarise(totcases = sum(cases)) 
mydat_QC$`Reported_Date` <- mydat_QC$reporteddate

test_QCspline  <- make_case_splines(mydat_QC)

mytest_QC = get_testprop(changedate = ymd("2021-12-21"), 
                         mysplines = test_QCspline, 
                         halftime = 35, steepness = 0.15)
glimpse(mytest_QC)
ggplot(mytest_QC, aes(x=date, y=test_prop))+geom_line()


################################  ON test Prop  ###########################
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



######################BC ##################################################

dat = get_british_columbia_case_data()

# going to start this whole thin in fall 2021 
agedat <- group_by(dat,`Reported Date`, `Age group`) %>%
  dplyr::summarise(cases = n()) %>%
  filter(`Reported Date` >= ymd("2021-09-01"))

# get a time series for <70 and 70+ 
lowerages = c("<10" , "10-19","20-29", "30-39",  "40-49" , "50-59" , "60-69")
agedat$under70 = "No" 
agedat$under70[which(agedat$`Age group` %in% lowerages)] = "Yes"

mydat =  group_by(agedat,`Reported Date`, under70) %>% 
  dplyr::summarise(totcases = sum(cases)) 

mydat$Reported_Date <- mydat$`Reported Date`

# make the splines 
splinetest = make_case_splines(mydat)

# check the quality 
ggplot(data = splinetest$pred, aes(x=Reported_Date, y=lcases))+geom_point(color="blue", alpha=0.5) +
  geom_line(data =splinetest$pred, aes(x=Reported_Date, y=predlcases))

# make the test proportion by adjusting from 70+ 

mytest = get_testprop(changedate = ymd("2021-12-21"), 
                      mysplines = splinetest, 
                      halftime = 20, steepness = 0.1)
# check : should be 1, then less than 1
ggplot(mytest, aes(x=date, y=test_prop))+geom_line()


### ######### -- begin junk leftovers ######## ########
# -- note - a lot of this is now available with the new functions
# make_case_splines and get_testprop. 
# the code below is a messy exploration done in development of those two functions


over70 = filter(mydat, under70=="No")
logover70 = over70 %>% mutate(lcases = log(totcases)) %>% select(Reported_Date, lcases)
ospline = smooth.spline(logover70$Reported_Date, logover70$lcases, df=15)
pred = data.frame(x = ospline$x, predlcases = ospline$y, 
                  Reported_Date = logover70$Reported_Date)
ggplot(data = logover70, aes(x=Reported_Date, y=lcases))+geom_point(color="blue", alpha=0.5) +
  geom_line(data =pred, aes(x=Reported_Date, y=predlcases))
# (2) under 70 
under70 = filter(mydat, under70=="Yes")
logunder70 = under70 %>% mutate(lcases = log(totcases)) %>% select(Reported_Date, lcases)
uspline = smooth.spline(logunder70$Reported_Date, logunder70$lcases, df=15)
upred = data.frame(x = uspline$x, predlcases = uspline$y, 
                   Reported_Date = logunder70$Reported_Date)
ggplot(data = logunder70, aes(x=Reported_Date, y=lcases))+geom_point(color="blue", alpha=0.5) +
  geom_line(data =upred, aes(x=Reported_Date, y=predlcases))

# three options for the offset 
# get the offset between 70+ and under 70 cases at a particular date, 
# which is the difference between the two above splines at a fixed date 
mydate = ymd("2021-12-21") 
# (option 1) offset at a particular date (Sally) 
offset = filter(upred,Reported_Date == mydate)$predlcases -
  filter(pred, Reported_Date == mydate)$predlcases 
# (option 2) offset time series - this varies quite a bit, actually 
seriesoffset = data.frame(date = upred$Reported_Date, soffset= upred$predlcases - 
                            pred$predlcases) 
ggplot(data = seriesoffset, aes(x=date, y = soffset))+geom_point()+
  geom_line() # yes , dec 21 was an absolute peak 

# (option 3) mean offset over the time in the fall before my chosen date 
offset = mean(filter(upred,Reported_Date <= mydate)$predlcases - 
                filter(pred,Reported_Date <= mydate)$predlcases)
# the problem is that this produces a model fit that is *under* the real reported cases in 
# those under 70 (because in Dec the offset was higher than the mean for the fall) 

# add whichever offset you have chosen to the 70+ cases to get a model for the total cases. 
pred$totmodel = pred$predlcases + offset 

# plot this with sensible y labels. note - would be better to make 
# one or two long data frames, and then there would be a legend that would make sense 
df = data.frame(date = logunder70$Reported_Date, 
                under70cases = filter(mydat, under70=="Yes")$totcases, 
                under70spline = exp(upred$predlcases), 
                over70cases = filter(mydat, under70 =="No")$totcases, 
                over70spline = exp(pred$predlcases), 
                totalsplines = exp(pred$totmodel)) 
ggplot(data = df, aes(x=date, y=under70cases))+
  geom_point(color="blue", alpha=0.5) +
  geom_line(data = df, aes(x=date, y=under70spline),color="blue") + 
  geom_point(data = df, aes(x=date, y =over70cases),color="darkgreen",alpha=0.5) + 
  geom_line(data = df, aes(x=date, y=over70spline),color="darkgreen") +
  geom_line(data=df, aes(x=date, y=totalsplines), color="black")+
  scale_y_continuous(trans = "log10")+annotation_logticks() +
  ylab("Cases - log scale") + ggtitle("Green - 70+. Blue - under 70. Black: combined spline") 



# ---- making test_prop using this stuff ---- 
# i want a function that transitions from the offset on dec 21 to the mean offset over the 
# fall period, in a time frame of about a month.
thalf=15
myoffset <- 3.73 - (3.73-2.68)/(1+ exp(-0.25*(1:60-thalf)))
plot(1:60, myoffset)

getoffset = function(startvalue = 3.73, endvalue=2.68, halftime=15, steepness=0.25, ndays=60) {
  return( startvalue - (startvalue-endvalue)/(1+exp(-steepness*(1:ndays-halftime))))  
}
myoffset= getoffset(halftime=100,steepness = 0.05) # really flat, like Sally's

# (option 4) 
L1 = ymd("2021-12-21")-min(pred$Reported_Date) # 
pred$totmodel = NA
pred$totmodel[1:L1] = upred$predlcases[1:L1] -pred$predlcases[1:L1]
L2 = nrow(pred)- L1 # how many values to fill in with the new offset? 
pred$totmodel[(L1+1):nrow(seriesoffset)] = myoffset[1:L2]



# pred$totmodel = pred$predlcases + seriesoffset$model
df = data.frame(date = logunder70$Reported_Date, 
                under70cases = filter(mydat, under70=="Yes")$totcases, 
                under70spline = exp(upred$predlcases), 
                over70cases = filter(mydat, under70 =="No")$totcases, 
                over70spline = exp(pred$predlcases), 
                totalsplines = exp(pred$totmodel)) 
ggplot(data = df, aes(x=date, y=under70cases))+
  geom_point(color="blue", alpha=0.5) +
  geom_line(data = df, aes(x=date, y=under70spline),color="blue") + 
  geom_point(data = df, aes(x=date, y =over70cases),color="darkgreen",alpha=0.5) + 
  geom_line(data = df, aes(x=date, y=over70spline),color="darkgreen") +
  geom_line(data=df, aes(x=date, y=totalsplines), color="black")+
  scale_y_continuous(trans = "log10")+annotation_logticks() +
  ylab("Cases - log scale") + ggtitle("Green - 70+. Blue - under 70. Black: combined spline") 

# --- here is how to get a smooth nice test_prop 
# we use  (model for Cases in 70+ )* ( 1 + exp(offset) ) for the 'would have beens'
# test_prop is now (spline model for total cases) / would have beens 
off1= getoffset(halftime=100,steepness = 0.05) # this one is really flat, like Sally's
offconst = upred$predlcases-pred$predlcases
L1 = ymd("2021-12-21")-min(pred$Reported_Date) # 
L2 = nrow(pred)- L1 # how many values to fill in with the new offset? 
offconst = upred$predlcases-pred$predlcases
offconst[(L1+1):nrow(seriesoffset)] = off1[1:L2]

testpropdf = data.frame(date = pred$Reported_Date, 
                        totalmodel = exp(pred$predlcases)+ exp(upred$predlcases), 
                        totaladjusted = exp(pred$predlcases+offconst), 
                        test_prop = pmin(1,(exp(pred$predlcases)+ exp(upred$predlcases))/ exp(pred$predlcases+offconst)))
plot(testpropdf$date, testpropdf$test_prop)

################################  SK test Prop  ###########################


dat = get_saskatchewan_case_data()

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



