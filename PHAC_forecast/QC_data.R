library(lubridate)
library(tidyverse)
library(tidyr)
library(dplyr)
location <- "CAN.csv"
linkRaw <- "https://raw.githubusercontent.com/ishaberry/Covid19Canada/master/timeseries_prov/cases_timeseries_prov.csv"
data <- readr::read_csv(linkRaw)
readr::write_csv(data, file.path("data", location))



#soureces of hospital data: https://resources-covid19canada.hub.arcgis.com/datasets/provincial-daily-totals/explore  

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

ii_jul1 <- which(dat$date == ymd("2022-02-11"))
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


###### make test_prop (testing correction) ################ 
#download data manually and save to data folder 
#https://www.inspq.qc.ca/covid-19/donnees/age-sexe/evolution-cas

dat_qc <- readr::read_csv("data/graph_1-1.a_page_age_et_sexe_evol_des_cas.csv")

dat_qc <- filter(dat_qc,  `Date de déclaration` >= ymd("2021-09-01")) 

# going to start this whole thin in fall 2021 
#agedat_qc <- group_by(dat_qc, `Date de déclaration`,agegroup10) %>%
 # dplyr::summarise(cases = n()) %>%
  #filter(reporteddate >= ymd("2021-09-01"))

totat_under70 <- dat_qc$`0-9 ans`+ dat_qc$`10-19 ans`+ dat_qc$`20-29 ans` + dat_qc$`30-39 ans` + dat_qc$`40-49 ans`+
               dat_qc$`50-59 ans` + dat_qc$`60-69 ans`

totalover70 <- dat_qc$`70-79 ans` + dat_qc$`80-89 ans`+ dat_qc$`50-59 ans` + dat_qc$`90 ans et plus`
                

dat_qc$Reported_Date <-  dat_qc$`Date de déclaration`
dat_qc$totat_under70 <- totat_under70
dat_qc$totalover70   <-   totalover70

dat_qc_pt1 <- dplyr::select(dat_qc, c("Reported_Date", "totat_under70"))
dat_qc_pt2 <- dplyr::select(dat_qc, c("Reported_Date", "totalover70"))

dat_qc_pt1$totcases <- dat_qc_pt1$totat_under70
dat_qc_pt1$under70 <- "Yes"
dat_qc_pt2$totcases <- dat_qc_pt2$totalover70 
dat_qc_pt2$under70 <- "No"



dat_qc_pt1 <-  dplyr::select(dat_qc_pt1, c("Reported_Date", "totcases", "under70"))
dat_qc_pt2 <-  dplyr::select(dat_qc_pt2, c("Reported_Date", "totcases", "under70"))

mydat_QC <- rbind(dat_qc_pt1,dat_qc_pt2)

test_QCspline <- make_case_splines(mydat_BC)

mytest_QC = get_testprop(changedate = ymd("2021-12-21"), 
                         mysplines = test_QCspline, 
                         halftime = 35, steepness = 0.15)
glimpse(mytest_QC)
ggplot(mytest_QC, aes(x=date, y=test_prop))+geom_line()

# 
# lowerages = c("0 to 9 ans" ,"20 to 29", "30 to 39", "40 to 49", "50 to 59", "60 to 69", "unknown")
# agedat_qc$under70 = "No" 
# agedat_qc$under70[which(agedat_qc$agegroup10 %in% lowerages)] = "Yes"
# 
# mydat_QC =  group_by(agedat_qc, reporteddate , under70) %>% 
#   dplyr::summarise(totcases = sum(cases)) 
# mydat_QC$`Reported_Date` <- mydat_QC$reporteddate
# 
# test_QCspline  <- make_case_splines(mydat_QC)


