

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



############ test_prop data #############