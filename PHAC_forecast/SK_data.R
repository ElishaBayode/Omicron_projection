
# Read and prepare data -----------------------------------------------------

#dat <- readr::read_csv("https://dashboard.saskatchewan.ca/export/cases/2097.csv")
location <- "CAN.csv"
linkRaw <- "https://raw.githubusercontent.com/ishaberry/Covid19Canada/master/timeseries_prov/cases_timeseries_prov.csv"
data <- readr::read_csv(linkRaw)
readr::write_csv(data, file.path("data", location))

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

ii_jul1 <- which(dat$date == ymd("2022-02-07"))
tot2d <- dat$adjust_cases[ii_jul1 + 1]
dat$adjust_cases[c(ii_jul1, ii_jul1 + 1)] <- realloc(tot2d, 2)

ggplot(dat, aes(date, adjust_cases)) +
  geom_point() +
  # x <- seq(0, 40, length.out = 200)
  # plot(x, dlnorm(x, log(12), 0.1), type = "l", xaxs = "i", yaxs = "i")
  scale_x_date(date_breaks = "2 weeks", date_labels = "%Y-%m-%d") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid.major = element_line(color = "grey"))

#
#dat$adjust_cases[nrow(dat)-1] <- round(dat$adjust_cases[nrow(dat)-2]/1.2,0)
#dat$adjust_cases[nrow(dat)-2] <- dat$adjust_cases[nrow(dat)-2]/2
tail(dat,20)
#dat$value <- dat$cases
dat$value <- dat$adjust_cases
dat <-dat[-nrow(dat)]  #dat[-c(nrow(dat)-1, nrow(dat)), ]

dat <- filter(dat, date >= "2020-09-22")
#dat <- filter(dat, date <= "2022-02-07")
##dat <- filter(dat, date <= "2021-12-23")
dat$day <- 1:length(dat$cases)
saveRDS(dat, file.path("data/SK-dat.rds"))
dat =readRDS(file.path("data/SK-dat.rds"))



#########SK test_prop #####################

#SK age structured data (????)

dat = get_saskatchewan_case_data()

# going to start this whole thin in fall 2021 
dat_sk <- filter(dat, Region == "Total")
dat_sk <- filter(dat_sk,  Date  >= ymd("2021-11-30")) 



totat_under70 <- diff(dat_sk$`Age 4 and Under`) + diff(dat_sk$`Age 5 to 11`) + diff(dat_sk$`Age 12 to 19`) +
                 diff(dat_sk$`Age 20 to 39`) + diff(dat_sk$`Ages 40 to 59`) + diff(dat_sk$`Age 60 to 79`)


