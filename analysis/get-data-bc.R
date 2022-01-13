library(readr)
get_british_columbia_case_data <- function(){
    path="http://www.bccdc.ca/Health-Info-Site/Documents/BCCDC_COVID19_Dashboard_Case_Details.csv"
    read_csv(path,col_types=cols(.default="c")) %>%
       dplyr::rename(`Reported Date`=Reported_Date,`Health Authority`=HA,`Age group`=Age_Group) %>%
        mutate(`Age group`=recode(`Age group`,"19-Oct"="10-19")) %>%
        mutate(Reported_Date=as.Date(`Reported Date`,tryFormats = c("%Y-%m-%d", "%m/%d/%Y")))
}
dat = get_british_columbia_case_data()

dat <- group_by(dat, Reported_Date) %>%
    summarise(cases = n()) %>%
    filter(Reported_Date >= ymd("2020-02-27"))
# still need to patch gaps in the first few days 
dat$Reported_Date[1:2] <- c(ymd("2020-03-01"), ymd("2020-03-02"))


dat <- dat[order(dat$Reported_Date), ]

dat$day <- seq(1, nrow(dat))
dat$value <- dat$cases
dat$date <- dat$Reported_Date

tail(dat)
# check if this is required 
# dat <- dat[-nrow(dat), ]


# Hospital admissions for Health Authorities scraped from dashboard. 
# Scrapes the dashboard 6pm every day and is updated shortly after.
hosp1path = "https://mountainmath.ca/bc_total_hospitalizations.csv"

# Sally posted 
# https://docs.google.com/spreadsheets/d/1KvX2bNs4hUYGY8Kk47c4SmFVHrKT_0vEYY0NqKQShZs/edit#gid=0

# Rob Dumont has : 
# https://docs.google.com/spreadsheets/d/1VxUA8Les5TK-Uarz-2Cw20B-r0xIXHYDX0eCU7u47NY/edit#gid=274465377

# ggplot(dat, aes(Reported_Date, cases)) +
#    geom_point() +
#    scale_x_date(date_breaks = "2 weeks", date_labels = "%Y-%m-%d") +
#    theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid.major = element_line(color = "grey"))

# was PUB-dat but moving to this
#
# dat = readRDS( file.path("data-generated/BC-dat.rds"))
# dat$day <- 1:length(dat$cases)

# saveRDS(dat, file.path("data-generated/BC-dat.rds"))
