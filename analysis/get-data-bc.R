library(readr)
library(CanCovidData)
library(dplyr)
library(lubridate)
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

# ---- HOSPITAL STUFF ---- 

# Hospital admissions for Health Authorities scraped from dashboard. 
# Scrapes the dashboard 6pm every day and is updated shortly after.
hosp1path = "https://mountainmath.ca/bc_total_hospitalizations.csv"

# Sally posted 
# https://docs.google.com/spreadsheets/d/1KvX2bNs4hUYGY8Kk47c4SmFVHrKT_0vEYY0NqKQShZs/edit#gid=0


# Rob Dumont has : 
# https://docs.google.com/spreadsheets/d/1VxUA8Les5TK-Uarz-2Cw20B-r0xIXHYDX0eCU7u47NY/edit#gid=274465377
# from this, I downloaded the hosp and the age outcomes and put them in the BC folder 
rdage = read_csv(file = "~/Dropbox/COVID/BC/BC COVID tracking - SR - Age Outcomes.csv", skip = 1) 
rdhos = read_csv(file = "~/Dropbox/COVID/BC/BC COVID tracking - SR - Hospitalizations.csv")

rdcases = rdage[,2:14]
colnames(rdcases)[3]  = "0-9"

rdhosp = rdage[,c(1,2,15:25)]
colnames(rdhosp) = colnames(rdcases)



# ggplot(dat, aes(Reported_Date, cases)) +
#    geom_point() +
#    scale_x_date(date_breaks = "2 weeks", date_labels = "%Y-%m-%d") +
#    theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid.major = element_line(color = "grey"))

# was PUB-dat but moving to this
#
# dat = readRDS( file.path("data-generated/BC-dat.rds"))
# dat$day <- 1:length(dat$cases)

# saveRDS(dat, file.path("data-generated/BC-dat.rds"))

# ---- WASTEWATER ---- 

# from Jens in the data_sources channel in a thread 
# Actually, script to get the data is at the bottom of this helper script file, function called “get_data_for_plant”. https://github.com/mountainMath/BCCovidSnippets/blob/main/R/helpers.R
library(CanCovidData) # jens 

# Usage:
    plants <- c("Annacis Island","Iona Island","Lions Gate","Lulu Island","Northwest Langley")
wwd = wastewater_data <- plants %>% 
    lapply(get_data_for_plant) %>%
    bind_rows()
wwd$Plant = as.factor(wwd$Plant)
ggplot(wwd, aes(x=Date, y=Value,color=Plant)) + geom_point() + geom_line() + 
    facet_grid(~Plant)


######### BEGIN SECTION FOR ADJUSTING REPORTED CASES 
# ---- adjust for under-reporting based on continued testing in the 70+s ----


dat = get_british_columbia_case_data()

# going to start this whole thin in fall 2021 
agedat <- group_by(dat, Reported_Date, `Age group`) %>%
    summarise(cases = n()) %>%
    filter(Reported_Date >= ymd("2021-09-01"))
# still need to patch gaps in the first few days 


# steps: 

# get a time series for <70 and 70+ 
lowerages = c("<10" , "10-19","20-29", "30-39",  "40-49" , "50-59" , "60-69")
agedat$under70 = "No" 
agedat$under70[which(agedat$`Age group` %in% lowerages)] = "Yes"

mydat =  group_by(agedat, Reported_Date, under70) %>% 
    summarise(totcases = sum(cases)) 



ggplot(data = mydat,  aes(x=Reported_Date, y=totcases, color=under70))+geom_point()

# get a smooth.spline for the log of 70+ and the log of under 70 

over70 = filter(mydat, under70=="No")
logover70 = over70 %>% mutate(lcases = log(totcases)) %>% select(Reported_Date, lcases)
ospline = smooth.spline(logover70$Reported_Date, logover70$lcases, df=15)
                        

pred = data.frame(x = ospline$x, predlcases = ospline$y, 
                  Reported_Date = logover70$Reported_Date)

ggplot(data = logover70, aes(x=Reported_Date, y=lcases))+geom_point(color="blue", alpha=0.5) +
    geom_line(data =pred, aes(x=Reported_Date, y=predlcases))

under70 = filter(mydat, under70=="Yes")
logunder70 = under70 %>% mutate(lcases = log(totcases)) %>% select(Reported_Date, lcases)
uspline = smooth.spline(logunder70$Reported_Date, logunder70$lcases, df=15)


upred = data.frame(x = uspline$x, predlcases = uspline$y, 
                  Reported_Date = logunder70$Reported_Date)
ggplot(data = logunder70, aes(x=Reported_Date, y=lcases))+geom_point(color="blue", alpha=0.5) +
    geom_line(data =upred, aes(x=Reported_Date, y=predlcases))


ggplot(data = logunder70, aes(x=Reported_Date, y=lcases))+geom_point(color="blue", alpha=0.5) +
    geom_line(data =upred, aes(x=Reported_Date, y=predlcases)) +
    geom_point(data = logover70, aes(x=Reported_Date, y=lcases)) + 
    geom_line(data =pred, aes(x=Reported_Date, y=predlcases))



# get the offset between 70+ and under 70 cases at a particular date, 
# which is the difference between the two above splines at a fixed date 
mydate = ymd("2021-12-21") 
# offset at a particular date (Sally) 
offset = filter(upred,Reported_Date == mydate)$predlcases -
    filter(pred, Reported_Date == mydate)$predlcases 
# offset time series - this varies quite a bit, actually 
seriesoffset = filter(upred,Reported_Date <= mydate)$predlcases - 
                  filter(pred,Reported_Date <= mydate)$predlcases
plot(seriesoffset)
# mean offset over the time in the fall before my chosen date 
offset = mean(filter(upred,Reported_Date <= mydate)$predlcases - 
                  filter(pred,Reported_Date <= mydate)$predlcases)
# the problem is that this produces a model fit that is *under* the real reported cases in 
# those under 70 (because in Dec the offset was higher than the mean for the fall) 

# add this offset to the 70+ cases to get a model for the total cases. 

pred$totmodel = pred$predlcases + offset 

ggplot(data = logunder70, aes(x=Reported_Date, y=lcases))+geom_point(color="blue", alpha=0.5) +
    geom_line(data =upred, aes(x=Reported_Date, y=predlcases)) +
    geom_point(data = logover70, aes(x=Reported_Date, y=lcases)) + 
    geom_line(data =pred, aes(x=Reported_Date, y=predlcases)) + 
    geom_line(data = pred, aes(x=Reported_Date, y = totmodel), color="grey")


# now plot this with sensible y labels 
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
