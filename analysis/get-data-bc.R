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
    dplyr::summarise(cases = n()) %>%
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
    dplyr::summarise(cases = n()) %>%
    filter(Reported_Date >= ymd("2021-09-01"))
# still need to patch gaps in the first few days 


# steps: 

# get a time series for <70 and 70+ 
lowerages = c("<10" , "10-19","20-29", "30-39",  "40-49" , "50-59" , "60-69")
agedat$under70 = "No" 
agedat$under70[which(agedat$`Age group` %in% lowerages)] = "Yes"

mydat =  group_by(agedat, Reported_Date, under70) %>% 
    dplyr::summarise(totcases = sum(cases)) 



ggplot(data = filter(mydat, Reported_Date > "2021-06-01"),  aes(x=Reported_Date, y=totcases, color=under70))+
    geom_point() + ylim(c(0,3000))

## get a smooth.spline for the log of 70+ and the log of under 70 
# (1) over 70
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


# ---- make a few different offsets, plot them together, and illustrate how 
# much the total cases depend on the offset assumptions


off1= getoffset(halftime=100,steepness = 0.05) # really flat, like Sally's
off2 = getoffset(halftime=15,steepness = 0.25) # changes over 60 days and not very steep
plot(off2)
L1 = ymd("2021-12-21")-min(pred$Reported_Date) # 
L2 = nrow(pred)- L1 # how many values to fill in with the new offset? 
offconst = upred$predlcases-pred$predlcases
offrevert = offconst
offconst[(L1+1):nrow(seriesoffset)] = off1[1:L2]
offrevert[(L1+1):nrow(seriesoffset)] = off2[1:L2]

# need a data frame with date, the cases by age band (mydat), and the model by age band 
# (pred and upred) 


df = mydat %>% ungroup() %>%  mutate(cases = totcases, 
                           date = Reported_Date, 
                           agegroup = ifelse(under70=="Yes", "Under 70", "Over 70")) %>% 
    select(date,cases,agegroup) 
bothmodels = rbind(pred %>% mutate(date = Reported_Date, model = exp(predlcases),
                                   adjmodel1=exp(predlcases),
                                   adjmodel2=exp(predlcases), 
                                   agegroup = "Over 70") %>% select(date, model, adjmodel1,
                                                                    adjmodel2, agegroup), 
                   upred  %>% mutate(date = Reported_Date, model = exp(predlcases),
                                     adjmodel2=exp(offrevert+pred$predlcases),
                                     adjmodel1 = exp(offconst + pred$predlcases),
                                     agegroup = "Under 70") %>% select(date, model,adjmodel1,
                                                                       adjmodel2, agegroup))
test  = merge(df, bothmodels)
ggplot(test, aes(x=date, y = cases, color=agegroup))+geom_point() +
    geom_line(aes(x=date, y=model)) + 
    geom_line(aes(x=date, y=adjmodel2),linetype= "dashed") + 
    geom_line(aes(x=date, y=adjmodel1),linetype= "dotdash")


# ---- just to show the offsets to the bc covid group 
library(tidyr)
tmp = data.frame(date= pred$Reported_Date, 
                 sally = offconst, 
                 rawdata = upred$predlcases-pred$predlcases,
                 anotherchoice = offrevert) %>% 
    pivot_longer(2:4, names_to = "type", values_to = "offset")
ggplot(filter(tmp, date > ymd("2021-12-01")), aes(x=date, y=offset, color=type))+geom_line()+ylim(c(0,4)) +
    ylab("Log of the scale factor between  <70 and 70+")






