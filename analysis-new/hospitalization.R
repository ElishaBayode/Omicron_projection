
# data description: each province is there with BC_pt etc 
#    'TotalCases':['pt'],
#    'DailyCases':['pd'],
#    'TotalDeaths':['dt'],
#    'DailyDeaths':['dd'],
#    'TotalHospitalized':['ht'],
#    'TotalICU':['it']
# Note Dean's github home for the pypm modelling is https://pypm.github.io/home/ 


deandata = read_csv("https://raw.githubusercontent.com/pypm/data/master/covid19/Canada/ca-pypm.csv")
glimpse(deandata)


# here, despite the complexities, using a dead simple 20 day lag and a factor of 
# 0.5 does a reasonble job, despite the huge testing limitations. 
# this may indicate that while both severity and testing have decreased a lot,
# we maintain testing in those at risk of hospitalization and it balances out? 


ggplot(deandata, aes(x=date, y= `BC-pd`))+geom_point(alpha=0.2)+
    geom_line(data=deandata, aes(x=date, y=`BC-ht`))+
    geom_point(data = deandata, aes(x=date, y = 0.5*lag(`BC-pd`,n=20)), color="blue", alpha=0.5)+
    ylim(c(0, 1000)) + xlim(c(ymd("2021-01-01"), ymd("2022-01-30")))



# anyway, I guess it provides a simple way to either (a) just project hosps forward in 
# teim from the model or (b) estimate something about the change in severity that is 
# implied by whatever we think about the changes in testing. 
# or, if we knew both the changes in testing and we were somewhat confident in changes 
# in severity from estimates in other places, we could fit the model using hosp data 
# as well as case data. 
