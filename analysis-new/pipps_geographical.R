
#################################################
###   Adding geographical distn. of BC waves  ###
#################################################


## Use proj_out in this function:
# project_HAs(total_out = proj_out, which_wave_match = 6)

# 
project_HAs <- function(total_out, which_wave_match = 6){
  # total_out = time series of cases you want to divide by HA
  # which_wave_match = which past covid wave in BC do you want to match the projections to (1:7)
  
  ##---- Load historical regional data, to get HA proportions over time
  gdat <- read_csv("data/BCCDC_COVID19_Regional_Summary_Data.csv")
  
  ##---- Wave start/end cutoff dates
  wave0 <- min(gdat$Date)
  wave1 <- "2020-06-01"
  wave2 <- "2020-09-24"
  wave3 <- "2021-02-28"
  wave4 <- "2021-07-01"
  wave5 <- "2021-12-01"
  wave6 <- "2022-03-01"
  wave7 <- max(gdat$Date)
  wave_cutoffs <- as.Date(c(wave0, wave1,wave2,wave3,   wave4,   wave5,wave6, wave7))
  #                                wuhan wuhan wuhan alpha/gamma delta  ba1   ba2

  
  ##---- What proportion of cases where in each HA over time?
  gdat <- gdat  %>% filter(HA!="All") %>% filter(HA!="Out of Canada") %>%
    group_by(Date, HA) %>%
    summarise(n = sum( Cases_Reported_Smoothed)) %>%
    mutate(percentage = n / sum(n))
  # Filter geo data to desired wave only:
  gdat <- gdat %>% filter(Date>wave_cutoffs[which_wave_match] & Date<=wave_cutoffs[which_wave_match+1])
  
  
  ##---- Adjust length of proportions to match length of new wave, by fitting a spline for each HA and interpolating
  split.dfs <- split(gdat, gdat$HA)
  spline.HAs <- lapply(split.dfs, function (x) spline(x$Date, x$percentage, n = nrow(total_out)))
    
  
  ##---- Multiply HA proportions by total cases
  total_out$Fraser <- total_out$Total*spline.HAs$Fraser$y
  total_out$Coastal <- total_out$Total*spline.HAs$`Vancouver Coastal`$y
  total_out$Island <- total_out$Total*spline.HAs$`Vancouver Island`$y
  total_out$Interior <- total_out$Total*spline.HAs$Interior$y
  total_out$Northern <- total_out$Total*spline.HAs$Northern$y
  
  
  ##---- Plot
  HAs <- c("Coastal", "Fraser", "Interior", "Northern", "Island")
  out_plot <- pivot_longer(total_out, c(Total,HAs), names_to = "HA", values_to = "count") %>%
    ggplot(aes(x=date, y=count, colour=HA)) + geom_line() + ylab("Incident cases") + xlab("Date") + theme_minimal()

  
  return(out_plot)
}



# Checking 

#testing <- total_out[27:31]
#testing <- testing %>% rowwise() %>% mutate(m = sum(c(Fraser, Coastal, Island, Interior, Northern)))
#testing$m - total_out$Total
# Ok, the sum of the HAs almost exactly matches the overall total, I think this is close enough

