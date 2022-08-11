
#############################################################
###   Dividing up waves geographically (by HA) or by age  ###
#############################################################


## Use objects like proj_out in these functions:
# project_HAs(total_out = proj_out, which_wave_match = 6)

# 
project_HAs <- function(total_out, which_wave_match = 6, facets = FALSE){
  # total_out = time series of cases you want to divide by HA
  # which_wave_match = which past covid wave in BC do you want to match the projections to (1:7), 7 = ba.2
  
  ##---- Load historical regional data, to get HA proportions over time
  gdat <- read_csv("data/BCCDC_COVID19_Regional_Summary_Data.csv")
   gdat %>% filter(HA!="All") %>% filter(HSDA=="All") %>% ggplot(aes(x=Date, y=Cases_Reported), group = HA) +
   geom_line(aes(color=HA)) #+ ylim(0, 750) + xlim(as.Date("2021-01-01"), as.Date("2021-05-29"))
  
  ##---- Wave start/end cutoff dates
  wave0 <- "2020-03-01"
  wave1 <- "2020-06-01"
  wave2 <- "2020-09-30"
  wave3 <- "2021-02-16"
  wave4 <- "2021-07-01"
  wave5 <- "2021-12-01"
  wave6 <- "2022-03-01"
  wave7 <- "2022-06-02"
  wave_cutoffs <- as.Date(c(wave0, wave1,wave2,wave3,   wave4,   wave5,wave6, wave7))
  #                                wuhan wuhan wuhan alpha/gamma delta  ba1   ba2
  # Filter geo data to desired wave only:
  gdat <- gdat %>% filter(Date>wave_cutoffs[which_wave_match] & Date<=wave_cutoffs[which_wave_match+1])
  
  ##---- What proportion of cases where in each HA over time?
  # 'complete' = fill in missing HA groups on each day with 0 cases
  gdat %<>% filter(HA!="All") %>% filter(HA!="Out of Canada") %>% filter(HSDA=="All") %>% complete(Date, HA, fill = list(Cases_Reported = 0, Cases_Reported_Smoothed = 0)) %>%
    group_by(Date, HA) %>%
    summarise(n = sum( Cases_Reported)) %>%
    mutate(percentage = n / sum(n))
  gdat %<>% group_by(Date) %>% mutate(allHA_total = sum(n))

  
  
  ##---- Adjust length of proportions to match length of new wave, by fitting a spline for each HA and interpolating
  split.dfs <- split(gdat, gdat$HA)
  spline.HAs <- lapply(split.dfs, function (x) spline(x$Date, x$percentage, n = nrow(total_out)))
    
  
  ##---- Multiply HA proportions by total cases
  total_out$Total.Fraser <- total_out$Total*spline.HAs$Fraser$y
  total_out$Total.Coastal <- total_out$Total*spline.HAs$`Vancouver Coastal`$y
  total_out$Total.Island <- total_out$Total*spline.HAs$`Vancouver Island`$y
  total_out$Total.Interior <- total_out$Total*spline.HAs$Interior$y
  total_out$Total.Northern <- total_out$Total*spline.HAs$Northern$y
  
  
  ##---- Plot
  HAs <- c("Total.Coastal", "Total.Fraser", "Total.Interior", "Total.Northern", "Total.Island")
  out_plot <- pivot_longer(total_out, c(Total,HAs), names_to = "HA", values_to = "count") %>%
    ggplot(aes(x=date, y=count, colour=HA)) + geom_line() + ylab("Incident cases") + xlab("Date") + theme_minimal() +  
    scale_colour_discrete(labels = c("Total", "Coastal", "Fraser", "Interior", "Island", "Northern"))
  if (facets) {
    out_plot <- pivot_longer(total_out, c(Total,HAs), names_to = "HA", values_to = "count") %>%
      ggplot(aes(x=date, y=count, colour=HA)) + geom_line() + facet_wrap(~HA, scales = "free") + ylab("Incident cases") + xlab("Date") + theme_minimal() +  
      scale_colour_discrete(labels = c("Total", "Coastal", "Fraser", "Interior", "Island", "Northern"))
  }
  return(list(df = total_out, plot = out_plot))
}



# Checking 

#testing <- total_out[27:31]
#testing <- testing %>% rowwise() %>% mutate(m = sum(c(Fraser, Coastal, Island, Interior, Northern)))
#testing$m - total_out$Total
# Ok, the sum of the HAs almost exactly matches the overall total, I think this is close enough




# 
project_ages <- function(total_out, which_wave_match = 6, facets = FALSE){
  # total_out = time series of cases you want to divide by age group
  # which_wave_match = which past covid wave in BC do you want to match the projections to (1:7), 7 = ba.2
  
  ##---- Load historical age based data, to get age proportions over time
  adat <- read_csv("data/BCCDC_COVID19_Dashboard_Case_Details.csv") %>%
      count(Reported_Date, Age_Group) %>% filter(Age_Group!="Unknown")
  ggplot(adat, aes(x=Reported_Date, y=n), group = Age_Group) +
      geom_line(aes(color=Age_Group)) #+ ylim(0, 250) + xlim(as.Date("2021-07-01"), as.Date("2021-12-01"))
  
  
  ##---- Wave start/end cutoff dates
  wave0 <- "2020-03-01"
  wave1 <- "2020-06-01"
  wave2 <- "2020-09-30"
  wave3 <- "2021-02-16"
  wave4 <- "2021-07-01"
  wave5 <- "2021-12-01"
  wave6 <- "2022-03-01"
  wave7 <- "2022-06-02"
  wave_cutoffs <- as.Date(c(wave0, wave1,wave2,wave3,   wave4,   wave5,wave6, wave7))
  #                                wuhan wuhan wuhan alpha/gamma delta  ba1   ba2
  # Filter age data to desired wave only:
  adat <- adat %>% filter(Reported_Date>wave_cutoffs[which_wave_match] & Reported_Date<=wave_cutoffs[which_wave_match+1])
  

  ##---- What proportion of cases where in each age group over time?
  # First, fill in missing age groups on each day with 0 cases
  adat <- complete(adat, Reported_Date, Age_Group, fill = list(n = 0))
  adat %<>% 
    group_by(Reported_Date)  %>%
    mutate(percentage = n / sum(n))
 adat %<>% group_by(Reported_Date) %>% mutate(allage_total = sum(n))
  

  
  ##---- Adjust length of proportions to match length of new wave, by fitting a spline for each age group and interpolating
  split.dfs <- split(adat, adat$Age_Group)
  spline.ages <- lapply(split.dfs, function (x) spline(x$ Reported_Date, x$percentage, n = nrow(total_out)))
  
  
  ##---- Multiply age proportions by total cases
  total_out$"Total.<10" <- total_out$Total*spline.ages$"<10"$y
  total_out$"Total.10-19" <- total_out$Total*spline.ages$"10-19"$y
  total_out$"Total.20-29" <- total_out$Total*spline.ages$"20-29"$y
  total_out$"Total.30-39" <- total_out$Total*spline.ages$"30-39"$y
  total_out$"Total.40-49" <- total_out$Total*spline.ages$"40-49"$y
  total_out$"Total.50-59" <- total_out$Total*spline.ages$"50-59"$y
  total_out$"Total.60-69" <- total_out$Total*spline.ages$"60-69"$y
  total_out$"Total.70-79" <- total_out$Total*spline.ages$"70-79"$y
  total_out$"Total.80-89" <- total_out$Total*spline.ages$"80-89"$y
  total_out$"Total.90+" <- total_out$Total*spline.ages$"90+"$y
  
  
  ##---- Plot
  Ages <- paste0("Total.", names(spline.ages))
  out_plot <- pivot_longer(total_out, c(Total,Ages), names_to = "Ages", values_to = "count") %>%
    ggplot(aes(x=date, y=count, colour=Ages)) + geom_line() + ylab("Incident cases") + xlab("Date") + theme_minimal() + labs(colour = "Age Group") +
    scale_colour_discrete(labels = c("Total", names(spline.ages)))
  if (facets) {
    out_plot <- pivot_longer(total_out, c(Total,Ages), names_to = "Ages", values_to = "count") %>%
      ggplot(aes(x=date, y=count, colour=Ages))+  geom_line() + facet_wrap(~Ages, scales = "free") + ylab("Incident cases") + xlab("Date") + theme_minimal() + labs(colour = "Age Group") +
      scale_colour_discrete(labels = c("Total", names(spline.ages)))
  }
  
  return(list(df = total_out, plot = out_plot))
}




