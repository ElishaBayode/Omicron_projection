
get_hosp_data <- function(intro_date, stop_date){
  hospdat <- readr::read_csv("data/all_hosps_byage1.csv") %>%
    mutate(week_of = case_when(week >= 44 ~ ymd("2021-10-31")+7*(week-44),
        week < 44 ~ ymd("2022-01-02") + 7*(week-1))) %>%
    filter(week_of >= intro_date & week_of <= stop_date) %>%
    group_by(week_of) %>%
    summarize(new = sum(new)) # combine ages for now
  return(hospdat)
}

get_IHR <- function(){
  hospdat <- get_hosp_data(ymd("2021-11-30"),ymd("2022-03-10")) # should be a lag but prob doesn't matter
  tot_hosp <- sum(hospdat$new)
  ## Seroprev for BA1 
  pop <- 5070000  
  prev <- 0.34 # from slides, approx. % infected between Dec 2021-Mar 2022
  IHR <- tot_hosp/(pop*prev)
return(IHR)
}

