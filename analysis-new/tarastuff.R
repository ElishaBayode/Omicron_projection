taram = readr::read_csv("~/Dropbox/COVID/BC/scenario-planning/Canadian COVID data_Infections_Line chart.csv")
glimpse(taram) 
taram = taram %>% mutate(chardate = date) %>% mutate(date = dmy(chardate))
ggplot(taram, aes(x=date, y = BC))+geom_point(alpha=0.5)+
  scale_x_date(date_breaks = "months", date_labels = "%b-%d") +theme_light()
# tara says broadly consistent wiht serology 
# but i get 
ba1tara = filter(taram, date > ymd("2021-12-01") & date < ymd("2022-03-30"))
sum(ba1tara$BC)/N
# so yes this is broadly consistent with the serology 