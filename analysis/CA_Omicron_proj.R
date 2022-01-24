library(erer)
library(gridExtra)#
library(data.table)
library(ggplot2)
AB_Omicron = as.data.frame(readRDS(file.path("data/AB_dat_sim.rds")))
BC_Omicron  = as.data.frame(readRDS(file.path("data/BC_dat_sim.rds")))
MB_Omicron  = as.data.frame(readRDS(file.path("data/MB_dat_sim.rds")))
ON_Omicron  = as.data.frame(readRDS(file.path("data/ON_dat_sim.rds")))
QC_Omicron  = as.data.frame(readRDS(file.path("data/QC_dat_sim.rds")))
SK_Omicron = as.data.frame(readRDS(file.path("data/SK_dat_sim.rds")))

#SK_Omicron <- filter(SK_Omicron, date <= ymd("2022-01-09"))
#BC_Omicron$day <- 1:length(BC_Omicron$value)
#plot(ON_Omicron$inc_tot)

QC_Omicron$date

AB_Omicron_rel = as.data.frame(readRDS(file.path("data/AB_dat_sim_rel.rds")))
BC_Omicron_rel  = as.data.frame(readRDS(file.path("data/BC_dat_sim_rel.rds")))
MB_Omicron_rel = as.data.frame(readRDS(file.path("data/MB_dat_sim_rel.rds")))
ON_Omicron_rel  = as.data.frame(readRDS(file.path("data/ON_dat_sim_rel.rds")))
QC_Omicron_rel  = as.data.frame(readRDS(file.path("data/QC_dat_sim_rel.rds")))
SK_Omicron_rel = as.data.frame(readRDS(file.path("data/SK_dat_sim_rel.rds")))






AB_Omicron_er = as.data.frame(readRDS(file.path("data/AB_forecast.rds")))
BC_Omicron_er  = as.data.frame(readRDS(file.path("data/BC_forecast.rds")))
MB_Omicron_er  = as.data.frame(readRDS(file.path("data/MB_forecast.rds")))
ON_Omicron_er  = as.data.frame(readRDS(file.path("data/ON_forecast.rds")))
QC_Omicron_er  = as.data.frame(readRDS(file.path("data/QC_forecast.rds")))
SK_Omicron_er  = as.data.frame(readRDS(file.path("data/SK_forecast.rds")))



AB_Omicron_er$day <- 1:30
BC_Omicron_er$day <- 1:30
MB_Omicron_er$day <- 1:30
ON_Omicron_er$day <- 1:30
QC_Omicron_er$day <- 1:30
SK_Omicron_er$day <- 1:30

AB_Omicron_int = as.data.frame(readRDS(file.path("data/AB_forecast_int.rds")))
BC_Omicron_int  = as.data.frame(readRDS(file.path("data/BC_forecast_int.rds")))
MB_Omicron_int  = as.data.frame(readRDS(file.path("data/MB_forecast_int.rds")))
ON_Omicron_int  = as.data.frame(readRDS(file.path("data/ON_forecast_int.rds")))
QC_Omicron_int  = as.data.frame(readRDS(file.path("data/QC_forecast_int.rds")))
SK_Omicron_int  = as.data.frame(readRDS(file.path("data/SK_forecast_int.rds")))
#------ quick fix ###### add day col in each script 
AB_Omicron_int$day <- 1:30
BC_Omicron_int$day <- 1:30
MB_Omicron_int$day <- 1:30
ON_Omicron_int$day <- 1:30
QC_Omicron_int$day <- 1:30
SK_Omicron_int$day <- 1:30


showOmicron <- list(AB=AB_Omicron,BC=BC_Omicron,MB=MB_Omicron,
                    ON=ON_Omicron,QC=QC_Omicron,SK=SK_Omicron
                    )  

showOmicron_rel <- list(AB=AB_Omicron_rel,BC=BC_Omicron_rel,MB=MB_Omicron_rel,
                    ON=ON_Omicron_rel,QC=QC_Omicron_rel,SK=SK_Omicron_rel
)  
showOmicron_er <- list(AB=AB_Omicron_er,BC=BC_Omicron_er,MB=MB_Omicron_er,
                    ON=ON_Omicron_er,QC=QC_Omicron_er,SK=SK_Omicron_er
                    )

showOmicron_int <- list(AB=AB_Omicron_int,BC=BC_Omicron_int,MB=MB_Omicron_int,
                       ON=ON_Omicron_int,QC=QC_Omicron_int,SK=SK_Omicron_int
                       )



write.list(z = showOmicron, file = "data/fits.csv")
write.list(z = showOmicron_er, file = "data/projection_current.csv")
write.list(z = showOmicron_int, file = "data/projection_adjusted.csv")
write.list(z = showOmicron_rel, file = "data/test_adjusted_fit.csv")
library(data.table)


#for plotting Canada combine 

showOmicron <- rbindlist(showOmicron, fill = TRUE)[,lapply(.SD, sum, na.rm = TRUE) , day]

showOmicron_rel <- rbindlist(showOmicron_rel, fill = TRUE)[,lapply(.SD, sum, na.rm = TRUE) , day]
showOmicron_rel$date <- QC_Omicron_rel$date

showOmicron_er <- rbindlist(showOmicron_er, fill = TRUE)[,lapply(.SD, sum, na.rm = TRUE) , day]
showOmicron_er$date <- QC_Omicron_er$date

showOmicron_int <- rbindlist(showOmicron_int, fill = TRUE)[,lapply(.SD, sum, na.rm = TRUE) , day]
showOmicron_int$date <- QC_Omicron_int$date






#plot(showOmicron$inc_tot)
#
#AB_Omicron_er_int = as.data.frame(readRDS(file.path("data/AB-incid_er_int.rds")))
#BC_Omicron_er_int  = as.data.frame(readRDS(file.path("data/BC-incid_er_int.rds")))
#MB_Omicron_er_int  = as.data.frame(readRDS(file.path("data/MB-incid_er_int.rds")))
#ON_Omicron_er_int  = as.data.frame(readRDS(file.path("data/ON-incid_er_int.rds")))
#QC_Omicron_er_int  = as.data.frame(readRDS(file.path("data/QC-incid_er_int.rds")))
#SK_Omicron_er_int  = as.data.frame(readRDS(file.path("data/SK-incid_er_int.rds")))


#showOmicron_int <- list(AB=AB_Omicron_int,BC=BC_Omicron_int,MB=MB_Omicron_int,
#                    ON=ON_Omicron_int,QC=QC_Omicron_int,SK=SK_Omicron_int)  

#write.list(z = showOmicron_int, file = "data/reducetransm.csv")

#showOmicron_er_int <- list(AB=AB_Omicron_er_int,BC=BC_Omicron_er_int,MB=MB_Omicron_er_int,
 #                      ON=ON_Omicron_er_int,QC=QC_Omicron_er_int,SK=SK_Omicron_er_int)  




#showOmicron_int <- rbindlist(showOmicron_int, fill = TRUE)[,lapply(.SD, sum, na.rm = TRUE) , time]
#showOmicron$date <- ON_Omicron$date
#showOmicron_int$date <- ON_Omicron$date
#showOmicron_er_int <- rbindlist(showOmicron_er_int, fill = TRUE)[,lapply(.SD, sum, na.rm = TRUE) , time]
#showOmicron$date <- ON_Omicron$date


showOmicron <- data.frame(showOmicron)
showOmicron_er <- data.frame(showOmicron_er)
showOmicron_int <- data.frame(showOmicron_int)



#to continue ###### 


newABshowSimple = readRDS(file.path("figs/AB-fig.rds"))
newBCshowSimple = readRDS(file.path("figs/BC-fig.rds"))
newMBshowSimple = readRDS(file.path("figs/MB-fig.rds"))
newONshowSimple = readRDS(file.path("figs/ON-fig.rds"))
newQCshowSimple = readRDS(file.path("figs/QC-fig.rds"))
newSKshowSimple = readRDS(file.path("figs/SK-fig.rds"))


combine_Omicron_freq <- grid.arrange(newONshowSimple, newQCshowSimple,  newMBshowSimple,nrow = 3, ncol = 1) #
ggsave(file="figs/Omicron_freqPAGE1.pdf", combine_Omicron_freq, width = 9, height = 12)
ggsave(file="figs/Omicron_freqPAGE1.png", combine_Omicron_freq, width = 9, height = 12)

combine_Omicron_freq <- grid.arrange( newABshowSimple, newBCshowSimple,newSKshowSimple, nrow = 3, ncol = 1) #
ggsave(file="figs/Omicron_freqPAGE2.pdf", combine_Omicron_freq, width = 9, height = 12)
ggsave(file="figs/Omicron_freqPAGE2.png", combine_Omicron_freq, width = 9, height = 12)

showOmicron$date <- BC_Omicron_rel$date

cols <- c("Current, with testing constraints" = "darkgreen", "Without testing constraints"="orange")

gg_combine  <- ggplot() + geom_line(data=showOmicron,aes(x=date,y=X50.),color='blue',size=1.2,alpha=0.4) +
   geom_ribbon(data=showOmicron,aes(x=date,ymin=X2.5.,ymax=X97.5.),fill='blue',alpha=0.1) +
   geom_point(data=showOmicron,aes(x=date, y=value),color='grey28', alpha=0.5) +
   geom_line(data=showOmicron_er,aes(x=date,y=X50.),color="darkgreen",size=1.2,alpha=0.5) +
   geom_ribbon(data=showOmicron_er,aes(x=date,ymin=X2.5.,ymax=X97.5.),fill='darkgreen',alpha=0.1) +
   geom_line(data=showOmicron_int,aes(x=date,y=X50., color="Without testing constraints"),size=1.2,alpha=0.5) +
   geom_ribbon(data=showOmicron_rel,aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill='orange',alpha=0.1) +
   geom_line(data=showOmicron_rel,aes(x=date,y=`50%`, color="Without testing constraints"),size=1.2,alpha=0.5) +
   geom_ribbon(data=showOmicron_int,aes(x=date,ymin=X2.5.,ymax=X97.5.),fill='orange',alpha=0.1)+
   #geom_line(aes(y=typical),color='blue') +
   labs(y="Reported cases",x="Date") + ylim(c(0,250000)) + 
   scale_x_date(date_breaks = "10 days", date_labels = "%b-%d-%y") +theme_light() +
   scale_color_manual(values = cols) +  theme(axis.text=element_text(size=15),
                                              plot.title = element_text(size=15, face="bold"),
                                              axis.text.x = element_text(angle = 45, hjust = 1),
                                              legend.position = "bottom", legend.title = element_text(size=15),
                                              legend.text = element_text(size=15),
                                              axis.title=element_text(size=15,face="bold")) +
   
   labs(color = " ",title="CA")

gg_combine

ggsave(file="figs/Canada.png", gg_combine, width = 10, height = 8)
ggsave(file="figs/Canada.pdf", gg_combine, width = 10, height = 8)
saveRDS(gg_combine, file.path("figs/Canada.rds"))
# # this version is good apparently - have checked.
# 
# g