library(erer)
library(gridExtra)
AB_Omicron = as.data.frame(readRDS(file.path("data/AB-incid.rds")))
BC_Omicron  = as.data.frame(readRDS(file.path("data/BC-incid.rds")))
MB_Omicron  = as.data.frame(readRDS(file.path("data/MB-incid.rds")))
ON_Omicron  = as.data.frame(readRDS(file.path("data/ON-incid.rds")))
QC_Omicron  = as.data.frame(readRDS(file.path("data/QC-incid.rds")))
SK_Omicron  = as.data.frame(readRDS(file.path("data/SK-incid.rds")))


plot(ON_Omicron$inc_tot)


AB_Omicron_er = as.data.frame(readRDS(file.path("data/AB-incid_er.rds")))
BC_Omicron_er  = as.data.frame(readRDS(file.path("data/BC-incid_er.rds")))
MB_Omicron_er  = as.data.frame(readRDS(file.path("data/MB-incid_er.rds")))
ON_Omicron_er  = as.data.frame(readRDS(file.path("data/ON-incid_er.rds")))
QC_Omicron_er  = as.data.frame(readRDS(file.path("data/QC-incid_er.rds")))
SK_Omicron_er  = as.data.frame(readRDS(file.path("data/SK-incid_er.rds")))


showOmicron <- list(AB=AB_Omicron,BC=BC_Omicron,MB=MB_Omicron,
                    ON=ON_Omicron,QC=QC_Omicron,SK=SK_Omicron)  

showOmicron_er <- list(AB=AB_Omicron_er,BC=BC_Omicron_er,MB=MB_Omicron_er,
                    ON=ON_Omicron_er,QC=QC_Omicron_er,SK=SK_Omicron_er)  

write.list(z = showOmicron, file = "data/current.csv")
library(data.table)

showOmicron <- rbindlist(showOmicron, fill = TRUE)[,lapply(.SD, sum, na.rm = TRUE) , time]
#showOmicron$date <- SK_Omicron$date

showOmicron_er <- rbindlist(showOmicron_er, fill = TRUE)[,lapply(.SD, sum, na.rm = TRUE) , time]
showOmicron$date <- ON_Omicron$date





AB_Omicron_int = as.data.frame(readRDS(file.path("data/AB-incid_int.rds")))
BC_Omicron_int  = as.data.frame(readRDS(file.path("data/BC-incid_int.rds")))
MB_Omicron_int  = as.data.frame(readRDS(file.path("data/MB-incid_int.rds")))
ON_Omicron_int  = as.data.frame(readRDS(file.path("data/ON-incid_int.rds")))
QC_Omicron_int  = as.data.frame(readRDS(file.path("data/QC-incid_int.rds")))
SK_Omicron_int  = as.data.frame(readRDS(file.path("data/SK-incid_int.rds")))


plot(showOmicron$inc_tot)

AB_Omicron_er_int = as.data.frame(readRDS(file.path("data/AB-incid_er_int.rds")))
BC_Omicron_er_int  = as.data.frame(readRDS(file.path("data/BC-incid_er_int.rds")))
MB_Omicron_er_int  = as.data.frame(readRDS(file.path("data/MB-incid_er_int.rds")))
ON_Omicron_er_int  = as.data.frame(readRDS(file.path("data/ON-incid_er_int.rds")))
QC_Omicron_er_int  = as.data.frame(readRDS(file.path("data/QC-incid_er_int.rds")))
SK_Omicron_er_int  = as.data.frame(readRDS(file.path("data/SK-incid_er_int.rds")))


showOmicron_int <- list(AB=AB_Omicron_int,BC=BC_Omicron_int,MB=MB_Omicron_int,
                    ON=ON_Omicron_int,QC=QC_Omicron_int,SK=SK_Omicron_int)  

write.list(z = showOmicron_int, file = "data/reducetransm.csv")

showOmicron_er_int <- list(AB=AB_Omicron_er_int,BC=BC_Omicron_er_int,MB=MB_Omicron_er_int,
                       ON=ON_Omicron_er_int,QC=QC_Omicron_er_int,SK=SK_Omicron_er_int)  


library(data.table)

showOmicron_int <- rbindlist(showOmicron_int, fill = TRUE)[,lapply(.SD, sum, na.rm = TRUE) , time]
showOmicron$date <- ON_Omicron$date
showOmicron_int$date <- ON_Omicron$date
showOmicron_er_int <- rbindlist(showOmicron_er_int, fill = TRUE)[,lapply(.SD, sum, na.rm = TRUE) , time]
#showOmicron$date <- ON_Omicron$date


showOmicron <- data.frame(showOmicron)
showOmicron_er <- data.frame(showOmicron_er)

showOmicron_int <- data.frame(showOmicron_int)
showOmicron_er_int <- data.frame(showOmicron_er_int)



#to continue ###### 


newABshowSimple = readRDS(file.path("figs/AB-fig.rds"))
newBCshowSimple = readRDS(file.path("figs/BC-fig.rds"))
newMBshowSimple = readRDS(file.path("figs/MB-fig.rds"))
newONshowSimple = readRDS(file.path("figs/ON-fig.rds"))
newQCshowSimple = readRDS(file.path("figs/QC-fig.rds"))
newSKshowSimple = readRDS(file.path("figs/SK-fig.rds"))


combine_Omicron_freq <- grid.arrange(newONshowSimple, newQCshowSimple, newMBshowSimple, nrow = 3, ncol = 1)
ggsave(file="figs/Omicron_freqPAGE1.pdf", combine_Omicron_freq, width = 9, height = 12)
ggsave(file="figs/Omicron_freqPAGE1.png", combine_Omicron_freq, width = 9, height = 12)

combine_Omicron_freq <- grid.arrange(newSKshowSimple, newABshowSimple, newBCshowSimple, nrow = 3, ncol = 1)
ggsave(file="figs/Omicron_freqPAGE2.pdf", combine_Omicron_freq, width = 9, height = 12)
ggsave(file="figs/Omicron_freqPAGE2.png", combine_Omicron_freq, width = 9, height = 12)









g <- ggplot(filter(myout,date < ymd("2022-02-28")), aes(x = date)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2, fill = "darkcyan") +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.2, alpha = 0.2, fill = "darkcyan") +
   geom_line(aes(x = date, y = mymean)) +
   
   geom_point(data = filter(candata, date<= ymd("2021-12-12")), aes(x = date, y = cases, alpha = 0.4),size=1.4) +
   ylab("Reported cases") + 
   scale_x_date(date_breaks = "months", date_labels = "%b") +
   theme_light() +
  geom_line(data=showOmicron, aes(x=date, y=inc_mut, col ="Omicron"), size=1.4) +
   geom_line(data=showOmicron, aes(x=date, y=inc_tot, col ="Total"),  size=1.4) +
   geom_point(data=showOmicron, aes(x=date, y=rcases), color="grey15", size=1.4, alpha=0.3) +
   theme(
     axis.text.x = element_text(size = 12),
    axis.title.x = element_blank(),
     legend.position = "bottom",
     legend.text = element_text(size=15),
#     axis.text.y = element_text(size = 12),
#     axis.title.y = element_text(size = 12),
#     #        panel.grid.major=element_line(color="grey"),
#     axis.line = element_line(colour = "grey", size = 1.5),
#     plot.margin = margin(1, 1, 1, 1, "cm")
#   ) +
#   scale_color_manual(values = cols) + 
#   guides(alpha = FALSE) + labs(colour ="") +
#   coord_cartesian(expand = FALSE, ylim = c(0, 50000))
# 
# # this version is good apparently - have checked.
# 
# g