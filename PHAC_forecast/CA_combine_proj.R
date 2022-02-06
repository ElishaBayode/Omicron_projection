library(erer)
library(gridExtra)#
library(data.table)
library(ggplot2)


intro_date <-  ymd("2021-11-30")
stop_date <- ymd("2022-01-30")
dat_ab <- readRDS("data/AB-dat.rds")
dat_ab <- dat_ab %>% filter(date >= intro_date & date <= stop_date) %>% select(c("value","day", "date"))
dat_bc  <- readRDS("data/BC-dat.rds")
dat_bc <- dat_bc %>% filter(date >= intro_date & date <= stop_date) %>% select(c("value","day", "date"))
dat_mb <- readRDS("data/MB-dat.rds")
dat_mb <- dat_mb %>% filter(date >= intro_date & date <= stop_date) %>% select(c("value","day", "date"))
dat_on <- readRDS("data/ON-dat.rds")
dat_on <- dat_on %>% filter(date >= intro_date & date <= stop_date) %>% select(c("value","day", "date"))
dat_qc <- readRDS("data/QC-dat.rds")
dat_qc <- dat_qc %>% filter(date >= intro_date & date <= stop_date) %>% select(c("value","day", "date"))
dat_sk <- readRDS("data/SK-dat.rds")
dat_sk <- dat_sk %>% filter(date >= intro_date & date <= stop_date) %>% select(c("value","day", "date"))



AB_Omicron = as.data.frame(readRDS(file.path("data/AB_test_constraints.rds"))) 

BC_Omicron  = as.data.frame(readRDS(file.path("data/BC_test_constraints.rds")))

MB_Omicron  = as.data.frame(readRDS(file.path("data/MB_test_constraints.rds")))

ON_Omicron  = as.data.frame(readRDS(file.path("data/ON_test_constraints.rds")))

QC_Omicron  = as.data.frame(readRDS(file.path("data/QC_test_constraints.rds")))

SK_Omicron = as.data.frame(readRDS(file.path("data/SK_test_constraints.rds")))


#SK_Omicron <- filter(SK_Omicron, date <= ymd("2022-01-09"))
#BC_Omicron$day <- 1:length(BC_Omicron$value)
#plot(ON_Omicron$inc_tot)

SK_Omicron$date

AB_Omicron_rel = as.data.frame(readRDS(file.path("data/AB_no_constraints.rds")))
BC_Omicron_rel  = as.data.frame(readRDS(file.path("data/BC_no_constraints.rds")))
MB_Omicron_rel = as.data.frame(readRDS(file.path("data/MB_no_constraints.rds")))
ON_Omicron_rel  = as.data.frame(readRDS(file.path("data/ON_no_constraints.rds")))
QC_Omicron_rel  = as.data.frame(readRDS(file.path("data/QC_no_constraints.rds")))
SK_Omicron_rel = as.data.frame(readRDS(file.path("data/SK_no_constraints.rds")))



AB_Omicron$day <- 1:92
BC_Omicron$day <- 1:92
MB_Omicron$day <- 1:92
ON_Omicron$day <- 1:92
QC_Omicron$day <- 1:92
SK_Omicron$day <- 1:92

dat_ab$day <- 1:62
dat_bc$day <- 1:62
dat_mb$day <- 1:62
dat_on$day <- 1:62
dat_qc$day <- 1:62
dat_sk$day <- 1:62

AB_Omicron_rel$day <- 1:92
BC_Omicron_rel$day <- 1:92
MB_Omicron_rel$day <- 1:92
ON_Omicron_rel$day <- 1:92
QC_Omicron_rel$day <- 1:92
SK_Omicron_rel$day <- 1:92


#------ quick fix ###### add day col in each script 
dat_combine <- list(dat_ab,dat_bc,dat_mb,dat_on,dat_qc,dat_sk)





showOmicron <- list(AB=AB_Omicron,BC=BC_Omicron,MB=MB_Omicron,
                    ON=ON_Omicron,QC=QC_Omicron,SK=SK_Omicron
)  

showOmicron_rel <- list(AB=AB_Omicron_rel,BC=BC_Omicron_rel,MB=MB_Omicron_rel,
                        ON=ON_Omicron_rel,QC=QC_Omicron_rel,SK=SK_Omicron_rel
)  



write.list(z = showOmicron, file = "data/adjusted.csv")
write.list(z = showOmicron_rel, file = "data/unadjusted.csv")
library(data.table)


#for plotting Canada combine 

showOmicron <- rbindlist(showOmicron, fill = TRUE)[,lapply(.SD, sum, na.rm = TRUE) , day]

showOmicron_rel <- rbindlist(showOmicron_rel, fill = TRUE)[,lapply(.SD, sum, na.rm = TRUE) , day]

dat_combine  <- rbindlist(dat_combine , fill = TRUE)[,lapply(.SD, sum, na.rm = TRUE) , day]







showOmicron <- data.frame(showOmicron)
showOmicron_rel <- data.frame(showOmicron_rel)

showOmicron$date <- AB_Omicron$date
showOmicron_rel$date <- AB_Omicron$date

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


dat_combine$date <-  dat_ab$date
cols <- c("Current, with testing constraints" = "darkgreen", "Without testing constraints"="orange")

  gg_combine  <- ggplot() + geom_line(data=showOmicron,aes(x=date,y=X50.),color='blue',size=1.2,alpha=0.4) +
  geom_ribbon(data=showOmicron,aes(x=date,ymin=X2.5.,ymax=X97.5.),fill='blue',alpha=0.1) +
  geom_point(data=dat_combine,aes(x=date, y=value),color='grey28', alpha=0.5) +
  geom_ribbon(data=showOmicron_rel, aes(x=date, ymin=X2.5.,ymax= X97.5.),fill='orange',alpha=0.1) +
  geom_line(data=showOmicron_rel,aes(x=date,y=X50., color="Without testing constraints"),size=1.2,alpha=0.5) +
  #geom_line(aes(y=typical),color='blue') +
  labs(y="Reported cases",x="Date") + ylim(c(0,200000)) + 
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

