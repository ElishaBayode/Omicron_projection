
AB_Omicron = as.data.frame(readRDS(file.path("data/AB-incid.rds")))
BC_Omicron  = as.data.frame(readRDS(file.path("data/BC-incid.rds")))
MB_Omicron  = as.data.frame(readRDS(file.path("data/MB-incid.rds")))
ON_Omicron  = as.data.frame(readRDS(file.path("data/ON-incid.rds")))
QC_Omicron  = as.data.frame(readRDS(file.path("data/QC-incid.rds")))
SK_Omicron  = as.data.frame(readRDS(file.path("data/SK-incid.rds")))



showOmicron <- list(AB=AB_Omicron,BC=BC_Omicron,MB=MB_Omicron,
                    ON=ON_Omicron,QC=QC_Omicron,SK=SK_Omicron)  
library(data.table)

showOmicron <- rbindlist(showOmicron, fill = TRUE)[,lapply(.SD, sum, na.rm = TRUE) , time]
showOmicron$date <- BC_Omicron$date
showOmicron <- data.frame(showOmicron)
cols <- c("Omicron" = "#009E73", "Total"="#D55E00")


#to continue ###### 














# g <- ggplot(filter(myout,date < ymd("2022-02-28")), aes(x = date)) +
#   geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2, fill = "darkcyan") +
#   geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.2, alpha = 0.2, fill = "darkcyan") +
#   geom_line(aes(x = date, y = mymean)) +
#   
#   geom_point(data = filter(candata, date<= ymd("2021-12-12")), aes(x = date, y = cases, alpha = 0.4),size=1.4) +
#   ylab("Reported cases") + 
#   scale_x_date(date_breaks = "months", date_labels = "%b") +
#   theme_light() +
#   geom_line(data=showOmicron, aes(x=date, y=inc_mut, col ="Omicron"), size=1.4) +
#   geom_line(data=showOmicron, aes(x=date, y=inc_tot, col ="Total"),  size=1.4) +
#   geom_point(data=showOmicron, aes(x=date, y=rcases), color="grey15", size=1.4, alpha=0.3) +
#   theme(
#     axis.text.x = element_text(size = 12),
#     axis.title.x = element_blank(),
#     legend.position = "bottom",
#     legend.text = element_text(size=15),
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