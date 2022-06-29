# ---- case data only for VCH ----

bcpub <- get_british_columbia_case_data() #readr::read_csv("data-raw/BCCDC_COVID19_Dashboard_Case_Details.csv")
bcpub$date <- lubridate::ymd(bcpub$`Reported Date`)

vch <-filter(bcpub, `Health Authority` == "Vancouver Coastal") %>%  
    group_by( date) %>%
  dplyr::summarise(cases = n()) 


# ---- read PHAC wastewater data and fix NAs that stop splines ----
bcww  = read_csv("~/Dropbox/COVID/PHAC-forecasts/wastewater/wastewaterdatanmlstatcan/MetroVancouver_20220131_NML_results.csv")
bcww= read_csv("~/Dropbox/COVID/PHAC-forecasts/wastewater/rewastewatersanitycheckquickquestion/MetroVancouver_20220207_NML_results.csv")

bcww$S_N1_1_cp = na.fill(as.numeric(bcww$S_N1_1_cp),"extend")
bcww$S_N1_2_cp = na.fill(as.numeric(bcww$S_N1_2_cp),"extend")
bcww$S_N2_1_cp = na.fill(as.numeric(bcww$S_N2_1_cp),"extend")
bcww$S_N2_2_cp = na.fill(as.numeric(bcww$S_N2_2_cp),"extend")

bcww$C_N1_1_cp = na.fill(as.numeric(bcww$C_N1_1_cp),"extend")
bcww$C_N1_2_cp = na.fill(as.numeric(bcww$C_N1_2_cp),"extend")
bcww$C_N2_1_cp = na.fill(as.numeric(bcww$C_N2_2_cp),"extend")
bcww$C_N2_2_cp = na.fill(as.numeric(bcww$C_N2_2_cp),"extend")

bcww$S_N1_1 = na.fill(as.numeric(bcww$S_N1_1_cp),"extend")
bcww$S_N1_2 = na.fill(as.numeric(bcww$S_N1_2_cp),"extend")
bcww$S_N2_1 = na.fill(as.numeric(bcww$S_N2_1_cp),"extend")
bcww$S_N2_2 = na.fill(as.numeric(bcww$S_N2_2_cp),"extend")

bcww$C_N1_1 = na.fill(as.numeric(bcww$C_N1_1_cp),"extend")
bcww$C_N1_2 = na.fill(as.numeric(bcww$C_N1_2_cp),"extend")
bcww$C_N2_1 = na.fill(as.numeric(bcww$C_N2_2_cp),"extend")
bcww$C_N2_2 = na.fill(as.numeric(bcww$C_N2_2_cp),"extend") # really should learn to do this with map


# ----explore what we have. Make some notes on the variables  ----
# eg do we really have flow or concentration or whatever
ggplot(filter(bcww,Date_sampled > ymd("2021-01-30") ), 
       aes(x=Date_sampled, C_N1_1_cp,color=Location)) +
  geom_line() + facet_wrap(~Location, scales="free") +
  geom_line(data = filter(vch,date >ymd("2021-01-30") & 
                            date < ymd("2022-02-01")) , 
            inherit.aes = F, aes(x=date, y = cases/5), color="black") # +

# add lots of the N2 and N1 data together , and again informally compare to cases
ggplot(filter(bcww,Date_sampled > ymd("2021-01-30") ), 
       aes(x=Date_sampled, C_N1_1_cp+C_N2_1_cp+S_N1_1_cp+S_N2_2_cp +
             C_N1_1+C_N2_1+S_N1_1+S_N2_2 ,
           color=Location)) +
  geom_line() + facet_wrap(~Location, scales="free") +
  geom_line(data = filter(vch,date >ymd("2021-01-30") & 
                            date < ymd("2022-02-01")) , 
            inherit.aes = F, aes(x=date, y = 2*cases), color="black")
# here you see the eventual conclusion from the statistical model: there is no substantial O wave 


# geom_line(aes(x=Date_sampled, y=S_PMMV_2_cp),color="grey")
# -- ** use _cp or the _avg for the actual virus copies. N1 and N2 have somewhat different signals as do the two replicates
# -- i doubt much will be lost by using the average of the two replicates
# I don't know whether to use the _cp variables for pmmv or the 'raw' ones 
# remove S_N1N2_PMMV and S_N1_PMMV etc before statistical model.
# these S_N1 etc are just 0 and 3. and infltemp, empty, and tsample
# the raw replicates are very similar. just use their average if using raw at all  
# same for the _cp. doesn't matter which one you use, mostly; could just average 
# the _PMMV_avg is the avg of the two _cp fields 
# PMMV: the _cp data have huge ranges. the non-cp data have reasonable ranges. But whatever they did 
# to get _cp for the virus part makes way more sense than the raw ones. so normalizing by the 
# raw pmmv seems wrong; normalizing by these huge numbers seems wrong too  esp in VLG where there 
# is a spike that is 10x as much as every other point
# rain_t  1 or 0 , as is  rain_y (but they arne't the same - mostly the same though) 
# temp: varies around 2-5. ph - changes, as exp around 7ish
# dinflvol: varies a lot between sites. varies by orders of magnitude. ASK ABOUT THIS 
# tss: varies over time. ASK ABOUT THIS . looks like a seasona temp but units are wrong? 
# conclusion: for statistical models, the pmmv varibles, the N1 and N2 obvs, 
# the temp , tss, rain_t and rain_y. Exclude the text fields and the above. 

# for spline work to normalize by PMMV I am a bit unclear because the PMMV signals are so strange
# QUESTIONS: -- dinflvol, tss, how do they get the _cp, why so variable PMMV in cp but not raw, 


# plot illustrating 4 data streams together, N1, N2, solid, concentrate 
ggplot(filter(bcww,Date_sampled > ymd("2021-01-30") & Date_sampled < ymd("2022-02-01")), 
       aes(x=Date_sampled, y=S_N1_1,color=Location))+
  geom_line() + facet_wrap(~Location) + 
  geom_line(aes(x=Date_sampled, y = S_N2_1)) + 
  geom_line(aes(x=Date_sampled, y = C_N2_1), linetype = "dotdash") +
  geom_line(aes(x=Date_sampled, y = C_N1_1), linetype = "dotdash")  +
  geom_line(data = filter(vch,date >ymd("2021-01-30") & 
                            date < ymd("2022-02-01")) , 
            inherit.aes = F, aes(x=date, y = cases/5), color="black")

ggplot(filter(bcww,Date_sampled > ymd("2021-01-30") & Date_sampled < ymd("2021-12-01")), 
       aes(x=Date_sampled, y=S_PMMV_avg,color=Location))+geom_line()+ facet_wrap(~Location)



# ideas: smoothed version, or average, of PMMV 1 and 2, and 1.10 (either S or C) 
# for normalization - to see if the normalized versions of S/C N1/2 replicates 1, 2 correlate 
# any better with cases in 2021 (pre-December) than the S/C things on their own 

# --- make some splines ---- 
# splines: i want splines for all of the virus count signals - given above, i want splines for the cp and the raw, 
# and i don't want to keep replicates separate, but i do want to keep N1 and N2 separate 
# i want to normalize each pmmv signal, and i guess I want a smoothed version too 
myspar=0.6
bcww = bcww %>% filter(Date_sampled > ymd("2021-02-05")) %>%  group_by(Location) %>%
  mutate(s_pmmv_avg_spl = smooth.spline(S_PMMV_avg, spar = myspar)$y,
         c_pmmv_avg_spl = smooth.spline(C_PMMV_avg, spar = myspar)$y,
         s_n1_spl = smooth.spline((S_N1_1+S_N1_2)/2, spar = myspar)$y,
         s_n2_spl = smooth.spline((S_N2_1+S_N2_2)/2, spar = myspar)$y,
         c_n1_spl = smooth.spline((C_N1_1+C_N1_2)/2, spar = myspar)$y,
         c_n2_spl = smooth.spline((C_N2_1+C_N2_2)/2, spar = myspar)$y,
         s_n1_spl_cp = smooth.spline((S_N1_1_cp+S_N1_2_cp)/2, spar = myspar)$y,
         s_n2_spl_cp = smooth.spline((S_N2_1_cp+S_N2_2_cp)/2, spar = myspar)$y,
         c_n1_spl_cp = smooth.spline((C_N1_1_cp+C_N1_2_cp)/2, spar = myspar)$y,
         c_n2_spl_cp = smooth.spline((C_N2_1_cp+C_N2_2_cp)/2, spar = myspar)$y) %>% 
  mutate(s_pmmv_avg_spl_norm = s_pmmv_avg_spl / mean(s_pmmv_avg_spl),
         c_pmmv_avg_spl_norm = c_pmmv_avg_spl / mean(c_pmmv_avg_spl)) %>% 
  mutate(s_n1_spl_norm = s_n1_spl / s_pmmv_avg_spl_norm, # SO suggest using the raw PMMV not exponentiated 
         s_n2_spl_norm = s_n2_spl / s_pmmv_avg_spl_norm, # and including the diluted version 
         c_n1_spl_norm = c_n1_spl / c_pmmv_avg_spl_norm,
         c_n2_spl_norm = c_n2_spl / c_pmmv_avg_spl_norm,
         s_n1_spl_cp_norm = s_n1_spl_cp / s_pmmv_avg_spl_norm, 
         s_n2_spl_cp_norm = s_n2_spl_cp / s_pmmv_avg_spl_norm,
         c_n1_spl_cp_norm = c_n1_spl_cp / c_pmmv_avg_spl_norm,
         c_n2_spl_cp_norm = c_n2_spl_cp / c_pmmv_avg_spl_norm)
# now i have solid/conc, n1/n2, cp/raw, and norm/non-norm. 
# this is too much. 

# but now i make some plots of these. let's put n1 and n2 together. i'll have 8 plots. 

# (1) solid, cp, non-norm . this is comparable to the S_N1_avg because it's not norm'd. 
mydf = bcww %>% mutate(N1 = s_n1_spl_cp, N2 = s_n2_spl_cp) %>% 
  select(Date_sampled, N1, N2, Location) %>% 
  pivot_longer(cols = c(2,3), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled) # splines
mydf2 = bcww %>% mutate(N1= S_N1_avg,N2= S_N2_avg) %>% 
  select(Date_sampled, N1, N2, Location) %>% 
  pivot_longer(cols = c(2,3), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled) # data 
ggplot(mydf, aes(x=date, y=count, color = target))+  geom_line()+facet_wrap(~Location) +
  geom_point(data=mydf2, aes(x=date, y=count, color=target)) +facet_wrap(~Location) #+ 
#  ylim(c(0,100))
# conclusion: the myspar parameter affects the apparent peak timing. Bumpy splines peak in early January. 
# smoother ones peak later but looks like still in January. 

# (2) concentrate, cp, non-norm 
mydf = bcww %>% mutate(N1 = c_n1_spl_cp, N2 = c_n2_spl_cp) %>% 
  select(Date_sampled, N1, N2, Location) %>% 
  pivot_longer(cols = c(2,3), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled) # splines
mydf2 = bcww %>% mutate(N1= C_N1_avg,N2= C_N2_avg) %>% 
  select(Date_sampled, N1, N2, Location) %>% 
  pivot_longer(cols = c(2,3), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled) # data 
ggplot(mydf, aes(x=date, y=count, color = target))+  geom_line()+facet_wrap(~Location) +
  geom_point(data=mydf2, aes(x=date, y=count, color=target)) +facet_wrap(~Location)
# conclusion is similar. early January. 

# (3) solid, cp, norm . this is comparable to the S_N1_avg/norm of s pmmv spline 
mydf = bcww %>% mutate(N1 = s_n1_spl_cp_norm, N2 = s_n2_spl_cp_norm) %>% 
  select(Date_sampled, N1, N2, Location) %>% 
  pivot_longer(cols = c(2,3), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled) # splines
mydf2 = bcww %>% mutate(N1= S_N1_avg/s_pmmv_avg_spl_norm,N2= S_N2_avg/s_pmmv_avg_spl_norm) %>% 
  select(Date_sampled, N1, N2, Location) %>% 
  pivot_longer(cols = c(2,3), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled) # data 
ggplot(mydf, aes(x=date, y=count, color = target))+  geom_line()+facet_wrap(~Location) +
  geom_point(data=mydf2, aes(x=date, y=count, color=target)) +facet_wrap(~Location) 
# i don't love this, because it doesn't have the delta wave in it. what does it matter if this peaks? 

# (4) concentrate, cp, norm . this is comparable to the C_N1_avg/norm of c pmmv spline 
mydf = bcww %>% mutate(N1 = c_n1_spl_cp_norm, N2 = c_n2_spl_cp_norm) %>% 
  select(Date_sampled, N1, N2, Location) %>% 
  pivot_longer(cols = c(2,3), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled) # splines
mydf2 = bcww %>% mutate(N1= C_N1_avg/c_pmmv_avg_spl_norm,N2= C_N2_avg/c_pmmv_avg_spl_norm) %>% 
  select(Date_sampled, N1, N2, Location) %>% 
  pivot_longer(cols = c(2,3), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled) # data 
ggplot(mydf, aes(x=date, y=count, color = target))+  geom_line()+facet_wrap(~Location) +
  geom_point(data=mydf2, aes(x=date, y=count, color=target)) +facet_wrap(~Location) 
# also misses the delta wave. peak time conclusions seem similar. 

# (5) as 1, but not with cp, with raw instead 
mydf = bcww %>% mutate(N1 = s_n1_spl, N2 = s_n2_spl) %>% 
  select(Date_sampled, N1, N2, Location) %>% 
  pivot_longer(cols = c(2,3), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled) # splines
mydf2 = bcww %>% mutate(N1= S_N1_1,N2= S_N2_1) %>% 
  select(Date_sampled, N1, N2, Location) %>% 
  pivot_longer(cols = c(2,3), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled) # data 
ggplot(mydf, aes(x=date, y=count, color = target))+  geom_line()+facet_wrap(~Location) +
  geom_point(data=mydf2, aes(x=date, y=count, color=target)) +facet_wrap(~Location)

# (6) concentrate, raw, non-norm 
mydf = bcww %>% mutate(N1 = c_n1_spl, N2 = c_n2_spl) %>% 
  select(Date_sampled, N1, N2, Location) %>% 
  pivot_longer(cols = c(2,3), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled) # splines
mydf2 = bcww %>% mutate(N1= C_N1_1,N2= C_N2_1) %>% 
  select(Date_sampled, N1, N2, Location) %>% 
  pivot_longer(cols = c(2,3), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled) # data 
ggplot(mydf, aes(x=date, y=count, color = target))+  geom_line()+facet_wrap(~Location) +
  geom_point(data=mydf2, aes(x=date, y=count, color=target)) +facet_wrap(~Location)


# (7) solid, raw, norm . this is comparable to the S_N1_avg/norm of s pmmv spline 
mydf = bcww %>% mutate(N1 = s_n1_spl_norm, N2 = s_n2_spl_norm) %>% 
  select(Date_sampled, N1, N2, Location) %>% 
  pivot_longer(cols = c(2,3), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled) # splines
mydf2 = bcww %>% mutate(N1= S_N1_1/s_pmmv_avg_spl_norm,N2= S_N2_1/s_pmmv_avg_spl_norm) %>% 
  select(Date_sampled, N1, N2, Location) %>% 
  pivot_longer(cols = c(2,3), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled) # data 
ggplot(mydf, aes(x=date, y=count, color = target))+  geom_line()+facet_wrap(~Location) +
  geom_point(data=mydf2, aes(x=date, y=count, color=target)) +facet_wrap(~Location) 
# again it doesn't have the delta wave in it. what does it matter if this peaks? 

# (8) concentrate,raw, norm . this is comparable to the C_N1_avg/norm of c pmmv spline 
mydf = bcww %>% mutate(N1 = c_n1_spl_norm, N2 = c_n2_spl_norm) %>% 
  select(Date_sampled, N1, N2, Location) %>% 
  pivot_longer(cols = c(2,3), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled) # splines
mydf2 = bcww %>% mutate(N1= C_N1_1/c_pmmv_avg_spl_norm,N2= C_N2_1/c_pmmv_avg_spl_norm) %>% 
  select(Date_sampled, N1, N2, Location) %>% 
  pivot_longer(cols = c(2,3), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled) # data 
ggplot(mydf, aes(x=date, y=count, color = target))+  geom_line()+facet_wrap(~Location) +
  geom_point(data=mydf2, aes(x=date, y=count, color=target)) +facet_wrap(~Location) 
# again - as in (4), misses the delta wave. peak time conclusions seem similar. 
# lower values of myspar make crazy things happen ... not sure why. 
# not sure whether i believe in splines with guesses for the smoothing parameter as any reasonable estimate
# for whether something has peaked yet or not. I think i don't believe them.
# however, the 'peak early Jan' is pretty robust from these 8 plots under this assumption. 


# get Jens to tell if we have 'location' based case data 

# ---- now make statistical models with selections of the variables ---- 
statbcww = filter(bcww,Location == "VLG") %>% select(Date_sampled,
                           S_N1_1, S_N1_2, S_N1_avg, 
                           C_N1_1, C_N1_2, C_N1_avg, 
                           S_PMMV_1, S_PMMV_2, 
                           C_PMMV_1, C_PMMV_2,
                           C_PMMV1.10_1, S_PMMV1.10_1,
                           C_PMMV1.10_2, S_PMMV1.10_2,  # %>% # same point about PMMV *not* avg 
                           temp, tss, rain_t,  rain_y, dinflvol) %>% 
  rename(date = Date_sampled) %>% left_join( select(vch, date, cases))

# regression model. Omit date as a predictor. use everything else. 

# NEXT UP -- NEED TO HANDLE LOCATION BETTER - i don't have cases for the regions
# may need both a ME model and a model with hidden variables c_location where sum(c) = cases
fitstop =  ymd("2021-10-01")
mod1 = lm(cases ~ ., data = filter(statbcww, date < fitstop)[,-1])
preds = predict(mod1, newdata = statbcww[,-1]) # predict for all the days 
ggplot(statbcww, aes(x=date, y=cases))+geom_point(color="blue", alpha = 0.5) + 
  geom_line(data = data.frame(date = statbcww$date, cases = preds)) + 
  annotate("rect", xmin = fitstop, xmax = max(statbcww$date), ymin = 0, ymax = 1400, 
            fill = "grey", alpha = 0.3) + theme_light()

# it is quite silly but I could do that separately for each location, and average the results 
allstats = lapply(unique(bcww$Location), function(x) {
  filter(bcww,Location == x) %>% select(Date_sampled, S_N1_1, S_N1_2, #S_N1_avg, 
                                        C_N1_1, C_N1_2, #C_N1_avg, 
                                        S_PMMV_1, S_PMMV_2, 
                                        C_PMMV_1, C_PMMV_2,
                                        C_PMMV1.10_1, S_PMMV1.10_1,
                                        C_PMMV1.10_2, S_PMMV1.10_2) %>% 
    rename(date = Date_sampled) %>% left_join( select(vch, date, cases))
})
allmods = lapply(allstats, function(x) { 
  lm(cases ~ ., data = filter(x, date < fitstop)[,-c(1,2)])})
allpreds = lapply(1:length(allmods), function(x) { 
  tmp = predict( allmods[[x]],  newdata = allstats[[x]][,-c(1,2)])
  return(data.frame(date = allstats[[x]]$date, pred = tmp))})
                  
tmp1 = merge(allpreds[[1]], allpreds[[2]], by = "date")
tmp2 = merge(allpreds[[3]], allpreds[[4]], by = "date")
tmp = merge(tmp2, allpreds[[5]], by = "date") 
colnames(tmp1) = c("date",unique(bcww$Location)[c(1,2)])
colnames(tmp) = c("date",unique(bcww$Location)[c(3,4,5)])
df = merge(tmp1, tmp, by="date")
df = df %>% mutate(avg = (VAI+VII+VLG+VLI+VNL )/5) %>% left_join(vch, by = "date")

ggplot(df, aes(x=date, y = avg))+geom_line() + geom_point()+
  geom_point(aes(x=date, y =cases), alpha=0.2, color="blue") + 
  annotate("rect", xmin = fitstop, xmax = max(statbcww$date), ymin = 0, ymax = 1250, 
           fill = "grey", alpha = 0.3) + theme_light()
# right. So - none of these statistical models from these locations, when fitted to the VCH case data, pick up 
# a large omicron wave of the kind we think we had . 


# ----- leftover plots 
ggplot(filter(bcww, Date_sampled < ymd("2021-12-01")), 
       aes(x=Date_sampled, y=pmmv_avg_spl,color=Location))+geom_line()+ facet_wrap(~Location)+
  geom_line(aes(x=Date_sampled, y=S_PMMV_avg,color=Location))

bcww = bcww %>% mutate( s_n1_norm = S_N1_avg / pmmv_avg_spl_norm) 
ggplot(filter(bcww, Date_sampled < ymd("2021-12-01")), 
       aes(x=Date_sampled, y=s_n1_norm,color=Location))+geom_line()+ facet_wrap(~Location)+
  geom_line(data = filter(vch,date >ymd("2021-01-30") & 
                            date < ymd("2021-12=01")) , 
            inherit.aes = F, aes(x=date, y = cases/5), color="black")

# the cp columns are really odd 
ggplot(bcww, aes(x=Date_sampled, y=as.numeric(S_N1_1_cp),color=Location))+
  geom_line()+facet_wrap(~Location) +
  geom_line(aes(x=Date_sampled, y = S_N1_1),color="gray")

# smooth and add all the comparable ww streams that i think should be related to case numbers
bcww = bcww %>% filter(Date_sampled > ymd("2021-02-05")) %>%  group_by(Location) %>%
  mutate(sn1_avg_spl = smooth.spline(S_N1_1_cp+ S_N1_2_cp, spar = 0.7)$y) %>% 
  mutate(sn1_avg_norm = sn1_avg_spl/ pmmv_avg_spl) %>% 
  mutate(sn2_avg_spl = smooth.spline(S_N2_1_cp + S_N2_2_cp, spar = 0.7)$y) %>%
  mutate(sn2_avg_norm = sn2_avg_spl/ pmmv_avg_spl) %>% 
  mutate(cn1_avg_spl = smooth.spline(C_N1_1_cp + C_N1_2_cp, spar = 0.7)$y) %>% 
  mutate(cn1_avg_norm = cn1_avg_spl/ pmmv_avg_spl) %>% 
  mutate(cn2_avg_spl = smooth.spline(C_N2_1_cp + C_N2_2_cp, spar = 0.7)$y) %>%
  mutate(cn2_avg_norm = cn2_avg_spl/ pmmv_avg_spl)

ggplot(bcww, aes(x=Date_sampled, y=sn1_avg_spl,color=Location))+
  geom_line()+facet_wrap(~Location) +
  geom_line(aes(x=Date_sampled, y =sn2_avg_spl))+
  geom_point(aes(x=Date_sampled, y = S_N1_1_cp+ S_N1_2_cp), alpha=0.5,color="black") + 
  geom_point(aes(x=Date_sampled, y = S_N2_1_cp+ S_N2_2_cp), alpha=0.5) 

mydf = bcww %>% mutate(N1=S_N1_1_cp + S_N1_2_cp, 
                       N2=S_N2_1_cp + S_N2_2_cp,
                       N1C = C_N1_1_cp + C_N1_2_cp,
                       N2C = C_N2_1_cp + C_N2_2_cp) %>% 
  select(Date_sampled, N1, N2, N1C, N2C, Location) %>% 
  pivot_longer(cols = c(2,3,4,5), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled)

mydf = bcww %>% mutate(N1=S_N1_1_cp + S_N1_2_cp, 
                       N2=S_N2_1_cp + S_N2_2_cp,
                       N1C = C_N1_1_cp + C_N1_2_cp,
                       N2C = C_N2_1_cp + C_N2_2_cp) %>% 
  mutate(n1_spl = smooth.spline(N1, spar = 0.7)$y,
         n2_spl = smooth.spline(N2, spar = 0.7)$y,
         n1c_spl = smooth.spline(N1C, spar = 0.7)$y,
         n2c_spl = smooth.spline(N2C, spar = 0.7)$y)

raws = mydf %>% 
  select(Date_sampled, N1, N2, N1C, N2C, Location) %>% 
  pivot_longer(cols = c(2,3,4,5), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled)
splines = mydf %>%  select(Date_sampled, n1_spl, n2_spl, n1c_spl, n2c_spl, Location) %>% 
  pivot_longer(cols = c(2,3,4,5), names_to = "target", values_to = "count") %>% 
  mutate(date = Date_sampled)


ggplot(filter(mydf, date> ymd("2021-09-01")), aes(x=date, y=count, color=target))+
  geom_line()+facet_wrap(~Location) +geom_point()

ggplot(filter(bcww, Date_sampled> ymd("2021-09-01")), 
       aes(x=Date_sampled, y = S_PMMV_avg))+
  geom_line()+facet_wrap(~Location)+geom_point()
