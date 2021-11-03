## Organize NM precipitation data ##
# AJW

#### read me ####
# the purpose of this script is to organize and clean VCNP precip data to get a daily average and monthly cumulative record for other analyses 

# all data was retreived from here: https://wrcc.dri.edu/vallescaldera/

#### libraries ####
library(tidyverse)
library(lubridate)
library(ggplot2)
library(scales)
library(plyr)
library(gridExtra)
library(zoo)
library(xts)
library(forecast)

#### load ppt data ####

header = c("datetime","ppt_in","wind_speed_mph","wind_dir_deg","air_temp_F","fuel_temp_F","rel_humidity_perc","battery_volts","maxgust_dir_deg","maxgust_speed_mph","solar_rad_w.m2")

redondo_ppt = read.table("data_raw/NM_precip/redondo_ppt_20040826_20211021.txt", skip=4, header=F, fill = T)
colnames(redondo_ppt) = header
redondo_ppt$station = "redondo"

hiddenvalley_ppt = read.table("data_raw/NM_precip/hiddenvalley_ppt_20110518_20211022.txt", skip=4, header=F, fill = T)
colnames(hiddenvalley_ppt) = header
hiddenvalley_ppt$station = "hiddenvalley"

losposos_ppt = read.table("data_raw/NM_precip/losposos_ppt_20040715_20211022.txt", skip=4, header=F, fill = T)
colnames(losposos_ppt) = header
losposos_ppt$station = "losposos"

sanantonio_ppt = read.table("data_raw/NM_precip/sanantonio_ppt_20040714_20211022.txt", skip=4, header=F, fill = T)
colnames(sanantonio_ppt) = header
sanantonio_ppt$station = "sanantonio"

vallegrande_ppt = read.table("data_raw/NM_precip/vallegrande_ppt_20031025_20211022.txt", skip=4, header=F, fill = T)
colnames(vallegrande_ppt) = header
vallegrande_ppt$station = "vallegrande"

valletoledo_ppt = read.table("data_raw/NM_precip/valletoledo_ppt_20050603_20211022.txt", skip=4, header=F, fill = T)
colnames(valletoledo_ppt) = header
valletoledo_ppt$station = "valletoledo"

ppt_all = rbind(redondo_ppt, hiddenvalley_ppt, losposos_ppt, sanantonio_ppt, vallegrande_ppt, valletoledo_ppt)

#### format datetime ####

# date/time is reported in header as "Hour of Day Ending is L.S.T.", which I ebelive means that it is in local standard time (LST = local standard time). I will therefore convert all dates/times to "Etc/GMT+7" which is Mountain Standard Time (winter) or UTC-7.
# Note that when I convert to US/Mountain, there are many NA date/time values, I beleive because this tz shifts for daylight savings time automatically, which is not accounted for in the dataset. When I use "Etc/GMT+7", there are no NAs

pptz = list(redondo_ppt, hiddenvalley_ppt, losposos_ppt, sanantonio_ppt, vallegrande_ppt, valletoledo_ppt, ppt_all)
names(pptz) = c("redondo_ppt", "hiddenvalley_ppt", "losposos_ppt", "sanantonio_ppt", "vallegrande_ppt", "valletoledo_ppt", "ppt_all")

for(i in names(pptz)){
  pptz[[i]][["datetime"]] = as.character(pptz[[i]][["datetime"]])
  pptz[[i]][["datetimeMT"]] = as.POSIXct(pptz[[i]][["datetime"]], "%Y%m%d%H%M", tz="Etc/GMT+7")
  # data 
}

list2env(pptz,envir=.GlobalEnv)

#### clean up values ####

# replace -9999 missing values with NA
for(i in names(pptz)){
  pptz[[i]][pptz[[i]]==-9999]=NA
}
list2env(pptz,envir=.GlobalEnv)


## remove unrealistically high per hour points and replace with NA
for(i in names(pptz)){
  pptz[[i]]$ppt_in[pptz[[i]]$ppt_in>3] = NA
}
list2env(pptz,envir=.GlobalEnv)


## round time stamps to nearest hour
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
## remove duplicates - consolidate duplicates by finding average
pptz_h = pptz
for(i in names(pptz_h)){
  pptz_h[[i]]$datetimeMT = lubridate::round_date(pptz_h[[i]]$datetimeMT, "60 minutes") 
  pptz_h[[i]] = 
    pptz_h[[i]] %>% 
    group_by(datetimeMT,station) %>% 
    dplyr::summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))
  pptz_h[[i]][is.nan(pptz_h[[i]])]=NA
}
# plot to check hrly against raw data
ppt_sites = c("redondo_ppt")
pdf("figures/NM_ppt_rawVShr_comp.pdf", width=20, height=10) 
par(mfrow=c(2,1))
for(i in ppt_sites){
  barplot(height = pptz[[i]][["ppt_in"]], names.arg = pptz[[i]][["datetimeMT"]],
          xlab="Date", ylab="Precipitation (in)", main=i, ylim = c(0,1.5))
  barplot(height = pptz_h[[i]][["ppt_in"]], names.arg = pptz_h[[i]][["datetimeMT"]],
          xlab="Date", ylab="Precipitation (in)", main=i, xaxt = "n")
  # xlim=c(as.Date("2003-10-24"), as.Date("2021-10-23"))  
  #axis.POSIXct(1, at=seq(min(pptz_h[[i]][["datetimeMT"]]), max(pptz_h[[i]][["datetimeMT"]]),  by="years"), format="%m-%Y")
  axis.POSIXct(1, at=pptz_h[[i]][["datetimeMT"]], labels=format(pptz_h[[i]][["datetimeMT"]], "%m/%d"))
}
dev.off() 
# save to individual dfs
list2env(pptz_h,envir=.GlobalEnv)

dat = pptz_h[[i]]
dat = dat[dat$datetimeMT>as.Date("2005-01-01") & dat$datetimeMT<as.Date("2005-02-01"),]
ggplot(data = dat, aes(x=datetimeMT, y=ppt_in)) +
  geom_col() + theme_classic()

# ## remove outliers routine
# # make new list for corrected data
# pptz_oc = pptz
# # replace pts over temp treshold with NA, iterating through dataframes in list by site
# for(i in names(pptz_oc)){
#   dat = pptz_oc[[i]]
#   ## add day column ##
#   dat$day = lubridate::yday(dat$datetimeMT)
#   ## add yr column ##
#   dat$yr = lubridate::year(dat$datetimeMT)
#   ## add day_yr column ##
#   dat$day_yr = paste(dat$day, dat$yr, sep="_")
#   ## get unique ids for each day
#   dayz= unique(dat$day_yr)
#   ## set temp threshold for each day ##
#   # threshold to identify outlier points is set here to be x1.5 the upper (75%) quantile of ppt data on each date. Can edit by changing either quantile probability or what it is divided by. Can also add a lower threshold.
#   daily.stats =
#     dat[dat$day_yr %in% dayz,] %>%
#     select(day_yr, ppt_in) %>%
#     group_by(day_yr) %>%
#     summarize_all(list(ppt_in.highQ = quantile), probs = 0.95, na.rm = TRUE) %>%
#     #summarize_all(list(ppt_in.highQ = sd), na.rm = TRUE) %>%
#     mutate_at(vars(-day_yr), ~(.*3))
#   daily.stats$ppt_in.highQ[daily.stats$ppt_in.highQ<=0] = NA
#   #hist(daily.stats$ppt_in.highQ)
#   # join threshold to data #
#   dat = left_join(dat, daily.stats, by="day_yr")
#   # replace pts above ppt_in theshold on service dates with NA #
#   dat$ppt_in[dat$ppt_in>dat$ppt_in.highQ & !is.na(dat$ppt_in.highQ)] = NA
#   #dat[which(dat$ppt_in > dat$ppt_in.highQ), c(5:15,17)] = NA
#   # replace in list
#   pptz_oc[[i]] = dat
# }


#### plot ####

range(ppt_all$ppt_in, na.rm = T)
range(ppt_all$datetimeMT)

ppt_sites = c("redondo_ppt", "hiddenvalley_ppt", "losposos_ppt", "sanantonio_ppt", "vallegrande_ppt", "valletoledo_ppt")
pdf("figures/NM_ppt_all_stations.pdf", width=20, height=12) 
par(mfrow=c(3,2))
for(i in ppt_sites){
  barplot(height = pptz[[i]][["ppt_in"]], names.arg = pptz[[i]][["datetimeMT"]],
          xlab="Date", ylab="Precipitation (in)", main=i)
  # xlim=c(as.Date("2003-10-24"), as.Date("2021-10-23"))  , xaxt = "n"
  #axis(1, pptz[[i]][["datetimeMT"]], format(pptz[[i]][["datetimeMT"]], "%m %Y"), cex.axis = .7)
}
dev.off() 
# 
# ppt_all_p =
#   ggplot(data = pptz[["ppt_all"]], aes(x=datetimeMT, y=ppt_in)) +
#   # geom_point(aes(colour=station)) +
#   # geom_line(aes(colour=station))+
#   geom_bar(stat="identity") +
#   facet_wrap(~station)
#   #scale_color_brewer(palette="RdBu", direction = -1, name="Station")

# # plot each
# par(mfrow=c(1,1))
# barplot(height = pptz[[1]][["ppt_in"]], names.arg = pptz[[1]][["datetimeMT"]],
#         xlab="Date", ylab="Precipitation (in)",main = names(pptz)[1])
# barplot(height = pptz[[1]][["ppt_in"]], names.arg = pptz[[1]][["datetimeMT"]],
#         xlab="Date", ylab="Precipitation (in)",main = names(pptz)[1])

#### average across sites ####

# check avaiable date ranges in each site
range(redondo_ppt$datetimeMT, na.rm = T)
range(hiddenvalley_ppt$datetimeMT, na.rm = T)
range(losposos_ppt$datetimeMT, na.rm = T)
range(sanantonio_ppt$datetimeMT, na.rm = T)
range(vallegrande_ppt$datetimeMT, na.rm = T)
range(valletoledo_ppt$datetimeMT, na.rm = T)
# there is significantly less data in hidden valley, so don't use that site if possible

# HOW TO AVERAGE:
# ID	    Grab_bgc_name	                            ppt_station	ppt_station_2	ppt_station_3
# IND_BB	Indios Creek - Post Fire (Below Burn)	    valletoledo	NA	          NA
# IND	    Indios Creek	                            valletoledo	NA	          NA
# IND_AB	Indios Creek - above burn	                valletoledo	NA	          NA
# RED	    Redondo Creek	                            redondo	    NA	          NA
# RSAW	  San Antonio - West	                      sanantonio	NA	          NA

# RSA	    San Antonio Creek - Toledo	              valletoledo	losposos	    NA
# SULF	  Sulfur Creek                            	redondo	    sanantonio	  NA
# EFJ	    East Fork Jemez River	                    losposos	  redondo	      vallegrande


ppt_IND_BB = valletoledo_ppt
ppt_IND_BB$Grab_bgc_name = "Indios Creek - Post Fire (Below Burn)"
ppt_IND_BB$ID = "IND_BB"

ppt_IND = valletoledo_ppt
ppt_IND$Grab_bgc_name = "Indios Creek"
ppt_IND$ID = "IND"

ppt_IND_AB = valletoledo_ppt
ppt_IND_AB$Grab_bgc_name = "Indios Creek - above burn"
ppt_IND_AB$ID = "IND_AB"

ppt_RED = redondo_ppt
ppt_RED$Grab_bgc_name = "Redondo Creek"
ppt_RED$ID = "RED"

ppt_RSAW = sanantonio_ppt
ppt_RSAW$Grab_bgc_name = "San Antonio - West"
ppt_RED$ID = "RSAW"

## create clean set of time stamps to join to ##
time <- data.frame(
  datetimeMT = seq.POSIXt(
    from = ISOdatetime(2004,07,15,12,0,0, tz = "Etc/GMT+7"),
    to = ISOdatetime(2021,10,22,15,00,0, tz= "Etc/GMT+7"),
    by = "15 min" ))

