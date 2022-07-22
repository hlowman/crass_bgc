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

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

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

# #### get hrly avg ####
# 
# ## plot monthly avg from all freq vs hrly - should be similar
# # round time stamps to nearest hour
# pptz_h = pptz
# for(i in names(pptz_h)){
#   pptz_h[[i]]$datetimeMT = lubridate::round_date(pptz_h[[i]]$datetimeMT, "hour") 
#   pptz_h[[i]] = 
#     pptz_h[[i]] %>% 
#     group_by(datetimeMT,station) %>% 
#     dplyr::summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))
#   pptz_h[[i]][is.nan(pptz_h[[i]])]=NA
# }
# # get monthly avg from all freq
# pptz_allfreq_to_m = pptz
# for(i in names(pptz_allfreq_to_m)){
#   pptz_allfreq_to_m[[i]]$datetimeMT = lubridate::round_date(pptz_allfreq_to_m[[i]]$datetimeMT, "month") 
#   pptz_allfreq_to_m[[i]] = 
#     pptz_allfreq_to_m[[i]] %>% 
#     group_by(datetimeMT,station) %>% 
#     dplyr::summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))
#   pptz_allfreq_to_m[[i]][is.nan(pptz_allfreq_to_m[[i]])]=NA
# }
# # get monthly avg from hrly
# pptz_hrly_to_m = pptz_h
# for(i in names(pptz_hrly_to_m)){
#   pptz_hrly_to_m[[i]]$datetimeMT = lubridate::round_date(pptz_hrly_to_m[[i]]$datetimeMT, "month") 
#   pptz_hrly_to_m[[i]] = 
#     pptz_hrly_to_m[[i]] %>% 
#     group_by(datetimeMT,station) %>% 
#     dplyr::summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))
#   pptz_hrly_to_m[[i]][is.nan(pptz_hrly_to_m[[i]])]=NA
# }
# # compare
# dat1 = pptz_allfreq_to_m[["valletoledo_ppt"]]
# p1 = ggplot(data = dat1, aes(x=datetimeMT, y=ppt_in)) +
#   geom_col() + theme_classic() + scale_y_continuous(limits = c(0,0.005))
# dat2 = pptz_hrly_to_m[["valletoledo_ppt"]]
# p2 = ggplot(data = dat2, aes(x=datetimeMT, y=ppt_in)) +
#   geom_col() + theme_classic() + scale_y_continuous(limits = c(0,0.005))
# grid.arrange(p1,p2)
# 
# # save to individual dfs
# list2env(pptz_h,envir=.GlobalEnv)
# 
# # ## remove outliers routine
# # # make new list for corrected data
# # pptz_oc = pptz
# # # replace pts over temp treshold with NA, iterating through dataframes in list by site
# # for(i in names(pptz_oc)){
# #   dat = pptz_oc[[i]]
# #   ## add day column ##
# #   dat$day = lubridate::yday(dat$datetimeMT)
# #   ## add yr column ##
# #   dat$yr = lubridate::year(dat$datetimeMT)
# #   ## add day_yr column ##
# #   dat$day_yr = paste(dat$day, dat$yr, sep="_")
# #   ## get unique ids for each day
# #   dayz= unique(dat$day_yr)
# #   ## set temp threshold for each day ##
# #   # threshold to identify outlier points is set here to be x1.5 the upper (75%) quantile of ppt data on each date. Can edit by changing either quantile probability or what it is divided by. Can also add a lower threshold.
# #   daily.stats =
# #     dat[dat$day_yr %in% dayz,] %>%
# #     select(day_yr, ppt_in) %>%
# #     group_by(day_yr) %>%
# #     summarize_all(list(ppt_in.highQ = quantile), probs = 0.95, na.rm = TRUE) %>%
# #     #summarize_all(list(ppt_in.highQ = sd), na.rm = TRUE) %>%
# #     mutate_at(vars(-day_yr), ~(.*3))
# #   daily.stats$ppt_in.highQ[daily.stats$ppt_in.highQ<=0] = NA
# #   #hist(daily.stats$ppt_in.highQ)
# #   # join threshold to data #
# #   dat = left_join(dat, daily.stats, by="day_yr")
# #   # replace pts above ppt_in theshold on service dates with NA #
# #   dat$ppt_in[dat$ppt_in>dat$ppt_in.highQ & !is.na(dat$ppt_in.highQ)] = NA
# #   #dat[which(dat$ppt_in > dat$ppt_in.highQ), c(5:15,17)] = NA
# #   # replace in list
# #   pptz_oc[[i]] = dat
# # }
# 
# 

#### get daily ppt sums ####

# round time stamps to nearest day, sum ppt by day, get mean of all other numerics
pptz_d = pptz
for(i in names(pptz_d)){
  pptz_d[[i]] = pptz_d[[i]][,c(1,2,12,13)]
  pptz_d[[i]]$date = lubridate::round_date(pptz_d[[i]]$datetimeMT, "day") 
  pptz_d[[i]] = 
    pptz_d[[i]] %>% 
    group_by(date,station) %>% 
    #dplyr::summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))
    dplyr::summarise(across(ppt_in, sum))
              #across(wind_speed_mph:solar_rad_w.m2, mean))
  pptz_d[[i]][is.nan(pptz_d[[i]])]=NA
}

# save to individual dfs
list2env(pptz_d,envir=.GlobalEnv)

#### average across ppt stations for each watershed ####

# check avaiable date ranges in each site
range(redondo_ppt$date, na.rm = T)
range(hiddenvalley_ppt$date, na.rm = T)
range(losposos_ppt$date, na.rm = T)
range(sanantonio_ppt$date, na.rm = T)
range(vallegrande_ppt$date, na.rm = T)
range(valletoledo_ppt$date, na.rm = T)
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
ppt_RSAW$ID = "RSAW"

ppt_RSA = 
  bind_rows(valletoledo_ppt, losposos_ppt) %>%
  group_by(date) %>%   
  dplyr::summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))
ppt_RSA$Grab_bgc_name = "San Antonio Creek - Toledo"
ppt_RSA$ID = "RSA"
ppt_RSA = tibble::add_column(ppt_RSA, station = "avg_valletoledo_losposos", .after = "date")
ppt_RSA[is.nan(ppt_RSA)]=NA

ppt_SULF = 
  bind_rows(redondo_ppt, sanantonio_ppt) %>%
  group_by(date) %>%   
  dplyr::summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))
ppt_SULF$Grab_bgc_name = "Sulfur Creek"
ppt_SULF$ID = "SULF"
ppt_SULF = tibble::add_column(ppt_SULF, station = "avg_redondo_sanantonio", .after = "date")
ppt_SULF[is.nan(ppt_SULF)]=NA

ppt_EFJ = 
  bind_rows(redondo_ppt, losposos_ppt, vallegrande_ppt) %>%
  group_by(date) %>%   
  dplyr::summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))
ppt_EFJ$Grab_bgc_name = "East Fork Jemez River"
ppt_EFJ$ID = "EFJ"
ppt_EFJ = tibble::add_column(ppt_EFJ, station = "avg_losposos_redondo_vallegrande", .after = "date")
ppt_EFJ[is.nan(ppt_EFJ)]=NA

## make df of all
# make columns match
cols_to_keep <- intersect(colnames(ppt_RED), colnames(ppt_EFJ))
ppt_IND_BB <- ppt_IND_BB[,cols_to_keep, drop=FALSE]
ppt_IND <- ppt_IND[,cols_to_keep, drop=FALSE]
ppt_IND_AB <- ppt_IND_AB[,cols_to_keep, drop=FALSE]
ppt_RSAW <- ppt_RSAW[,cols_to_keep, drop=FALSE]
ppt_RSA <- ppt_RSA[,cols_to_keep, drop=FALSE]
ppt_SULF <- ppt_SULF[,cols_to_keep, drop=FALSE]
ppt_EFJ <- ppt_EFJ[,cols_to_keep, drop=FALSE]
ppt_all = rbind(ppt_IND_BB,ppt_IND,ppt_IND_AB,ppt_RED,ppt_RSAW,ppt_RSA,ppt_SULF,ppt_EFJ)

# pptz_h_ws = list(ppt_IND_BB,ppt_IND,ppt_IND_AB,ppt_RSAW,ppt_RSA,ppt_SULF,ppt_EFJ,ppt_all)
# names(pptz_h_ws) = c("ppt_IND_BB","ppt_IND","ppt_IND_AB","ppt_RSAW","ppt_RSA","ppt_SULF","ppt_EFJ","ppt_all")

pptz_d_ws = list(ppt_IND_BB,ppt_IND,ppt_IND_AB,ppt_RSAW,ppt_RSA,ppt_SULF,ppt_EFJ,ppt_all)
names(pptz_d_ws) = c("ppt_IND_BB","ppt_IND","ppt_IND_AB","ppt_RSAW","ppt_RSA","ppt_SULF","ppt_EFJ","ppt_all")

#
#### get monthly cumulative (sum) ####

ppt_all[which(is.na(ppt_all$ppt_in)),]
# 6 NAs... don't seem worth worrying about for calculating sums

pptz_m = pptz_d_ws
for(i in names(pptz_m)){
  pptz_m[[i]]$month = lubridate::month(pptz_m[[i]]$date) 
  pptz_m[[i]]$year = lubridate::year(pptz_m[[i]]$date) 
  pptz_m[[i]] = 
    pptz_m[[i]] %>% 
    group_by(month, year,ID) %>% 
    dplyr::summarize(cumulative_precip_in = sum(ppt_in, na.rm = TRUE)) %>%
    ungroup()
  pptz_m[[i]][is.nan(pptz_m[[i]])]=NA
  # convert inches to mm
  pptz_m[[i]]$cumulative_precip_mm = pptz_m[[i]]$cumulative_precip_in*25.4
}

# pull out df of all sites
ppt_cum_m = pptz_m[["ppt_all"]]

# plot 
ppt_cum_m_p =
  ggplot(data = ppt_cum_m, aes(x=month, y=cumulative_precip_mm)) +
  geom_bar(stat="identity") +
  facet_wrap(~ID)

# save df of all sites
# saveRDS(ppt_cum_m, "data_working/VCNPprecip_m_cum_edited_110321.rds")
saveRDS(ppt_cum_m, "data_working/VCNPprecip_m_cum_edited_20220721.rds")

#### get daily average ####

ppt_all[which(is.na(ppt_all$ppt_in)),]
# 4 NAs... don't seem worth worrying about for calculating sums

pptz_d = pptz_h_ws
for(i in names(pptz_d)){
  pptz_d[[i]]$datetimeMT = lubridate::round_date(pptz_d[[i]]$datetimeMT, "day") 
  pptz_d[[i]] = 
    pptz_d[[i]] %>% 
    group_by(datetimeMT,ID) %>% 
    dplyr::summarize(cumulative_precip_in = mean(ppt_in, na.rm = TRUE)) %>%
    ungroup()
  pptz_d[[i]][is.nan(pptz_d[[i]])]=NA
  # convert inches to mm
  pptz_d[[i]]$cumulative_precip_mm = pptz_d[[i]]$cumulative_precip_in*25.4
}

# pull out df of all sites
ppt_avg_d = pptz_d[["ppt_all"]]

# plot 
ppt_avg_dp =
  ggplot(data = ppt_avg_d, aes(x=datetimeMT, y=cumulative_precip_mm)) +
  geom_bar(stat="identity") +
  facet_wrap(~ID)

# save df of all sites
saveRDS(ppt_avg_d, "data_working/VCNPprecip_d_avg_edited_110321.rds")

#### get monthly average ####

ppt_all[which(is.na(ppt_all$ppt_in)),]
# 4 NAs... don't seem worth worrying about for calculating sums

pptz_m_avg = pptz_h_ws
for(i in names(pptz_m_avg)){
  pptz_m_avg[[i]]$datetimeMT = lubridate::round_date(pptz_m_avg[[i]]$datetimeMT, "month") 
  pptz_m_avg[[i]] = 
    pptz_m_avg[[i]] %>% 
    group_by(datetimeMT,ID) %>% 
    dplyr::summarize(cumulative_precip_in = mean(ppt_in, na.rm = TRUE)) %>%
    ungroup()
  pptz_m_avg[[i]][is.nan(pptz_m_avg[[i]])]=NA
  # convert inches to mm
  pptz_m_avg[[i]]$cumulative_precip_mm = pptz_m_avg[[i]]$cumulative_precip_in*25.4
}

# pull out df of all sites
ppt_avg_m = pptz_m_avg[["ppt_all"]]

# plot 
ppt_avg_m_p =
  ggplot(data = ppt_avg_m, aes(x=datetimeMT, y=cumulative_precip_mm)) +
  geom_bar(stat="identity") +
  facet_wrap(~ID)

# save df of all sites
saveRDS(ppt_avg_m, "data_working/VCNPprecip_m_avg_edited_110321.rds")
