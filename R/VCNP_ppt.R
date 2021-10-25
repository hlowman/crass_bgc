## Explore NM precipitation data ##

#### read me ####
# CRASS_bgc_ws	              ppt_station	ppt_station_2	ppt_station_3
# Indios Creek - above burn	  valletoledo	NA	NA
# Indios Creek	              valletoledo	NA	NA
# San Antonio Creek - Toledo	valletoledo	losposos	NA
# East Fork-Jemez River	      losposos	  redondo	vallegrande
# Redondo Creek	              redondo	    NA	NA
# Sulfur Creek	              redondo	    sanantonio	NA
# San Antonio - West	        sanantonio	 NA	NA

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

pptz = list(redondo_ppt, hiddenvalley_ppt, losposos_ppt, sanantonio_ppt, vallegrande_ppt, valletoledo_ppt, ppt_all)
names(pptz) = c("redondo_ppt", "hiddenvalley_ppt", "losposos_ppt", "sanantonio_ppt", "vallegrande_ppt", "valletoledo_ppt", "ppt_all")

for(i in names(pptz)){
  pptz[[i]][["datetime"]] = as.character(pptz[[i]][["datetime"]])
  pptz[[i]][["datetimeMT"]] = as.POSIXct(pptz[[i]][["datetime"]], "%Y%m%d%H%M", tz="US/Mountain")
}

list2env(pptz,envir=.GlobalEnv)

#### clean values ####

for(i in names(pptz)){
  pptz[[i]][pptz[[i]]==-9999]=NA
}

list2env(pptz,envir=.GlobalEnv)



#### plot ####

range(ppt_all$ppt_in, na.rm = T)
range(ppt_all$datetimeMT, na.rm = T)

ppt_sites = c("redondo_ppt", "hiddenvalley_ppt", "losposos_ppt", "sanantonio_ppt", "vallegrande_ppt", "valletoledo_ppt")
pdf("data_raw/NM_precip/plot_all_stations.pdf", width=12, height=12) 
par(mfrow=c(3,2))
for(i in ppt_sites){
  barplot(height = pptz[[i]][["ppt_in"]], names.arg = pptz[[i]][["datetimeMT"]],
          xlab="Date", ylab="Precipitation (in)", xaxt = "n")
  # xlim=c(as.Date("2003-10-24"), as.Date("2021-10-23"))  
  axis(1, pptz[[i]][["datetimeMT"]], format(pptz[[i]][["datetimeMT"]], "%m %Y"), cex.axis = .7)
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
