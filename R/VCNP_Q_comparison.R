#### read me ####

# the purpose of this script is to compare available discharge records for watersheds in the Valles Caldera, NM
# this is only possible for the NM sites

#### libraries ####

library(tidyverse)

#### DP - load data ####

VALL_EFJR_HV_DP = read.csv(file = "data_raw/Discharge/from_DaveP/JemezHV_Discharge.csv", skip = 14)
names(VALL_EFJR_HV_DP) = c("datetime_UTC","datetime_UTC0600","Q_cfs","approval","grade","qualifiers")
table(VALL_EFJR_HV_DP$approval)
table(VALL_EFJR_HV_DP$grade)
table(VALL_EFJR_HV_DP$qualifiers)

VALL_RSA_VT_DP = read.csv(file = "data_raw/Discharge/from_DaveP/RSA_Toledo_Discharge.csv", skip = 14)
names(VALL_RSA_VT_DP) = c("datetime_UTC","datetime_UTC0700","Q_cfs","approval","grade","qualifiers")
table(VALL_RSA_VT_DP$approval)
table(VALL_RSA_VT_DP$grade)
table(VALL_RSA_VT_DP$qualifiers)

VALL_RSA_W_DP = read.csv(file = "data_raw/Discharge/from_DaveP/RSA_West_Discharge.csv", skip = 14)
names(VALL_RSA_W_DP) = c("datetime_UTC","datetime_UTC0700","Q_cfs","approval","grade","qualifiers")
table(VALL_RSA_W_DP$approval)
table(VALL_RSA_W_DP$grade)
table(VALL_RSA_W_DP$qualifiers)

#### DP - format date/time and join sites ####

VALL_EFJR_HV_DP$datetime_UTC = gsub(pattern = "T", replacement = " ", VALL_EFJR_HV_DP$datetime_UTC)
VALL_EFJR_HV_DP$datetime_UTC = gsub(pattern = "Z", replacement = "", VALL_EFJR_HV_DP$datetime_UTC)
VALL_EFJR_HV_DP$datetime_UTC = as.POSIXct(VALL_EFJR_HV_DP$datetime_UTC,, format = "%Y-%m-%d %H:%M:%S", tz="UTC")
VALL_EFJR_HV_DP$datetime_MT = as.POSIXct(format(VALL_EFJR_HV_DP$datetime_UTC, tz="US/Mountain"))
VALL_EFJR_HV_DP$site = "VALL_EFJR_HV"

VALL_RSA_VT_DP$datetime_UTC = gsub(pattern = "T", replacement = " ", VALL_RSA_VT_DP$datetime_UTC)
VALL_RSA_VT_DP$datetime_UTC = gsub(pattern = "Z", replacement = "", VALL_RSA_VT_DP$datetime_UTC)
VALL_RSA_VT_DP$datetime_UTC = as.POSIXct(VALL_RSA_VT_DP$datetime_UTC,, format = "%Y-%m-%d %H:%M:%S", tz="UTC")
VALL_RSA_VT_DP$datetime_MT = as.POSIXct(format(VALL_RSA_VT_DP$datetime_UTC, tz="US/Mountain"))
VALL_RSA_VT_DP$site = "VALL_RSA_VT"

VALL_RSA_W_DP$datetime_UTC = gsub(pattern = "T", replacement = " ", VALL_RSA_W_DP$datetime_UTC)
VALL_RSA_W_DP$datetime_UTC = gsub(pattern = "Z", replacement = "", VALL_RSA_W_DP$datetime_UTC)
VALL_RSA_W_DP$datetime_UTC = as.POSIXct(VALL_RSA_W_DP$datetime_UTC,, format = "%Y-%m-%d %H:%M:%S", tz="UTC")
VALL_RSA_W_DP$datetime_MT = as.POSIXct(format(VALL_RSA_W_DP$datetime_UTC, tz="US/Mountain"))
VALL_RSA_W_DP$site = "VALL_RSA_W"

Q_DP = rbind(VALL_EFJR_HV_DP[,c(1,3,7,8)],VALL_RSA_VT_DP[,c(1,3,7,8)],VALL_RSA_W_DP[,c(1,3,7,8)])

#### DP - plot data ####

# ggplot(Q_DP, aes(x=datetime_MT, y=Q_cfs, fill=site)) + 
#   geom_path() +
#   theme_classic()

par(mfrow=c(3,1))
for( i in c("VALL_EFJR_HV", "VALL_RSA_VT", "VALL_RSA_W")){
  temp = Q_DP[Q_DP$site==i,]
  plot(x=temp$datetime_MT, y=temp$Q_cfs, type="l")
}

# plot(x=Q_DP$datetime_MT[Q_DP$site=="VALL_RSA_VT"], 
#      y=Q_DP$Q_cfs[Q_DP$site=="VALL_RSA_VT"], 
#      type="l")

#### BS - load data ####

VALL_EFJR_BS = read.csv("data_raw/Discharge/from_BetsyS/EFJ_2009_2013.csv")
VALL_EFJR_BS = VALL_EFJR_BS[,1:4]

#### BS - format date/time ####

VALL_EFJR_BS$datetime_MT = as.POSIXct(VALL_EFJR_BS$date_time, format="%m/%d/%y %H:%M", tz="US/Mountain")


#### USGS - load data ####



#### plot all data ####

par(mfrow=c(4,1))

plot(x=VALL_EFJR_BS$datetime_MT, y=VALL_EFJR_BS$Q_cfs, 
     type="l", xlab="", ylab="Q (csf)", ylim=c(0,420), 
     xlim=c(as.POSIXct("2008-01-01"),as.POSIXct("2021-01-01")))
title(main = "VALL_EFJR from PT, Betsy S.")

for( i in c("VALL_EFJR_HV", "VALL_RSA_VT", "VALL_RSA_W")){
  temp = Q_DP[Q_DP$site==i,]
  plot(x=temp$datetime_MT, y=temp$Q_cfs, 
       type="l", xlab="", ylab="Q (csf)", ylim=c(0,420), 
       xlim=c(as.POSIXct("2008-01-01"),as.POSIXct("2021-01-01")))
  title(main = paste(i, "from PT, Dave P."))
}

par(mfrow=c(1,1))
plot(x=VALL_EFJR_BS$datetime_MT, y=VALL_EFJR_BS$Q_cfs, 
     type="l", xlab="", ylab="Q (csf)", ylim=c(0,150), 
     xlim=c(as.POSIXct("2012-08-01"),as.POSIXct("2013-10-01")),
     col="blue", lwd=3)
temp = Q_DP[Q_DP$site=="VALL_EFJR_HV",]
lines(x=temp$datetime_MT, y=temp$Q_cfs, 
     type="l", xlab="", ylab="Q (csf)", ylim=c(0,150), 
     xlim=c(as.POSIXct("2012-08-01"),as.POSIXct("2013-10-01")),
     col="red")

