# MARSS model @ Santa Barbara Coastal LTER Sites
# October 24, 2021
# Heili Lowman, Alex Webster

# The following script will run a MARSS analysis at SBC LTER sites for the CRASS project.
# M - Multi-variate
# AR - Auto-regressive
# S - State
# S - Space

# Note: Following a discussion with John Melack, the following sites have both enough chemistry
# and precip data, as well as fires that occurred within that timeframe to be included in the
# MARSS analysis performed below.

# Arroyo Burro - AB00
# Atascadero - AT07
# Gaviota - GV01
# Arroyo Hondo - HO00
# Mission Creek (at Rocky Nook) - MC06
# Refugio - RG01
# Rattlesnake - RS02
# San Pedro - SP02

# The following site may be added later, but a longer precipitation record needs to be
# identified for them - Bell Canyon (BC02).

#### Setup ####

# Load packages
library(tidyverse)
library(lubridate)
library(here)
library(MARSS)

# load fxn to replace NaNs with NAs
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))

# Load datasets
# Stream Chemistry - all sites
chem <- readRDS("data_working/SBchem_edited_110721.rds")
# Precipitation - all sites
precip <- readRDS("data_working/SBprecip_edited_110721.rds")
# Fire Events - all sites
fire <- readRDS("data_working/SBfire_edited_111721.rds")
# Site Location information
location <- read_csv("data_raw/sbc_sites_stream_hydro.csv")

#### Filter & Joining Data ####

# I first need to identify the matching stream sites for the precip data
precip_ed <- precip %>%
  mutate(sitecode_match = factor(case_when(
    sitecode == "GV202" ~ "GV01",
    sitecode == "BARA" ~ "HO00", # HO201 doesn't start until 2002
    sitecode == "RG202" ~ "RG01",
    sitecode == "SMPA" ~ "SP02",
    sitecode == "GORY" ~ "AT07",
    sitecode == "CAWTP" ~ "AB00",
    sitecode == "STFS" ~ "MC06", # BOGA doesn't start until 2005
    sitecode == "ELDE" ~ "RS02"))) %>%
  dplyr::rename(site_precip = site,
         sitecode_precip = sitecode) %>%
  mutate(Day = 1) %>% # new column of "days"
  mutate(Date = make_date(Year, Month, Day))

## Timeframe Selection

# Examine precip data coverage for MARSS timescale delineation
precip_ed %>%
  drop_na(sitecode_match) %>%
  ggplot(aes(x = Year, y = sitecode_match, color = sitecode_match)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none")

# So, covariate data, which cannot be missing, can run from a maximum of 9/2002 to 7/2016.

# To avoid strange gaps in data, I'm going to start with the fire data,
# since I know it extends from 9/1/2002 to 7/1/2016. This will ensure all other datasets are 
# joined to these dates in full (which was causing problems earlier).
fire_precip <- left_join(fire, precip_ed, by = c("site" = "sitecode_match", "date" = "Date"))

# So, for some reason we're missing 1 years data (11/2013-9/2014) for the AB00, MC06, and RS02 
# precip sites, so I'm going to fill them in as zeros for now. Eventually, I will go back to the 
# county website and fill them in by hand. (https://www.countyofsb.org/pwd/dailyrain.sbc)
fire_precip <- fire_precip %>%
  mutate_at(vars(cumulative_precip_mm), ~replace(., is.na(.), 0)) %>%
  mutate(year = year(date),
         month = month(date)) # and the Year/Month didn't populate, so adding in new columns

# Then, left join with chemistry so as not to lose any data.
dat <- left_join(fire_precip, chem, by = c("site", "year" = "Year", "month" = "Month"))

# Adding in dummy covariates by season
n_months <- dat$date %>%
  unique() %>%
  length()

seas_1 <- sin(2 * pi * seq(n_months) / 12)
seas_2 <- cos(2 * pi * seq(n_months) / 12)

dat <- dat %>%
  mutate(Season1 = rep(seas_1, 8),
         Season2 = rep(seas_2, 8))

# AJW: replace NaNs with NAs
dat[is.nan(dat)] = NA

# And, inspect dataset for missing covariate data.
sum(is.na(dat$fire)) # 0
sum(is.na(dat$cumulative_precip_mm)) # 0
# Great!

# Note for future me - be VERY careful with the joining above. Something weird was happening previously where precip data that IS present was simply dropping off.

#### Model fit ####

# Data: Stream Chemistry analytes (NH4, NO3, TDN, TPN, PO4, TDP, TP, TPP, TPC, TSS, SCond)
# Covariates: Month, Precip, Fire

# Starting with NH4 for test run.

# Note: Not scaling for now, but this should also be added in later.
dat_nh4 <- dat %>%
  select(site, year, month, Season1, Season2, mean_nh4_uM, cumulative_precip_mm, AB00_Tea, AB00_Jesusita, AT07_Jesusita, GV01_Gaviota, HO00_Gaviota, HO00_Sherpa, MC06_Tea, MC06_Jesusita, RG01_Gaviota, RG01_Sherpa, RS02_Tea, RS02_Jesusita, SP02_Gap) %>% # keeping year and month in so I don't get the "duplicates arise" error
  pivot_wider(names_from = site, values_from = c(mean_nh4_uM, cumulative_precip_mm, AB00_Tea, AB00_Jesusita, AT07_Jesusita, GV01_Gaviota, HO00_Gaviota, HO00_Sherpa, MC06_Tea, MC06_Jesusita, RG01_Gaviota, RG01_Sherpa, RS02_Tea, RS02_Jesusita, SP02_Gap)) %>%
  select(year, month, Season1, Season2, 
         mean_nh4_uM_AB00, mean_nh4_uM_AT07, mean_nh4_uM_GV01, mean_nh4_uM_HO00, mean_nh4_uM_MC06, mean_nh4_uM_RG01, mean_nh4_uM_RS02, mean_nh4_uM_SP02,
         cumulative_precip_mm_AB00, cumulative_precip_mm_AT07, cumulative_precip_mm_GV01, cumulative_precip_mm_HO00, cumulative_precip_mm_MC06, cumulative_precip_mm_RG01, cumulative_precip_mm_RS02, cumulative_precip_mm_SP02,
         AB00_Tea_AB00, AB00_Jesusita_AB00, AT07_Jesusita_AT07, GV01_Gaviota_GV01, HO00_Gaviota_HO00, HO00_Sherpa_HO00, MC06_Tea_MC06, MC06_Jesusita_MC06, RG01_Gaviota_RG01, RG01_Sherpa_RG01, RS02_Tea_RS02, RS02_Jesusita_RS02, SP02_Gap_SP02)

par(mfrow=c(2,1))
plot(dat_nh4$mean_nh4_uM_HO00, type="b")
plot(dat_nh4$mean_nh4_uM_RG01, type="b")
par(mfrow=c(1,1))
  
# Pull out only NH4 data
dat_dep <- t(dat_nh4[,5:12])

# Make covariate inputs
dat_cov <- dat_nh4[,c(3:4,13:33)]
dat_cov <- t(scale(dat_cov))

# Replace NaNs with 0, which don't play nice with the scale()
# dat_cov[is.nan(dat_cov)] <- 0

# make C matrix
# this matrix controls what covars predict what response vars; in contrast to...
# unconstrained: seperately estimates ALL correlations of predictor vars with x, i.e. how predictor vars drive ts dynamics. I think this structure is ok given that covars are unique to each site, but it would not be ok if there is a mix of shared and unique covars (are year and month used? bc if so, these are shared)
CC <- matrix(list("Season1", "Season1", "Season1", "Season1", "Season1", "Season1", "Season1", "Season1",
                  "Season2", "Season2", "Season2", "Season2", "Season2", "Season2", "Season2", "Season2", 
                  "AB00_precip", 0, 0, 0, 0, 0, 0, 0,
                  0, "AT07_precip", 0, 0, 0, 0, 0, 0,
                  0, 0, "GV01_precip", 0, 0, 0, 0, 0,
                  0, 0, 0, "HO00_precip", 0, 0, 0, 0,
                  0, 0, 0, 0, "MC06_precip", 0, 0, 0,
                  0, 0, 0, 0, 0, "RG01_precip", 0, 0,
                  0, 0, 0, 0, 0, 0, "RS02_precip", 0,
                  0, 0, 0, 0, 0, 0, 0, "SP02_precip",
                  "AB00_Tea_fire", 0, 0, 0, 0, 0, 0, 0,
                  "AB00_Jesusita_fire", 0, 0, 0, 0, 0, 0, 0,
                  0, "AT07_Jesusita_fire", 0, 0, 0, 0, 0, 0,
                  0, 0, "GV01_Gaviota_fire", 0, 0, 0, 0, 0,
                  0, 0, 0, "HO00_Gaviota_fire", 0, 0, 0, 0,
                  0, 0, 0, "HO00_Sherpa_fire", 0, 0, 0, 0,
                  0, 0, 0, 0, "MC06_Tea_fire", 0, 0, 0,
                  0, 0, 0, 0, "MC06_Jesusita_fire", 0, 0, 0,
                  0, 0, 0, 0, 0, "RG01_Gaviota_fire", 0, 0,
                  0, 0, 0, 0, 0, "RG01_Sherpa_fire", 0, 0,
                  0, 0, 0, 0, 0, 0, "RS01_Tea_fire", 0,
                  0, 0, 0, 0, 0, 0, "RS01_Jesusita_fire", 0,
                  0, 0, 0, 0, 0, 0, 0, "SP02_Gap_fire"),8,23)


# Model setup
mod_list <- list(
  # tinitx = "zero", # setting initial state value to time = 0
  B = "diagonal and unequal",
  U = "zero", # zero: does NOT allow a drift term in process model to be estimated # removing due to lack of anticipated monotonic trend
  C = CC, # see Alex's matrix above
  c = dat_cov, # we should probably de-mean and scale covariates to units of sd 
  Q = "diagonal and unequal", # diagonal and unequal: allows for and estimates the covariance matrix of process errors
  Z = "identity", # identity: estimated state processes are unique to each site
  A = "zero",
  R = "zero" # diagonal and equal: allows for and estimates the covariance matrix of observations errors (may want to provide a number for this from method precision etc if possible) - changed to "zero" on 11/22 to "turn off" observation error
)

# Fit model
# fit <- MARSS(y = dat_dep, model = mod_list,
#                 control = list(maxit = 5000), method = "BFGS")
# see pg 5 in MARSS manual for notes on method BFGS vs method EM: EM algorithm gives more robust estimation for datasets replete with missing values and for high-dimensional models with various constraints. BFGS is faster and is good enough for some datasets. Typically, both should be tried.

fit <- MARSS(y = dat_dep, model = mod_list,
             control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE) #default method = "EM"

# export model fit
saveRDS(fit, file = "data_working/marss_test_run/fit_112221_R0.rds")

#### Plotting Results ####

## check for hidden errors
fit[["errors"]]
# no hidden errors

### Plot coef and coef estimates ###

## estimates

# non-parametric:
#est = MARSSparamCIs(fit, method = "innovations", alpha = 0.05, nboot = 70, silent=F) 
# ^ Errors were caught in MARSSboot :Innovations bootstrapping uses the innovations resampling and can only be done if there are no missing values in the data.
# Error: Stopped in MARSSboot() due to problem(s) with function arguments.

# parametric:
est_fit = MARSSparamCIs(fit, method = "parametric", alpha = 0.05, nboot = 10, silent=F) # nboot should be ~ 2000 for final results
# note - this code takes a while to run, so can import "CIs_fit_112221_R0.csv"
# from the data_working/marss_test_run folder to bypass this and the 
# next few steps

# hessian method is much fast but not ideal for final results
# est_fit = MARSSparamCIs(fit)

CIs_fit = cbind(
  est_fit$par$U,
  est_fit$par.lowCI$U,
  est_fit$par.upCI$U)
CIs_fit = as.data.frame(CIs_fit)
names(CIs_fit) = c("Est.", "Lower", "Upper")
CIs_fit$parm = rownames(CIs_fit)
CIs_fit[,1:3] = round(CIs_fit[,1:3], 3)

## save CI table ##
write_csv(CIs_fit, "data_working/marss_test_run/CIs_fit_112221_R0.csv")

# Creating quick panel plots of covariate data to try and diagnose
# potential convergence issues:
(nh4plot <- ggplot(dat, aes(x = date, y = mean_nh4_uM)) +
  geom_point(aes(color = site)) +
  labs(x = "Date",
       y = "NH4+ (uM)") +
  facet_wrap(.~site) +
  theme_bw() +
  theme(legend.position = "none"))

(no3plot <- ggplot(dat, aes(x = date, y = mean_no3_uM)) +
    geom_point(aes(color = site)) +
    labs(x = "Date",
         y = "NO3- (uM)") +
    facet_wrap(.~site) +
    theme_bw() +
    theme(legend.position = "none"))

(tssplot <- ggplot(dat, aes(x = date, y = mean_tss_mgL)) +
    geom_point(aes(color = site)) +
    labs(x = "Date",
         y = "TSS (mg/L)") +
    facet_wrap(.~site) +
    theme_bw() +
    theme(legend.position = "none"))

# Stopped here on 11/22/21, H.L.

### PLOT HO00 ###
CIs_HO00 = rbind(CIs_fit[1:2,], CIs_fit[grepl("HO00", CIs_fit$parm),])
CIs_HO00$parm_name = c("Season1","Season2","Precip. (monthly cumm.)", "Gaviota Fire", "Sherpa Fire")
# plot
RESULTS_HO00 = 
  ggplot(CIs_HO00, aes(parm_name, Est.)) + 
  geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
  geom_point(position=position_dodge(width=0.3), size=2) + 
  theme_bw()+
  theme(plot.title = element_text(size = 8)) +
  theme(axis.text = element_text(size = 8)) +
  geom_hline(aes(yintercept=0), linetype="dashed")+
  coord_flip()+ ggtitle("HO00 NH4+")+ 
  ylab("") +
  theme(plot.margin=unit(c(.2,-.2,.05,.01),"cm"))
RESULTS_HO00

### PLOT RG01 ###
CIs_RG01 = rbind(CIs_fit[1:2,], CIs_fit[grepl("RG01", CIs_fit$parm),])
CIs_RG01$parm_name = c("Season1","Season2","Precip. (monthly cumm.)", "Gaviota Fire", "Sherpa Fire")
# plot
RESULTS_RG01 = 
  ggplot(CIs_RG01, aes(parm_name, Est.)) + 
  geom_errorbar(aes(ymin=Lower, ymax=Upper),position=position_dodge(width=0.25), width=0.25) +
  geom_point(position=position_dodge(width=0.3), size=2) + 
  theme_bw()+
  theme(plot.title = element_text(size = 8)) +
  theme(axis.text = element_text(size = 8)) +
  geom_hline(aes(yintercept=0), linetype="dashed")+
  coord_flip()+ ggtitle("RG01 NH4+")+ 
  ylab("") +
  theme(plot.margin=unit(c(.2,-.2,.05,.01),"cm"))
RESULTS_RG01

gridExtra::grid.arrange(RESULTS_HO00, RESULTS_RG01, nrow=1)
g <- gridExtra::arrangeGrob(RESULTS_HO00, RESULTS_RG01, nrow=1)
ggsave("figures/MARSS_NH4_HO00_RG01.pdf", g)


#### Script for diagnoses ####

dat = dat_dep
time = c(1:ncol(dat_dep))
resids <- residuals(fit)
kf=print(fit, what="kfs") # Kalman filter and smoother output

### Compare to null model ###
mod_list_null <- list(
  B = "diagonal and unequal",
  U = "zero", 
  Q = "diagonal and unequal", 
  Z = "identity", 
  A = "zero",
  R = "diagonal and equal" 
)
mod.null <- MARSS(y = dat_dep, model = mod_list_null,
                  control = list(maxit= 2000, allow.degen=TRUE, trace=1), fit=TRUE)

bbmle::AICtab(fit, mod.null)

#           dAIC df
# fit       0   15
# mod.null 11   7 
# RESULT: covariate model is better than null model - suggests covars add information

### Do resids have temporal autocorrelation? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
forecast::Acf(resids$model.residuals[1,], main="HO00 model residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$state.residuals[1,], main="HO00 state residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$model.residuals[2,], main="RG01 model residuals", na.action=na.pass, lag.max = 24)
forecast::Acf(resids$state.residuals[2,], main="RG01 state residuals", na.action=na.pass, lag.max = 24)
mtext("Do resids have temporal autocorrelation?", outer = TRUE, cex = 1.5)
# RESULT: yes, at lag 1 for RG01 and in multiple locations for HO00

### Do resids have temporal trend? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
plot(resids$model.residuals[1,], ylab="model residual", xlab="", main="HO00 model residuals")
abline(h=0)
plot(resids$state.residuals[1,], ylab="state residual", xlab="", main="HO00 state residuals")
abline(h=0)
plot(resids$model.residuals[2,], ylab="model residual", xlab="", main="RC01 model residuals")
abline(h=0)
plot(resids$state.residuals[2,], ylab="state residual", xlab="", main="RC01 state residuals")
abline(h=0)
mtext("Do resids have temporal trend?", outer = TRUE, cex = 1.5)

### Are resids normal? ###
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
qqnorm(resids$model.residuals[1,], main="HO00 model residuals", pch=16, 
       xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[1,])[1]))
qqline(resids$model.residuals[1,])
qqnorm(resids$state.residuals[1,], main="HO00 state residuals", pch=16, 
       xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[1,])[1]))
qqline(resids$state.residuals[1,])
qqnorm(resids$model.residuals[2,], main="RG01 model residuals", pch=16, 
       xlab=paste("shapiro test: ", shapiro.test(resids$model.residuals[2,])[1]))
qqline(resids$model.residuals[2,])
qqnorm(resids$state.residuals[2,], main="RG01 state residuals", pch=16, 
       xlab=paste("shapiro test: ", shapiro.test(resids$state.residuals[2,])[1]))
qqline(resids$state.residuals[2,])
mtext("Are resids normal?", outer = TRUE, cex = 1.5)

### residuals vs fitted ###
par(mfrow=c(2,2),oma = c(2, 0, 2, 0))
scatter.smooth(as.vector(kf$xtT[1,]) ~ resids$model.residuals[1,],
               main = "HO00 model resids vs fitted (y conditioned)")
abline(lm(as.vector(kf$xtT[1,]) ~ resids$model.residuals[1,]), col="blue")
scatter.smooth(as.vector(kf$xtT[1,]) ~ resids$state.residuals[1,],
               main = "HO00 state resids vs fitted (y conditioned)")
abline(lm(as.vector(kf$xtT[1,]) ~ resids$state.residuals[1,]), col="blue")
scatter.smooth(as.vector(kf$xtT[2,]) ~ resids$model.residuals[2,],
               main = "RC01 model resids vs fitted (y conditioned)")
abline(lm(as.vector(kf$xtT[2,]) ~ resids$model.residuals[2,]), col="blue")
scatter.smooth(as.vector(kf$xtT[2,]) ~ resids$state.residuals[2,],
               main = "RC01 state resids vs fitted (y conditioned)")
abline(lm(as.vector(kf$xtT[2,]) ~ resids$state.residuals[2,]), col="blue")
mtext("Trends in fitted values vs residuals?", outer = TRUE, cex = 1.5)

### residuals vs observed ###
par(mfrow=c(2,2),oma = c(2, 0, 2, 0))
scatter.smooth(as.vector(dat[1,]) ~ resids$model.residuals[1,],
               main = "HO00 model resids vs observed")
abline(lm(as.vector(kf$xtT[1,]) ~ resids$model.residuals[1,]))
scatter.smooth(as.vector(dat[1,]) ~ resids$state.residuals[1,], 
               main = "HO00 state resids vs observed") 
abline(lm(as.vector(kf$xtT[1,]) ~ resids$state.residuals[1,]))
scatter.smooth(as.vector(dat[2,]) ~ resids$model.residuals[2,],
               main = "RC01 model resids vs observed")
abline(lm(as.vector(kf$xtT[2,]) ~ resids$model.residuals[2,]))
scatter.smooth(as.vector(dat[2,]) ~ resids$state.residuals[2,], 
               main = "RC01 state resids vs observed") 
abline(lm(as.vector(kf$xtT[2,]) ~ resids$state.residuals[2,]))
mtext("Trends in observed values vs residuals?", outer = TRUE, cex = 1.5)

# reset plotting window
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

# End of script.

