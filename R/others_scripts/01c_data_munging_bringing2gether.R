# updated: JMH 14 Sept. 21
##-------------------
## required packages
##-------------------
library(here)
library(readr)
library(tidyverse)
library(readxl)
library(GGally)
library(openxlsx)

# Common time period
CTstart <- as.POSIXct(paste0("1985-10-31", format = "%Y-%m-%d"))
# changed this to include more data 26 Jul 21
CTend <- as.POSIXct(paste0("2021-10-31", format = "%Y-%m-%d"))
# 13,149 - updated 24 sept 21


##-------------------
# Get data
##-------------------
# Detection limit data, prepared by David and Tamara 
# from info provided by data providors
LOQ1 <- readxl::read_xlsx(file.path(here::here("data/JMHnewMungedDat"), 
                                   "maximum detection limits.xlsx"),1, 
                         range = "A1:H12") 

names(LOQ1) <- c("Site", "SiteCode", "DL_Ca_mgL", "DL_DOC_mgL", "DL_NH4_N_mgL", "DL_NO3_N_mgL", "DL_TP_TDP_PO4_P_mgL", "DL_SO4_S_mgL")

LOQ <- LOQ1 %>% 
  mutate(across(DL_NH4_N_mgL:DL_TP_TDP_PO4_P_mgL, as.numeric)) %>% 
  # calculate value to replace with if < DL
  # RV = replacement value, which is 1/2 of the detection limit
  # DL = "Detection limit"
  mutate(RV_Ca_mgL = DL_Ca_mgL/2,
         RV_DOC_mgL = DL_DOC_mgL/2,
         RV_NH4_N_mgL = DL_NH4_N_mgL/2,
         RV_NO3_N_mgL = DL_NO3_N_mgL/2,
         RV_TP_TDP_PO4_P_mgL = DL_TP_TDP_PO4_P_mgL/2,
         RV_SO4_S_mgL = DL_SO4_S_mgL/2) %>% 
  select(-Site) %>% 
  #drop sites not used
  filter(!(SiteCode %in% c("COW", "LUQ", "SAN")))
  

BBWM <- read.csv(file.path(here::here("data/JMHnewMungedDat"), 
                  "01_BBWMcomb.csv"))

DOR <- read.csv(file.path(here::here("data/JMHnewMungedDat"), 
                           "01_DORcomb.csv"))

ELA <- read.csv(file.path(here::here("data/JMHnewMungedDat"), 
                           "01_ELAcomb.csv"))

HBEF <- read.csv(file.path(here::here("data/JMHnewMungedDat"), 
                           "01_HBEFcomb.csv"))

HJA <- read.csv(file.path(here::here("data/JMHnewMungedDat"), 
                           "01_HJAcomb.csv"))

MEF <- read.csv(file.path(here::here("data/JMHnewMungedDat"), 
                           "01_MEFcomb.csv"))

SEF <- read.csv(file.path(here::here("data/JMHnewMungedDat"), 
                           "01_SEFcomb.csv"))

SLP <- read.csv(file.path(here::here("data/JMHnewMungedDat"), 
                           "01_SLPcomb.csv"))

TLW <- read.csv(file.path(here::here("data/JMHnewMungedDat"), 
                           "01_TLWcomb.csv"))

# combine solute data
sol2 <- rbind(BBWM, DOR, ELA, HBEF, HJA, MEF, SEF, SLP, TLW) %>% 
          mutate_at(vars(Site, WS), funs(factor)) %>% 
          mutate(Date = as.Date(Date, format = "%Y-%m-%d"),
                 SiteWs = as.factor(paste0(Site,"_",WS))) %>% 
          # average to day to get rid of a couple replicate samples
          group_by(SiteWs, Date) %>% 
          summarise(across(Q_Ls:SO4_mgL, ~ mean(.x, na.rm = TRUE))) %>% 
          # truncate to focal dates, but not really triming at end
          filter(Date >= CTstart & Date <= CTend) %>% 
          separate(SiteWs, sep =  "_", into = c("Site", "WS"), remove = FALSE) %>% 
          mutate(Site = as.factor(Site),
                 WS = as.factor(WS)) %>% 
          # need to drop SEF
          filter(Site != "SEF") %>% 
          # What TDP really is: BBWM = N/A; DOR = TP; ELA = TDP; HBEF = SRP; HJA = TDP; MEF = TP (DA: little PP so TP = TDP)
          # SLP = N/A (no SRP/TDP data in window); TLW = TP (DA: little PP so TP = TDP)
          rename(TDP_mgL = "SRP_mgL")

sol1 <- sol2 %>% 
        # join with detection limit info
        left_join(LOQ, by = c("Site" = "SiteCode")) %>% 
        # if below DL then use DL/2
        mutate(Ca_mgL.u = ifelse(Ca_mgL < DL_Ca_mgL, RV_Ca_mgL, Ca_mgL),
               DOC_mgL.u = ifelse(DOC_mgL < DL_DOC_mgL, RV_DOC_mgL, DOC_mgL),
               NH4_mgL.u = ifelse(NH4_mgL < DL_NH4_N_mgL, RV_NH4_N_mgL, NH4_mgL),
               NO3_mgL.u = ifelse(NO3_mgL < DL_NO3_N_mgL, RV_NO3_N_mgL, NO3_mgL),
               TDP_mgL.u = ifelse(TDP_mgL < DL_TP_TDP_PO4_P_mgL, RV_TP_TDP_PO4_P_mgL, TDP_mgL),
               SO4_mgL.u = ifelse(SO4_mgL < DL_SO4_S_mgL, RV_SO4_S_mgL, SO4_mgL)) 

# %>%
#         # drop DL and RV's
#         select(-c(DL_Ca_mgL:RV_SO4_S_mgL))

# NEED TO RUN THESE CHECKS
  # Calcium
  summary(sol1$DL_Ca_mgL)
  ggplot() +
    geom_point(data = sol1, aes(y = Ca_mgL, x = Date)) +
    geom_hline(data = sol1, aes(yintercept = DL_Ca_mgL), color = "red") +
    facet_wrap(vars(SiteWs), scales = "free_y") 
  
  # DOC
  summary(sol1$DL_DOC_mgL)
  ggplot() +
    geom_point(data = sol1, aes(y = DOC_mgL, x = Date)) +
    geom_hline(data = sol1, aes(yintercept = DL_DOC_mgL), color = "red") +
    facet_wrap(vars(SiteWs), scales = "free_y") 
  
  # NH4 - matters a lot!
  summary(sol1$DL_NH4_N_mgL)
  ggplot() +
    geom_point(data = sol1, aes(y = NH4_mgL, x = Date)) +
    geom_hline(data = sol1, aes(yintercept = DL_NH4_N_mgL), color = "red") +
    facet_wrap(vars(SiteWs), scales = "free_y") 
  
  ggplot() +
    geom_point(data = sol1, aes(y = NH4_mgL.u, x = Date)) +
    geom_hline(data = sol1, aes(yintercept = DL_NH4_N_mgL), color = "red") +
    scale_y_log10()+
    facet_wrap(vars(SiteWs), scales = "free_y") 
  
  # NO3 - also matters a lot
  summary(sol1$DL_NO3_N_mgL)
  ggplot() +
    geom_point(data = sol1, aes(y = NO3_mgL, x = Date)) +
    geom_hline(data = sol1, aes(yintercept = DL_NO3_N_mgL), color = "red") +
    facet_wrap(vars(SiteWs), scales = "free_y") 
  
  ggplot() +
    geom_point(data = sol1, aes(y = NO3_mgL.u, x = Date)) +
    geom_hline(data = sol1, aes(yintercept = DL_NO3_N_mgL), color = "red") +
    scale_y_log10()+
    facet_wrap(vars(SiteWs), scales = "free_y") 
  
  # TDP - also matters a lot for HBEF and TLW
  summary(sol1$DL_TP_TDP_PO4_P_mgL)
  ggplot() +
    geom_point(data = sol1, aes(y = TDP_mgL, x = Date)) +
    geom_hline(data = sol1, aes(yintercept = DL_TP_TDP_PO4_P_mgL), color = "red") +
    facet_wrap(vars(SiteWs), scales = "free_y") 
  
  ggplot() +
    geom_point(data = sol1, aes(y = TDP_mgL.u, x = Date)) +
    geom_hline(data = sol1, aes(yintercept = DL_TP_TDP_PO4_P_mgL), color = "red") +
    scale_y_log10()+
    facet_wrap(vars(SiteWs), scales = "free_y") 
  
  # SO4 - also matters a lot
  summary(sol1$DL_SO4_S_mgL)
  ggplot() +
    geom_point(data = sol1, aes(y = SO4_mgL.u, x = Date)) +
    geom_hline(data = sol1, aes(yintercept = DL_SO4_S_mgL), color = "red") +
    facet_wrap(vars(SiteWs), scales = "free_y") 
    
# still duplicates?
  dim(sol1[duplicated(sol1[,c("SiteWs", "Date")]) == TRUE,])
  
# remove raw data from dataframe
  sol <-  sol1 %>% 
    select(SiteWs:Q_Ls, Ca_mgL.u:SO4_mgL.u) %>% 
    # remove u - so these values are now adjusted for the site-specific DL
    rename(Ca_mgL = Ca_mgL.u, DOC_mgL = DOC_mgL.u, NH4_mgL = NH4_mgL.u, NO3_mgL = NO3_mgL.u, 
           TDP_mgL = TDP_mgL.u, SO4_mgL = SO4_mgL.u)

  
# Make conc plots
# all sites
  # NOTE: NOT ALL TDP VALUES ARE TDP SEE WATERSHED DATA NOTES SPREADSHEET
pdf(file.path(here::here("plots"),
             "RawConcPlotAllSites20210915.pdf"), width = 25, height = 10)
ggplot(sol %>%
         pivot_longer(col = c(Q_Ls:SO4_mgL), names_to = "solute", values_to = "conc") %>%
         mutate(solute = fct_relevel(solute,
                                     c("Q_Ls", "Ca_mgL", "DOC_mgL", "NH4_mgL", "NO3_mgL", "TDP_mgL", "SO4_mgL"))), aes(y = conc, x = Date, color = Site)) +
          geom_point(size = 0.5) +
        facet_grid(solute ~ SiteWs, scales = "free_y")

dev.off()

# TLW
pdf(file.path(here::here("plots"),
              "RawConcPlotTLW_20210915.pdf"), width = 15, height = 10)
ggplot(sol %>%
         # rename(TDP_mgL = "SRP_mgL") %>% 
         pivot_longer(col = c(Q_Ls:SO4_mgL), names_to = "solute", values_to = "conc") %>%
         mutate(solute = fct_relevel(solute,
                                     c("Q_Ls", "Ca_mgL", "DOC_mgL", "NH4_mgL", "NO3_mgL", "TDP_mgL", "SO4_mgL"))) %>% 
         filter(Site == "TLW"), aes(y = conc, x = Date, color = Site)) +
  geom_point(size = 0.5) +
  facet_grid(solute ~ SiteWs, scales = "free_y") 
dev.off()

# SLP ------- DROPPED
# SRP here is very limited, not sure if SRP or TDP

# pdf(file.path(here::here("plots"),
#               "RawConcPlotSLP_20210725.pdf"), width = 5, height = 10)
# ggplot(sol %>%
#          pivot_longer(col = c(Q_Ls:SO4_mgL), names_to = "solute", values_to = "conc") %>%
#          mutate(solute = fct_relevel(solute,
#                                      c("Q_Ls", "Ca_mgL", "DOC_mgL", "NH4_mgL", "NO3_mgL", "TDP_mgL", "SO4_mgL"))) %>% 
#          filter(Site == "SLP"), aes(y = conc, x = Date, color = Site)) +
#   geom_point(size = 0.5) +
#   facet_grid(solute ~ SiteWs, scales = "free_y")
# dev.off()

# MEF
pdf(file.path(here::here("plots"),
              "RawConcPlotMEF_20210915.pdf"), width = 10, height = 10)
ggplot(sol %>%
         # rename(TDP_mgL = "SRP_mgL") %>% 
         pivot_longer(col = c(Q_Ls:SO4_mgL), names_to = "solute", values_to = "conc") %>%
         mutate(solute = fct_relevel(solute,
                                     c("Q_Ls", "Ca_mgL", "DOC_mgL", "NH4_mgL", "NO3_mgL", "TDP_mgL", "SO4_mgL"))) %>% 
         filter(Site == "MEF"), aes(y = conc, x = Date, color = Site)) +
  geom_point(size = 0.5) +
  facet_grid(solute ~ SiteWs, scales = "free_y")
dev.off()

# HJA
pdf(file.path(here::here("plots"),
              "RawConcPlotHJA_20210915.pdf"), width = 10, height = 10)
ggplot(sol %>%
         # rename(TDP_mgL = "SRP_mgL") %>% 
         pivot_longer(col = c(Q_Ls:SO4_mgL), names_to = "solute", values_to = "conc") %>%
         mutate(solute = fct_relevel(solute,
                                     c("Q_Ls", "Ca_mgL", "DOC_mgL", "NH4_mgL", "NO3_mgL", "TDP_mgL", "SO4_mgL"))) %>% 
         filter(Site == "HJA"), aes(y = conc, x = Date, color = Site)) +
  geom_point(size = 0.5) +
  facet_grid(solute ~ SiteWs, scales = "free_y")
dev.off()

# HBEF
pdf(file.path(here::here("plots"),
              "RawConcPlotHBEF_20210915.pdf"), width = 20, height = 10)
ggplot(sol %>%
         pivot_longer(col = c(Q_Ls:SO4_mgL), names_to = "solute", values_to = "conc") %>%
         mutate(solute = fct_relevel(solute,
                                     c("Q_Ls", "Ca_mgL", "DOC_mgL", "NH4_mgL", "NO3_mgL", "TDP_mgL", "SO4_mgL"))) %>% 
         filter(Site == "HBEF"), aes(y = conc, x = Date, color = Site)) +
  geom_point(size = 0.5) +
  facet_grid(solute ~ SiteWs, scales = "free_y")
dev.off()

# ELA
pdf(file.path(here::here("plots"),
              "RawConcPlotELA_20210915.pdf"), width = 12, height = 10)
ggplot(sol %>%
         # rename(TDP_mgL = "SRP_mgL") %>% 
         pivot_longer(col = c(Q_Ls:SO4_mgL), names_to = "solute", values_to = "conc") %>%
         mutate(solute = fct_relevel(solute,
                                     c("Q_Ls", "Ca_mgL", "DOC_mgL", "NH4_mgL", "NO3_mgL", "TDP_mgL", "SO4_mgL"))) %>% 
         filter(Site == "ELA"), aes(y = conc, x = Date, color = Site)) +
  geom_point(size = 0.5) +
  facet_grid(solute ~ SiteWs, scales = "free_y")
dev.off()

# DOR
pdf(file.path(here::here("plots"),
              "RawConcPlotDOR_20210915.pdf"), width = 25, height = 10)
ggplot(sol %>%
         # rename(TP_mgL = "SRP_mgL") %>% 
         pivot_longer(col = c(Q_Ls:SO4_mgL), names_to = "solute", values_to = "conc") %>%
         mutate(solute = fct_relevel(solute,
                                     c("Q_Ls", "Ca_mgL", "DOC_mgL", "NH4_mgL", "NO3_mgL", "TDP_mgL", "SO4_mgL"))) %>% 
         filter(Site == "DOR"), aes(y = conc, x = Date, color = Site)) +
  geom_point(size = 0.5) +
  facet_grid(solute ~ SiteWs, scales = "free_y")
dev.off()

# BBWM
pdf(file.path(here::here("plots"),
              "RawConcPlotBBWM_20210915.pdf"), width = 5, height = 10)
ggplot(sol %>%
         pivot_longer(col = c(Q_Ls:SO4_mgL), names_to = "solute", values_to = "conc") %>%
         mutate(solute = fct_relevel(solute,
                                     c("Q_Ls", "Ca_mgL", "DOC_mgL", "NH4_mgL", "NO3_mgL", "TDP_mgL", "SO4_mgL"))) %>% 
         filter(Site == "BBWM"), aes(y = conc, x = Date, color = Site)) +
  geom_point(size = 0.5) +
  facet_grid(solute ~ SiteWs, scales = "free_y")
dev.off()


##-------------------
# Calculate FWMC
##-------------------



# using eq here: https://ncwqr.files.wordpress.com/2017/06/d-time-weighted-and-flow-weighted-mean-concentrations.pdf
# Campbell w/ hubbard brook infilled the conc values with means
# going to have to get sol in the right form
# HJA IS FWMC ALREADY
solFW <- sol %>% 
  ungroup() %>% 
  # add Interval col
  mutate(Interval = as.numeric("NA")) %>% 
  # wide to long
  pivot_longer(cols = c(Ca_mgL:SO4_mgL), names_to = "solute", values_to = "mgL") %>% 
  mutate(SiteWsSol = paste0(SiteWs, "_", solute)) %>% 
  select(SiteWsSol, Date, Q_Ls, mgL, Interval) %>% 
  # need to remove NA's
  filter(!is.na(Q_Ls)) %>% 
  filter(!is.na(mgL)) %>% 
  # need to arrange
  arrange(SiteWsSol, Date)

# This calculates the interval between samples
SiteWsSolnames <-  unique(solFW$SiteWsSol)

# duplicate for loop
solFW2 <- solFW

for(w in 1:length(SiteWsSolnames)){
  # w <- 1
  SiteWsSolnames_w <- SiteWsSolnames[w]
  sol2_w <- solFW[solFW$SiteWsSol == SiteWsSolnames_w,]
  sol2_wDates <- sol2_w$Date
  
  for(t in 2:length(sol2_wDates)){
    # t=2
    t_i <- sol2_wDates[t]
    t_im1 <- sol2_wDates[t-1]
    Interval_ti <- as.numeric(t_i - t_im1)
    sol2_w[t,]$Interval <- Interval_ti # this will be in days
  }
  # put back in new df
  solFW2[solFW2$SiteWsSol == SiteWsSolnames_w,]$Interval <-  sol2_w$Interval
}

# calculate components of FWMC
solFW3 <- solFW2 %>% 
  mutate(Interval = as.numeric(Interval)) %>% 
  separate(SiteWsSol, sep = "_", into = c("Site", "WS", "Solute", "conc")) %>% 
  # don't really need conc indicator since always mgL
  select(-conc) %>% 
  # make cols a factor
  mutate(across(c(Site, WS, Solute), factor)) %>% 
  mutate(SiteWs = as.factor(paste0(Site,"_", WS))) %>% 
  # calculate FWMC
  # There's a wide range in this interval 1-5649 days
  # <60 = 117653; <30 = 116448; <14 = 101975 <- numbers from the last round of this, but general point stands
  # dim(sol3[sol3$Interval <14,])
  # I'd rather do < 10
  mutate(FWMCtop = ifelse(Interval <= 30, 
                          mgL * Interval * (Q_Ls*60*60*24),
                          as.numeric("NA")), # convert Ls to L/d, as interval is in day
         FWMCbottom = ifelse(Interval <= 30,
                             Interval * (Q_Ls*60*60*24),
                             as.numeric("NA")),
  # this is only for daily stuff - summing for monthly
         FWMC = FWMCtop/FWMCbottom) %>% 
  # HJA is already FWMC - this also circumvents that calculation
  mutate(FWMC = ifelse(Site == "HJA", mgL, FWMC))


pdf(file.path(here::here("plots"),
              "IntervalsBWconcSmp_20210915.pdf"), width = 10, height = 25)
ggplot(solFW3, aes(y = log10(Interval), x = Date)) +
  geom_point()+
  facet_grid(SiteWs ~ Solute, scales = "free_y") +
  geom_hline(yintercept = log10(14), color = "red") +
  geom_hline(yintercept = log10(30), color = "green") 
dev.off()

ggplot(solFW3, aes(y = FWMC, x = Date)) +
  geom_point()+
  facet_grid(Solute ~ SiteWs, scales = "free_y")

hist(log10(solFW3$Interval))


# summarize to monthly
# avg conc and and Q data
solM_Qconc <- solFW3 %>% 
  # make Y and Month cols
  mutate(Y = strftime(Date, format = "%Y"),
         M = strftime(Date, format = "%m")) %>% 
  # average to month
  group_by(SiteWs, Y, M, Solute) %>% 
  summarise(across(c(Q_Ls, mgL), ~mean(.x, na.rm = TRUE)))

# sum FWMC to monthly
solM_FWMC <- solFW3 %>% 
  # make Y and Month cols
  mutate(Y = strftime(Date, format = "%Y"),
         M = strftime(Date, format = "%m")) %>% 
  # sum to month
  group_by(SiteWs, Y, M, Solute) %>% 
  summarise(across(c(FWMCtop, FWMCbottom), ~sum(.x, na.rm = TRUE))) %>% 
  mutate(FWMC = FWMCtop/FWMCbottom)
 
# combine conc, Q, and FWMC
solM <- solM_Qconc %>% 
  full_join(solM_FWMC, by = c("SiteWs", "Y", "M", "Solute")) %>% 
  mutate(Date = as.Date(paste0(Y,"-",M, "-01"), format = "%Y-%m-%d")) %>% 
  separate(SiteWs, sep = "_", into = c("Site", "WS"), remove = FALSE) %>% 
  # HJA conc = FWMC, so we have to do this again
  mutate(FWMC = ifelse(Site == "HJA", mgL, FWMC)) %>% 
  ungroup() %>% 
  select(-FWMCtop, -FWMCbottom,-Y,-M) %>% 
  mutate(Site = as.factor(Site),
         WS = as.factor(WS)) %>% 
  # remove SEF
  filter(Site != "SEF") %>% 
  # removing site/solute that are all (or close to all) NAs
  # very limited SRP data in SLP - remove
  filter(!(Site == "SLP" & Solute == "SRP")) %>% 
  # MEF NO3 and NH4 data unreliable
  filter(!(Site == "MEF" & Solute == "NO3")) %>% 
  filter(!(Site == "MEF" & Solute == "NH4")) %>% 
  #no SRP data
  filter(!(Site == "BBWM" & Solute == "SRP")) 
  
  
  
  
   
ggplot(solM, aes(y = FWMC, x = mgL, color = Solute)) +
  geom_point() +
  facet_wrap(vars(Solute), scales = "free") +
  geom_abline(intercept = 0, slope = 1)

# export FWMC plots
# all site
pdf(file.path(here::here("plots"),
             "FWMCPlotAllSites_20210915.pdf"), width = 25, height = 10)
ggplot(solM, aes(y = FWMC, x = Date, color = Site)) +
  geom_point()+
  facet_grid(Solute ~ SiteWs, scales = "free_y") +
  ylab("Flow-weighted mean conc (mg solute/L)")
dev.off()

#BBWM
pdf(file.path(here::here("plots"),
              "FWMCPlotBBWM_20210915.pdf"), width = 5, height = 10)
ggplot(solM %>% 
         filter(Site == "BBWM"), aes(y = FWMC, x = Date)) +
  geom_point()+
  facet_grid(Solute ~ SiteWs, scales = "free_y") +
  ylab("Flow-weighted mean conc (mg solute/L)")
dev.off()

#DOR
pdf(file.path(here::here("plots"),
              "FWMCPlotDOR_20210915.pdf"), width = 25, height = 10)
ggplot(solM %>% 
         filter(Site == "DOR"), aes(y = FWMC, x = Date)) +
  geom_point()+
  facet_grid(Solute ~ SiteWs, scales = "free_y") +
  ylab("Flow-weighted mean conc (mg solute/L)")
dev.off()

#ELA
pdf(file.path(here::here("plots"),
              "FWMCPlotELA_20210915.pdf"), width = 15, height = 10)
ggplot(solM %>% 
         filter(Site == "ELA"), aes(y = FWMC, x = Date)) +
  geom_point()+
  facet_grid(Solute ~ SiteWs, scales = "free_y") +
  ylab("Flow-weighted mean conc (mg solute/L)")
dev.off()

#HBEF
pdf(file.path(here::here("plots"),
              "FWMCPlotHBEF_20210915.pdf"), width = 25, height = 10)
ggplot(solM %>% 
         filter(Site == "HBEF"), aes(y = FWMC, x = Date)) +
  geom_point()+
  facet_grid(Solute ~ SiteWs, scales = "free_y") +
  ylab("Flow-weighted mean conc (mg solute/L)")
dev.off()

#HJA
pdf(file.path(here::here("plots"),
              "FWMCPlotHJA_20210915.pdf"), width = 10, height = 10)
ggplot(solM %>% 
         filter(Site == "HJA"), aes(y = FWMC, x = Date)) +
  geom_point()+
  facet_grid(Solute ~ SiteWs, scales = "free_y") +
  ylab("Flow-weighted mean conc (mg solute/L)")
dev.off()

#MEF
pdf(file.path(here::here("plots"),
              "FWMCPlotMEF_20210915.pdf"), width = 10, height = 10)
ggplot(solM %>% 
         filter(Site == "MEF"), aes(y = FWMC, x = Date)) +
  geom_point()+
  facet_grid(Solute ~ SiteWs, scales = "free_y") +
  ylab("Flow-weighted mean conc (mg solute/L)")
dev.off()

#SEF - no data
# pdf(file.path(here::here("data/plots"),
#               "FWMCPlotSEF.pdf"), width = 25, height = 10)
# ggplot(solM %>% 
#          filter(Site == "SEF"), aes(y = FWMC, x = Date, color = Site)) +
#   geom_point()+
#   facet_grid(Solute ~ SiteWs, scales = "free_y") +
#   ylab("Flow-weighted mean conc (mg solute/L)")
# dev.off()

#SLP
pdf(file.path(here::here("plots"),
              "FWMCPlotSLP_20210915.pdf"), width = 5, height = 10)
ggplot(solM %>% 
         filter(Site == "SLP"), aes(y = FWMC, x = Date)) +
  geom_point()+
  facet_grid(Solute ~ SiteWs, scales = "free_y") +
  ylab("Flow-weighted mean conc (mg solute/L)")
dev.off()

#TLW
pdf(file.path(here::here("plots"),
              "FWMCPlotTLW_20210915.pdf"), width = 15, height = 10)
ggplot(solM %>% 
         filter(Site == "TLW"), aes(y = FWMC, x = Date)) +
  geom_point()+
  facet_grid(Solute ~ SiteWs, scales = "free_y") +
  ylab("Flow-weighted mean conc (mg solute/L)")
dev.off()

# notes on WS/solutes to drop
#BBWM: lots of forced values for NH4
# DOR: lots of forced values for SO4 - DROP SO4
# ELA looks fine
# HBEF - Lots of forced values for NH4 and NO3- might want to drop NH4 and TDP for WS7-9
    # HBEF why do there seem to be two lines in NH4? Seems like the values for WS7-9 have 2 sig units. Also they seem to have forced a lot to zero.
# HJA- drop TDP
# MEF - fine
# TLW - fine



# export dataframe
write.csv(solM, file.path(here::here("data/JMHnewMungedDat"), 
                            "02_MonthlyQConcFWMCallSites_20210915.csv"))


# STILL NEED TOMAKE DATAFRAME FOR MARS, but have to get Tamara's thoughts on what to drop.


###########
# make dataframe for MARs
##########
# create blank df with every month

TSxts <- as.character(seq(as.Date("1985-11-01"), length = 428, by = "months"))
# add in site_WS
SiteWs <- levels(solM$SiteWs)
TSxts2 <- rep(TSxts, times = length(SiteWs))
SiteWS2 <- rep(SiteWs, each = length(TSxts))
BlankTS <- as.data.frame(cbind(SiteWS2, TSxts2))
names(BlankTS) <-  c("SiteWs", "Date")



MARSdf <- solM %>% 
  ungroup() %>% 
  mutate(SiteWs = paste0(Site,"_",WS)) %>% 
  select(SiteWs, Date, Solute, FWMC) %>% 
  pivot_wider(names_from = "Solute", values_from = "FWMC")

MARSdf2 <- BlankTS %>% 
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>% 
  full_join(MARSdf, by = c("SiteWs", "Date")) %>% 
  separate(SiteWs, sep = "_", into= c("Site", "WS"), remove = FALSE) %>% 
  rename(Ca_fwmc_mgL = "Ca", DOC_fwmc_mgL = "DOC", NH4_fwmc_mgL = "NH4",
         NO3_fwmc_mgL = "NO3", SO4_fwmc_mgL = "SO4", 
         # TDP is target
         # exceptions: DOR = TP, HBEF = SRP
         TDP_fwmc_mgL = "TDP")

# How much data is missing for each WS?
MARSmissingDataByWS <- MARSdf2 %>% 
  group_by(Site, WS) %>% 
  summarise(across(Date:TDP_fwmc_mgL, ~round(sum(is.na(.))/301*100,0))) #there should be 301 rows for each WS

MARSmissingDataBySite <- MARSmissingDataByWS %>% 
  group_by(Site) %>% 
  summarise(across(Date:TDP_fwmc_mgL, mean))

# export MARS data
write.csv(MARSdf2, file.path(here::here("data/JMHnewMungedDat"), 
                   "02_Dat4Mars_2mostRecent.csv"))

# NA's by WS
write.csv(MARSmissingDataByWS, file.path(here::here("data/JMHnewMungedDat"), 
                             "MARSmissingDataByWS_2mostRecent.csv"))

# NA's by site
write.csv(MARSmissingDataBySite, file.path(here::here("data/JMHnewMungedDat"), 
                             "MARSmissingDataBySite_2mostRecent.csv"))


save.image(file.path(here::here("analysis"),
                     "01c_data_munging_bringing2getherRdat"))


# load(file.path(here::here("analysis"),
#                "01c_data_munging_bringing2getherRdat"))







