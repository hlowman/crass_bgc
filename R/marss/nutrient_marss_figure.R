# Nutrient MARSS models
# Compiled results figure
# Script started October 30, 2023 by Heili Lowman

#### Setup ####

# Load packages.
library(tidyverse)
library(lubridate)
library(MARSS)
library(naniar) 
library(here)
library(bbmle)
library(broom)
library(stats4)

# Load in data.

# NH4 #
# no legacy, 1 state
noleg_nh4 <- readRDS(file = "data_working/marss_fits/fit_092123_1state_nh4_mBFGS.rds")
# 1y legacy, 1 state
leg1_nh4 <- readRDS(file = "data_working/marss_fits/fit_092123_1state_nh4_1ylegacy_mBFGS.rds")
# 2y legacy, 1 state
leg2_nh4 <- readRDS(file = "data_working/marss_fits/fit_092123_1state_nh4_2ylegacy_mBFGS.rds")
# 3y legacy, 1 state
leg3_nh4 <- readRDS(file = "data_working/marss_fits/fit_092123_1state_nh4_3ylegacy_mBFGS.rds")
# 4y legacy, 1 state
leg4_nh4 <- readRDS(file = "data_working/marss_fits/fit_092123_1state_nh4_4ylegacy_mBFGS.rds")
# 5y legacy, 1 state
leg5_nh4 <- readRDS(file = "data_working/marss_fits/fit_092123_1state_nh4_5ylegacy_mBFGS.rds")

# NO3 #
# no legacy, 1 state
noleg_no3 <- readRDS(file = "data_working/marss_fits/fit_092123_1state_no3_mBFGS.rds")
# 1y legacy, 1 state
leg1_no3 <- readRDS(file = "data_working/marss_fits/fit_092123_1state_no3_1ylegacy_mBFGS.rds")
# 2y legacy, 1 state
leg2_no3 <- readRDS(file = "data_working/marss_fits/fit_092123_1state_no3_2ylegacy_mBFGS.rds")
# 3y legacy, 1 state
leg3_no3 <- readRDS(file = "data_working/marss_fits/fit_092123_1state_no3_3ylegacy_mBFGS.rds")
# 4y legacy, 1 state
leg4_no3 <- readRDS(file = "data_working/marss_fits/fit_092123_1state_no3_4ylegacy_mBFGS.rds")
# 5y legacy, 1 state
leg5_no3 <- readRDS(file = "data_working/marss_fits/fit_092123_1state_no3_5ylegacy_mBFGS.rds")

# PO4 #
# no legacy, 1 state
noleg_po4 <- readRDS(file = "data_working/marss_fits/fit_092123_1state_po4_mBFGS.rds")
# 1y legacy, 1 state
leg1_po4 <- readRDS(file = "data_working/marss_fits/fit_092123_1state_po4_1ylegacy_mBFGS.rds")
# 2y legacy, 1 state
leg2_po4 <- readRDS(file = "data_working/marss_fits/fit_092123_1state_po4_2ylegacy_mBFGS.rds")
# 3y legacy, 1 state
leg3_po4 <- readRDS(file = "data_working/marss_fits/fit_092123_1state_po4_3ylegacy_mBFGS.rds")
# 4y legacy, 1 state
leg4_po4 <- readRDS(file = "data_working/marss_fits/fit_092123_1state_po4_4ylegacy_mBFGS.rds")
# 5y legacy, 1 state
leg5_po4 <- readRDS(file = "data_working/marss_fits/fit_092123_1state_po4_5ylegacy_mBFGS.rds")

#### Extract CIs ####

# Extract necessary confidence interval info
noleg_nh499 <- MARSSparamCIs(noleg_nh4, alpha = 0.01)
leg1y_nh499 <- MARSSparamCIs(leg1_nh4, alpha = 0.01)
leg2y_nh499 <- MARSSparamCIs(leg2_nh4, alpha = 0.01)
leg3y_nh499 <- MARSSparamCIs(leg3_nh4, alpha = 0.01)
leg4y_nh499 <- MARSSparamCIs(leg4_nh4, alpha = 0.01)
leg5y_nh499 <- MARSSparamCIs(leg5_nh4, alpha = 0.01)
noleg_nh495 <- MARSSparamCIs(noleg_nh4, alpha = 0.05)
leg1y_nh495 <- MARSSparamCIs(leg1_nh4, alpha = 0.05)
leg2y_nh495 <- MARSSparamCIs(leg2_nh4, alpha = 0.05)
leg3y_nh495 <- MARSSparamCIs(leg3_nh4, alpha = 0.05)
leg4y_nh495 <- MARSSparamCIs(leg4_nh4, alpha = 0.05)
leg5y_nh495 <- MARSSparamCIs(leg5_nh4, alpha = 0.05)

noleg_no399 <- MARSSparamCIs(noleg_no3, alpha = 0.01)
leg1y_no399 <- MARSSparamCIs(leg1_no3, alpha = 0.01)
leg2y_no399 <- MARSSparamCIs(leg2_no3, alpha = 0.01)
leg3y_no399 <- MARSSparamCIs(leg3_no3, alpha = 0.01)
leg4y_no399 <- MARSSparamCIs(leg4_no3, alpha = 0.01)
leg5y_no399 <- MARSSparamCIs(leg5_no3, alpha = 0.01)
noleg_no395 <- MARSSparamCIs(noleg_no3, alpha = 0.05)
leg1y_no395 <- MARSSparamCIs(leg1_no3, alpha = 0.05)
leg2y_no395 <- MARSSparamCIs(leg2_no3, alpha = 0.05)
leg3y_no395 <- MARSSparamCIs(leg3_no3, alpha = 0.05)
leg4y_no395 <- MARSSparamCIs(leg4_no3, alpha = 0.05)
leg5y_no395 <- MARSSparamCIs(leg5_no3, alpha = 0.05)

noleg_po499 <- MARSSparamCIs(noleg_po4, alpha = 0.01)
leg1y_po499 <- MARSSparamCIs(leg1_po4, alpha = 0.01)
leg2y_po499 <- MARSSparamCIs(leg2_po4, alpha = 0.01)
leg3y_po499 <- MARSSparamCIs(leg3_po4, alpha = 0.01)
leg4y_po499 <- MARSSparamCIs(leg4_po4, alpha = 0.01)
leg5y_po499 <- MARSSparamCIs(leg5_po4, alpha = 0.01)
noleg_po495 <- MARSSparamCIs(noleg_po4, alpha = 0.05)
leg1y_po495 <- MARSSparamCIs(leg1_po4, alpha = 0.05)
leg2y_po495 <- MARSSparamCIs(leg2_po4, alpha = 0.05)
leg3y_po495 <- MARSSparamCIs(leg3_po4, alpha = 0.05)
leg4y_po495 <- MARSSparamCIs(leg4_po4, alpha = 0.05)
leg5y_po495 <- MARSSparamCIs(leg5_po4, alpha = 0.05)

# Format confidence intervals into dataframes
noleg_nh499 = data.frame(
  "Est." = noleg_nh499$par$U,
  "Lower99" = noleg_nh499$par.lowCI$U,
  "Upper99" = noleg_nh499$par.upCI$U)
noleg_nh499$Parameter = rownames(noleg_nh499)
noleg_nh499[,1:3] = round(noleg_nh499[,1:3], 3)
noleg_nh499$Model = "immediate duration"
noleg_nh499$Nutrient = "NH4"

noleg_nh495 = data.frame(
  "Est." = noleg_nh495$par$U,
  "Lower95" = noleg_nh495$par.lowCI$U,
  "Upper95" = noleg_nh495$par.upCI$U)
noleg_nh495$Parameter = rownames(noleg_nh495)
noleg_nh495[,1:3] = round(noleg_nh495[,1:3], 3)
noleg_nh495$Model = "immediate duration"
noleg_nh495$Nutrient = "NH4"

noleg_no399 = data.frame(
  "Est." = noleg_no399$par$U,
  "Lower99" = noleg_no399$par.lowCI$U,
  "Upper99" = noleg_no399$par.upCI$U)
noleg_no399$Parameter = rownames(noleg_no399)
noleg_no399[,1:3] = round(noleg_no399[,1:3], 3)
noleg_no399$Model = "immediate duration"
noleg_no399$Nutrient = "NO3"

noleg_no395 = data.frame(
  "Est." = noleg_no395$par$U,
  "Lower95" = noleg_no395$par.lowCI$U,
  "Upper95" = noleg_no395$par.upCI$U)
noleg_no395$Parameter = rownames(noleg_no395)
noleg_no395[,1:3] = round(noleg_no395[,1:3], 3)
noleg_no395$Model = "immediate duration"
noleg_no395$Nutrient = "NO3"

noleg_po499 = data.frame(
  "Est." = noleg_po499$par$U,
  "Lower99" = noleg_po499$par.lowCI$U,
  "Upper99" = noleg_po499$par.upCI$U)
noleg_po499$Parameter = rownames(noleg_po499)
noleg_po499[,1:3] = round(noleg_po499[,1:3], 3)
noleg_po499$Model = "immediate duration"
noleg_po499$Nutrient = "PO4"

noleg_po495 = data.frame(
  "Est." = noleg_po495$par$U,
  "Lower95" = noleg_po495$par.lowCI$U,
  "Upper95" = noleg_po495$par.upCI$U)
noleg_po495$Parameter = rownames(noleg_po495)
noleg_po495[,1:3] = round(noleg_po495[,1:3], 3)
noleg_po495$Model = "immediate duration"
noleg_po495$Nutrient = "PO4"

leg1y_nh499 = data.frame(
  "Est." = leg1y_nh499$par$U,
  "Lower99" = leg1y_nh499$par.lowCI$U,
  "Upper99" = leg1y_nh499$par.upCI$U)
leg1y_nh499$Parameter = rownames(leg1y_nh499)
leg1y_nh499[,1:3] = round(leg1y_nh499[,1:3], 3)
leg1y_nh499$Model = "1 year duration"
leg1y_nh499$Nutrient = "NH4"

leg1y_nh495 = data.frame(
  "Est." = leg1y_nh495$par$U,
  "Lower95" = leg1y_nh495$par.lowCI$U,
  "Upper95" = leg1y_nh495$par.upCI$U)
leg1y_nh495$Parameter = rownames(leg1y_nh495)
leg1y_nh495[,1:3] = round(leg1y_nh495[,1:3], 3)
leg1y_nh495$Model = "1 year duration"
leg1y_nh495$Nutrient = "NH4"

leg1y_no399 = data.frame(
  "Est." = leg1y_no399$par$U,
  "Lower99" = leg1y_no399$par.lowCI$U,
  "Upper99" = leg1y_no399$par.upCI$U)
leg1y_no399$Parameter = rownames(leg1y_no399)
leg1y_no399[,1:3] = round(leg1y_no399[,1:3], 3)
leg1y_no399$Model = "1 year duration"
leg1y_no399$Nutrient = "NO3"

leg1y_no395 = data.frame(
  "Est." = leg1y_no395$par$U,
  "Lower95" = leg1y_no395$par.lowCI$U,
  "Upper95" = leg1y_no395$par.upCI$U)
leg1y_no395$Parameter = rownames(leg1y_no395)
leg1y_no395[,1:3] = round(leg1y_no395[,1:3], 3)
leg1y_no395$Model = "1 year duration"
leg1y_no395$Nutrient = "NO3"

leg1y_po499 = data.frame(
  "Est." = leg1y_po499$par$U,
  "Lower99" = leg1y_po499$par.lowCI$U,
  "Upper99" = leg1y_po499$par.upCI$U)
leg1y_po499$Parameter = rownames(leg1y_po499)
leg1y_po499[,1:3] = round(leg1y_po499[,1:3], 3)
leg1y_po499$Model = "1 year duration"
leg1y_po499$Nutrient = "PO4"

leg1y_po495 = data.frame(
  "Est." = leg1y_po495$par$U,
  "Lower95" = leg1y_po495$par.lowCI$U,
  "Upper95" = leg1y_po495$par.upCI$U)
leg1y_po495$Parameter = rownames(leg1y_po495)
leg1y_po495[,1:3] = round(leg1y_po495[,1:3], 3)
leg1y_po495$Model = "1 year duration"
leg1y_po495$Nutrient = "PO4"

leg2y_nh499 = data.frame(
  "Est." = leg2y_nh499$par$U,
  "Lower99" = leg2y_nh499$par.lowCI$U,
  "Upper99" = leg2y_nh499$par.upCI$U)
leg2y_nh499$Parameter = rownames(leg2y_nh499)
leg2y_nh499[,1:3] = round(leg2y_nh499[,1:3], 3)
leg2y_nh499$Model = "2 year duration"
leg2y_nh499$Nutrient = "NH4"

leg2y_nh495 = data.frame(
  "Est." = leg2y_nh495$par$U,
  "Lower95" = leg2y_nh495$par.lowCI$U,
  "Upper95" = leg2y_nh495$par.upCI$U)
leg2y_nh495$Parameter = rownames(leg2y_nh495)
leg2y_nh495[,1:3] = round(leg2y_nh495[,1:3], 3)
leg2y_nh495$Model = "2 year duration"
leg2y_nh495$Nutrient = "NH4"

leg2y_no399 = data.frame(
  "Est." = leg2y_no399$par$U,
  "Lower99" = leg2y_no399$par.lowCI$U,
  "Upper99" = leg2y_no399$par.upCI$U)
leg2y_no399$Parameter = rownames(leg2y_no399)
leg2y_no399[,1:3] = round(leg2y_no399[,1:3], 3)
leg2y_no399$Model = "2 year duration"
leg2y_no399$Nutrient = "NO3"

leg2y_no395 = data.frame(
  "Est." = leg2y_no395$par$U,
  "Lower95" = leg2y_no395$par.lowCI$U,
  "Upper95" = leg2y_no395$par.upCI$U)
leg2y_no395$Parameter = rownames(leg2y_no395)
leg2y_no395[,1:3] = round(leg2y_no395[,1:3], 3)
leg2y_no395$Model = "2 year duration"
leg2y_no395$Nutrient = "NO3"

leg2y_po499 = data.frame(
  "Est." = leg2y_po499$par$U,
  "Lower99" = leg2y_po499$par.lowCI$U,
  "Upper99" = leg2y_po499$par.upCI$U)
leg2y_po499$Parameter = rownames(leg2y_po499)
leg2y_po499[,1:3] = round(leg2y_po499[,1:3], 3)
leg2y_po499$Model = "2 year duration"
leg2y_po499$Nutrient = "PO4"

leg2y_po495 = data.frame(
  "Est." = leg2y_po495$par$U,
  "Lower95" = leg2y_po495$par.lowCI$U,
  "Upper95" = leg2y_po495$par.upCI$U)
leg2y_po495$Parameter = rownames(leg2y_po495)
leg2y_po495[,1:3] = round(leg2y_po495[,1:3], 3)
leg2y_po495$Model = "2 year duration"
leg2y_po495$Nutrient = "PO4"

leg3y_nh499 = data.frame(
  "Est." = leg3y_nh499$par$U,
  "Lower99" = leg3y_nh499$par.lowCI$U,
  "Upper99" = leg3y_nh499$par.upCI$U)
leg3y_nh499$Parameter = rownames(leg3y_nh499)
leg3y_nh499[,1:3] = round(leg3y_nh499[,1:3], 3)
leg3y_nh499$Model = "3 year duration"
leg3y_nh499$Nutrient = "NH4"

leg3y_nh495 = data.frame(
  "Est." = leg3y_nh495$par$U,
  "Lower95" = leg3y_nh495$par.lowCI$U,
  "Upper95" = leg3y_nh495$par.upCI$U)
leg3y_nh495$Parameter = rownames(leg3y_nh495)
leg3y_nh495[,1:3] = round(leg3y_nh495[,1:3], 3)
leg3y_nh495$Model = "3 year duration"
leg3y_nh495$Nutrient = "NH4"

leg3y_no399 = data.frame(
  "Est." = leg3y_no399$par$U,
  "Lower99" = leg3y_no399$par.lowCI$U,
  "Upper99" = leg3y_no399$par.upCI$U)
leg3y_no399$Parameter = rownames(leg3y_no399)
leg3y_no399[,1:3] = round(leg3y_no399[,1:3], 3)
leg3y_no399$Model = "3 year duration"
leg3y_no399$Nutrient = "NO3"

leg3y_no395 = data.frame(
  "Est." = leg3y_no395$par$U,
  "Lower95" = leg3y_no395$par.lowCI$U,
  "Upper95" = leg3y_no395$par.upCI$U)
leg3y_no395$Parameter = rownames(leg3y_no395)
leg3y_no395[,1:3] = round(leg3y_no395[,1:3], 3)
leg3y_no395$Model = "3 year duration"
leg3y_no395$Nutrient = "NO3"

leg3y_po499 = data.frame(
  "Est." = leg3y_po499$par$U,
  "Lower99" = leg3y_po499$par.lowCI$U,
  "Upper99" = leg3y_po499$par.upCI$U)
leg3y_po499$Parameter = rownames(leg3y_po499)
leg3y_po499[,1:3] = round(leg3y_po499[,1:3], 3)
leg3y_po499$Model = "3 year duration"
leg3y_po499$Nutrient = "PO4"

leg3y_po495 = data.frame(
  "Est." = leg3y_po495$par$U,
  "Lower95" = leg3y_po495$par.lowCI$U,
  "Upper95" = leg3y_po495$par.upCI$U)
leg3y_po495$Parameter = rownames(leg3y_po495)
leg3y_po495[,1:3] = round(leg3y_po495[,1:3], 3)
leg3y_po495$Model = "3 year duration"
leg3y_po495$Nutrient = "PO4"

leg4y_nh499 = data.frame(
  "Est." = leg4y_nh499$par$U,
  "Lower99" = leg4y_nh499$par.lowCI$U,
  "Upper99" = leg4y_nh499$par.upCI$U)
leg4y_nh499$Parameter = rownames(leg4y_nh499)
leg4y_nh499[,1:3] = round(leg4y_nh499[,1:3], 3)
leg4y_nh499$Model = "4 year duration"
leg4y_nh499$Nutrient = "NH4"

leg4y_nh495 = data.frame(
  "Est." = leg4y_nh495$par$U,
  "Lower95" = leg4y_nh495$par.lowCI$U,
  "Upper95" = leg4y_nh495$par.upCI$U)
leg4y_nh495$Parameter = rownames(leg4y_nh495)
leg4y_nh495[,1:3] = round(leg4y_nh495[,1:3], 3)
leg4y_nh495$Model = "4 year duration"
leg4y_nh495$Nutrient = "NH4"

leg4y_no399 = data.frame(
  "Est." = leg4y_no399$par$U,
  "Lower99" = leg4y_no399$par.lowCI$U,
  "Upper99" = leg4y_no399$par.upCI$U)
leg4y_no399$Parameter = rownames(leg4y_no399)
leg4y_no399[,1:3] = round(leg4y_no399[,1:3], 3)
leg4y_no399$Model = "4 year duration"
leg4y_no399$Nutrient = "NO3"

leg4y_no395 = data.frame(
  "Est." = leg4y_no395$par$U,
  "Lower95" = leg4y_no395$par.lowCI$U,
  "Upper95" = leg4y_no395$par.upCI$U)
leg4y_no395$Parameter = rownames(leg4y_no395)
leg4y_no395[,1:3] = round(leg4y_no395[,1:3], 3)
leg4y_no395$Model = "4 year duration"
leg4y_no395$Nutrient = "NO3"

leg4y_po499 = data.frame(
  "Est." = leg4y_po499$par$U,
  "Lower99" = leg4y_po499$par.lowCI$U,
  "Upper99" = leg4y_po499$par.upCI$U)
leg4y_po499$Parameter = rownames(leg4y_po499)
leg4y_po499[,1:3] = round(leg4y_po499[,1:3], 3)
leg4y_po499$Model = "4 year duration"
leg4y_po499$Nutrient = "PO4"

leg4y_po495 = data.frame(
  "Est." = leg4y_po495$par$U,
  "Lower95" = leg4y_po495$par.lowCI$U,
  "Upper95" = leg4y_po495$par.upCI$U)
leg4y_po495$Parameter = rownames(leg4y_po495)
leg4y_po495[,1:3] = round(leg4y_po495[,1:3], 3)
leg4y_po495$Model = "4 year duration"
leg4y_po495$Nutrient = "PO4"

leg5y_nh499 = data.frame(
  "Est." = leg5y_nh499$par$U,
  "Lower99" = leg5y_nh499$par.lowCI$U,
  "Upper99" = leg5y_nh499$par.upCI$U)
leg5y_nh499$Parameter = rownames(leg5y_nh499)
leg5y_nh499[,1:3] = round(leg5y_nh499[,1:3], 3)
leg5y_nh499$Model = "5 year duration"
leg5y_nh499$Nutrient = "NH4"

leg5y_nh495 = data.frame(
  "Est." = leg5y_nh495$par$U,
  "Lower95" = leg5y_nh495$par.lowCI$U,
  "Upper95" = leg5y_nh495$par.upCI$U)
leg5y_nh495$Parameter = rownames(leg5y_nh495)
leg5y_nh495[,1:3] = round(leg5y_nh495[,1:3], 3)
leg5y_nh495$Model = "5 year duration"
leg5y_nh495$Nutrient = "NH4"

leg5y_no399 = data.frame(
  "Est." = leg5y_no399$par$U,
  "Lower99" = leg5y_no399$par.lowCI$U,
  "Upper99" = leg5y_no399$par.upCI$U)
leg5y_no399$Parameter = rownames(leg5y_no399)
leg5y_no399[,1:3] = round(leg5y_no399[,1:3], 3)
leg5y_no399$Model = "5 year duration"
leg5y_no399$Nutrient = "NO3"

leg5y_no395 = data.frame(
  "Est." = leg5y_no395$par$U,
  "Lower95" = leg5y_no395$par.lowCI$U,
  "Upper95" = leg5y_no395$par.upCI$U)
leg5y_no395$Parameter = rownames(leg5y_no395)
leg5y_no395[,1:3] = round(leg5y_no395[,1:3], 3)
leg5y_no395$Model = "5 year duration"
leg5y_no395$Nutrient = "NO3"

leg5y_po499 = data.frame(
  "Est." = leg5y_po499$par$U,
  "Lower99" = leg5y_po499$par.lowCI$U,
  "Upper99" = leg5y_po499$par.upCI$U)
leg5y_po499$Parameter = rownames(leg5y_po499)
leg5y_po499[,1:3] = round(leg5y_po499[,1:3], 3)
leg5y_po499$Model = "5 year duration"
leg5y_po499$Nutrient = "PO4"

leg5y_po495 = data.frame(
  "Est." = leg5y_po495$par$U,
  "Lower95" = leg5y_po495$par.lowCI$U,
  "Upper95" = leg5y_po495$par.upCI$U)
leg5y_po495$Parameter = rownames(leg5y_po495)
leg5y_po495[,1:3] = round(leg5y_po495[,1:3], 3)
leg5y_po495$Model = "5 year duration"
leg5y_po495$Nutrient = "PO4"

#### Join data ####

# Bind all together
nh4CIs99 = rbind(noleg_nh499, leg1y_nh499, leg2y_nh499, 
              leg3y_nh499, leg4y_nh499,leg5y_nh499)

nh4CIs95 = rbind(noleg_nh495, leg1y_nh495, leg2y_nh495, 
              leg3y_nh495, leg4y_nh495,leg5y_nh495)

no3CIs99 = rbind(noleg_no399, leg1y_no399, leg2y_no399, 
                 leg3y_no399, leg4y_no399,leg5y_no399)

no3CIs95 = rbind(noleg_no395, leg1y_no395, leg2y_no395, 
                 leg3y_no395, leg4y_no395,leg5y_no395)

po4CIs99 = rbind(noleg_po499, leg1y_po499, leg2y_po499, 
                 leg3y_po499, leg4y_po499,leg5y_po499)

po4CIs95 = rbind(noleg_po495, leg1y_po495, leg2y_po495, 
                 leg3y_po495, leg4y_po495,leg5y_po495)

# Add column for site names
nh4CIs99$Stream = gsub("_","",str_sub(nh4CIs99$Parameter, start= -4))
nh4CIs95$Stream = gsub("_","",str_sub(nh4CIs95$Parameter, start= -4))
no3CIs99$Stream = gsub("_","",str_sub(no3CIs99$Parameter, start= -4))
no3CIs95$Stream = gsub("_","",str_sub(no3CIs95$Parameter, start= -4))
po4CIs99$Stream = gsub("_","",str_sub(po4CIs99$Parameter, start= -4))
po4CIs95$Stream = gsub("_","",str_sub(po4CIs95$Parameter, start= -4))

# Join data
nh4CIs <- full_join(nh4CIs99, nh4CIs95)
no3CIs <- full_join(no3CIs99, no3CIs95)
po4CIs <- full_join(po4CIs99, po4CIs95)

nCIs <- full_join(nh4CIs, no3CIs)
nutCIs <- full_join(nCIs, po4CIs)

# Simplify parameter names
nutCIs$Parm_simple = c(rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4),
                    
                    rep("Ppt",4),
                    rep("Perc. burn",4),
                    rep("Ppt x Perc. burn",4))

#### Plot all together ####

# Add column to designate those sites at which effects are significant.
nutCIs <- nutCIs %>%
  mutate(sig = factor(case_when(`Est.` > 0 & 
                                  Lower95 > 0 & Lower99 > 0 & 
                                  Upper95 > 0 & Upper99 > 0 ~ "sig_pos",
                                `Est.` < 0 & 
                                  Lower95 < 0 & Lower99 < 0 &
                                  Upper95 < 0 & Upper99 < 0 ~ "sig_neg",
                                `Est.` > 0 & 
                                  Lower95 > 0 & 
                                  Upper95 > 0  ~ "weak_sig_pos",
                                `Est.` < 0 & 
                                  Lower95 < 0 &
                                  Upper95 < 0 ~ "weak_sig_neg",
                                TRUE ~ "not_sig"), 
                      levels = c("weak_sig_pos", "sig_pos", 
                                 "not_sig", 
                                 "sig_neg", "weak_sig_neg"))) %>%
  mutate(model = factor(Model, levels = c("immediate duration",
                                          "1 year duration",
                                          "2 year duration",
                                          "3 year duration",
                                          "4 year duration",
                                          "5 year duration")))

(all_nut_fig <- ggplot(nutCIs, aes(x = factor(Parm_simple, 
                                       levels = c("Ppt x Perc. burn",
                                                  "Perc. burn",
                                                  "Ppt")),
                            y = Est., fill = sig, shape = Stream)) + 
    geom_errorbar(aes(ymin = Lower99, ymax = Upper99),
                  position=position_dodge(width = 0.5), width = 0) +
    geom_point(position=position_dodge(width = 0.5), 
               alpha = 0.8, size = 8) + 
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    scale_fill_manual(values = c("gray50", "black", "white", "black", "gray50")) +
    theme_bw()+
    theme(axis.title = element_text(size = 24),
          axis.text = element_text(size = 20),
          strip.text.x = element_text(size = 24),
          strip.text.y = element_text(size = 24),
          legend.title=element_text(size = 20), 
          legend.text=element_text(size = 20)) +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    coord_flip(ylim = c(-1.0, 1.0)) + 
    labs(y = "Effect Size", 
         x = "Covariates",
         fill = "Significance") +
    theme(plot.margin = unit(c(.2,.2,.05,.05),"cm")) + 
    guides(shape = guide_legend("Stream"), fill = "none") +
    facet_grid(Nutrient~model))

# Export plot.
# ggsave(("MARSS_nuts_111523.png"),
#        path = "figures",
#        width = 65,
#        height = 36,
#        units = "cm"
# )

# End of script.
