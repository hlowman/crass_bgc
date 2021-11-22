########################################################################################################################################################################################
#Author: Adrianne Smits
#Date: October 18, 2017
#Stream Resiliency RCN, Large Scale Drivers Group
#Title: 4Basins models
#Description: This script runs multivariate state-space models on annual flow and solute load time series from the four major
#tributaries in the Mississippi-Atchafalaya River Basin: the upper Mississippi, Missouri, Ohio-Tennessee, and Arkansas-Red
# All data used in this analysis are publically available
############################################################################################################################################################################################

##Load MARSS package
library(MARSS)

##Set working directory (tell R where your files are)

#######################################################################################################################################################################################################
##Load clean, z-score standardized data:

#Z-scored observation time series data for response variable of interest (flow or a solute)
###INPUT DESIRED RESPONSE VARIABLE
dat.z <- as.matrix(read.csv('02_TN_z.csv', row.names=1))###INPUT DESIRED RESPONSE VARIABLE
dim(dat.z)#five time series, 36 years of data each
dat.z
#number of different time series
nsites <- nrow(dat.z)

##Load z-scored covariate data, check to see if it's z-scored and correct
covariates.z <-as.matrix(read.csv('02_Covariates_z.csv', row.names=1))
row.names(covariates.z)
#variance should be 1 for all the covariates
apply(covariates.z,1,var)

##Load the estimates of observation error for each site (R matrix),
#These estimates are derived from LOADEST 95% confidence intervals and scaled to the total variance of each time series
R.est <- read.csv('02_4Basins_R_estimates_usingMedian.csv',header=TRUE)
R.est <- R.est[,-1]

######################################################################################################################################################################################################
##Create inputs to MARSS function (matrices):

#Create names of models to run: (list) ###UPDATE THIS WHEN TRYING NEW SETS OF MODELS!!!!!!
#Name Format: 'mod', first number, '.','second number','first letter','_','second letters'

#Note: First number indicates number of states (Example: mod1.1 estimates one underyling state process, mod3.1 indicates three state processes)
#Note: Second number indicates which covariates are included in model (Example: mod4.2 has four estimated state processes, and includes SAT as a covariate)
#       Covariate data that is area-weighted (temperature, precip, PDSI, crop acreage) will be used to estimate a single effect
#       for covariate data that is very large scale (SAT, MEI, NAO), we will test models that include either uniform effects of the covariate
#       or site-specific effects (ex. the effect of MEI on total N loads will be estimated separately for UpperMiss, OhioTenn, ArkRed, and Missouri) 
  #Covariate Codes:
  #1 = no covariates
  #2 = SAT (annual land surface air temperature anomaly for the Northern Hemisphere)
  #3 = NAO (annual North Atlantic Oscillation index)
  #4 = Max MEI (annual maximum value of multivariate ENSO index)
  #5 = Crop acres (wheat + corn, area-weighted)
  #6 = PDSI (Palmer Drought Severity Index,area weighted)
  #7 = Precipitation (annual total, area weighted)
  #8 = Temperature (annual average, area-weighted)
  #9 = Max MEI lagged by one year 
  #10 = NAO lagged by one year
  
#Note: First Letter indicates the process error structure in the model (environmental variability):
  #Process error codes
  #e = equal variance and covariance
  #du= diagonal and unequal
  #de = diagonal and equal

#Note: The second set of letters indicates whether to estimate a single effect of a covariate for all the state processes, 
#or whether to estimate state-specific covariate effects
  #Effect codes
  #sh= shared covariate effect among all states
  #sep= separate covariate effects estimated for each state

#UPDATE THIS
mod.names <- c("mod1.1", "mod1.2", "mod1.3","mod1.4",
               "mod1.5",
               "mod1.6", "mod1.7", "mod1.8", "mod1.9", "mod1.10",
               
               "mod2.1e","mod2.1du","mod2.1de", 
               "mod2.2e_sh","mod2.2du_sh","mod2.2de_sh", 
               "mod2.2e_sep","mod2.2du_sep","mod2.2de_sep", 
               "mod2.3e_sh","mod2.3du_sh","mod2.3de_sh", 
               "mod2.3e_sep","mod2.3du_sep","mod2.3de_sep",
               "mod2.4e_sh","mod2.4du_sh","mod2.4de_sh",  
               "mod2.4e_sep","mod2.4du_sep","mod2.4de_sep",  
                "mod2.5e_sh","mod2.5du_sh","mod2.5de_sh", 
               "mod2.5e_sep","mod2.5du_sep","mod2.5de_sep", 
               "mod2.6e_sh","mod2.6du_sh","mod2.6de_sh", 
               "mod2.7e_sh","mod2.7du_sh","mod2.7de_sh", 
               "mod2.8e_sh","mod2.8du_sh","mod2.8de_sh", 
               "mod2.9e_sh","mod2.9du_sh","mod2.9de_sh", 
               "mod2.9e_sep","mod2.9du_sep","mod2.9de_sep",
               "mod2.10e_sh","mod2.10du_sh","mod2.10de_sh", 
               "mod2.10e_sep","mod2.10du_sep","mod2.10de_sep",
               
               "mod3.1e","mod3.1du","mod3.1de", 
               "mod3.2e_sh","mod3.2du_sh","mod3.2de_sh", 
               "mod3.2e_sep","mod3.2du_sep","mod3.2de_sep", 
               "mod3.3e_sh","mod3.3du_sh","mod3.3de_sh", 
               "mod3.3e_sep","mod3.3du_sep","mod3.3de_sep",
               "mod3.4e_sh","mod3.4du_sh","mod3.4de_sh",  
               "mod3.4e_sep","mod3.4du_sep","mod3.4de_sep",  
               "mod3.5e_sh","mod3.5du_sh","mod3.5de_sh", 
              "mod3.5e_sep","mod3.5du_sep","mod3.5de_sep", 
               "mod3.6e_sh","mod3.6du_sh","mod3.6de_sh", 
               "mod3.7e_sh","mod3.7du_sh","mod3.7de_sh", 
               "mod3.8e_sh","mod3.8du_sh","mod3.8de_sh",
               "mod3.9e_sh","mod3.9du_sh","mod3.9de_sh", 
               "mod3.9e_sep","mod3.9du_sep","mod3.9de_sep",
               "mod3.10e_sh","mod3.10du_sh","mod3.10de_sh", 
               "mod3.10e_sep","mod3.10du_sep","mod3.10de_sep",
               
               "mod4.1e","mod4.1du","mod4.1de", 
               "mod4.2e_sh","mod4.2du_sh","mod4.2de_sh", 
               "mod4.2e_sep","mod4.2du_sep","mod4.2de_sep", 
               "mod4.3e_sh","mod4.3du_sh","mod4.3de_sh", 
               "mod4.3e_sep","mod4.3du_sep","mod4.3de_sep",
               "mod4.4e_sh","mod4.4du_sh","mod4.4de_sh",  
               "mod4.4e_sep","mod4.4du_sep","mod4.4de_sep", 
               "mod4.5e_sh","mod4.5du_sh","mod4.5de_sh", 
               "mod4.5e_sep","mod4.5du_sep","mod4.5de_sep",
               "mod4.6e_sh","mod4.6du_sh","mod4.6de_sh", 
               "mod4.7e_sh","mod4.7du_sh","mod4.7de_sh", 
               "mod4.8e_sh","mod4.8du_sh","mod4.8de_sh",
               "mod4.9e_sh","mod4.9du_sh","mod4.9de_sh", 
               "mod4.9e_sep","mod4.9du_sep","mod4.9de_sep",
               "mod4.10e_sh","mod4.10du_sh","mod4.10de_sh", 
               "mod4.10e_sep","mod4.10du_sep","mod4.10de_sep")

#Number of models
nmods <- length(mod.names)

#Names of matrix inputs to MARSS equation:
mat.names <- c('B', 'Q','Z', 'R','A','U','c','C')

#Empty list to put model outputs into:
mod.output <- list()

##Create vectors, matrices, or lists of parameter matrices:

##B (degree of mean reversion)
B.list <- c('diagonal and equal' ,'diagonal and unequal' )

##Q (process error matrix)
#equalvarcov = all states have equal variance, and equal covariance with one another
#diagonal and unequal = all states have different variances, no covariance among states
#diagonal and equal = all states have the same variance, but no covariance
Q.list <- c('equalvarcov', 'diagonal and unequal','diagonal and equal')
Q.list[[3]]

##Z (number of estimated state processes)
#one state:
Z.1 <- matrix(1,4,1)
#Two states: OhioTenn is its own state, all others are grouped into one state
Z.2 <- factor(c("upperMiss","OhioTenn","upperMiss","upperMiss"))

#Three states: UpperMiss and Missouri basins are lumped into one state (u),
#but OhioTenn and ArkRed remain separate states
Z.3 <- factor(c("upperMiss","OhioTenn", "upperMiss", "ArkRed"))

#Four states (each basin is its own state)
Z.4 <- 'identity'

#list of Z matrix structures
Z.list <- list(Z.1, Z.2, Z.3, Z.4)
Z.list[[3]]

##R (observation error matrix)
#set obs. error to a fixed value for each basin. Choose the value based on solute!!!!Change this!!!!!!!!!
R.fixed <- matrix(list(0),nsites,nsites)
diag(R.fixed) <-R.est[which(R.est$solutes=='TN'),2:5] #MAKE SURE TO CHANGE THE SOLUTE!!!!
R.list <- list('diagonal and equal', 'diagonal and unequal', R.fixed)
R.list[[3]]

##A (not estimated in this analysis)
A.list <- 'zero'

##U (not estimated in this analysis)
U.list <- 'zero'

#xO (initial conditions)
x0.model <- 'zero'

#VO (intitial conditions)
V0.model <- 'zero'

##c: UPDATE THIS!!!
#different sets of covariate matrices
row.names(covariates.z)

#no covariates
nocovar <- matrix(0)
#agricultural variables
  #basin-specific area-weighted crop data (corn + wheat acres planted)
basin_crop_acres <- matrix(covariates.z[59:68,],nrow=10)#ArkRed, Missouri, OhioTenn, UpperMiss corn+wheat acres
row.names(basin_crop_acres) <- row.names(covariates.z[59:68,])
  #different groupings of crop data
All4_tog <- matrix(basin_crop_acres[6,],nrow=1)
OT_alone <- matrix(basin_crop_acres[c(8,4),],nrow=2)
OT_Ark_alone <- matrix(basin_crop_acres[c(9,4,1),],nrow=3)
All4_sep <- matrix(basin_crop_acres[c(5,4,3,1),],nrow=4)
  #Climate variables (no lags)
temp_anomaly <- matrix(covariates.z[11,],nrow=1)# Mean annual air surface temperature anomaly over Northern Hemisphere
NAO <- matrix(covariates.z[1,],nrow=1)#North Atlantic Oscillation index annual
Max_MEI <- matrix(covariates.z[3,],nrow=1)#Yearly maximum MEI- large positives= El Nino
PDSI_all4tog <- matrix(covariates.z[81,],nrow=1) #PDSI
PDSI_OT_alone <- matrix(covariates.z[c(84,75),],nrow=2)
PDSI_OT_ARK_alone <- matrix(covariates.z[c(87,75,73),],nrow=3)
PDSI_All4sep <- matrix(covariates.z[c(76,75,74,73),],nrow=4)
PPT_all4tog <- matrix(covariates.z[82,],nrow=1) #precipitation area-weighted
PPT_OT_alone <- matrix(covariates.z[c(85,71),],nrow=2)
PPT_OT_ARK_alone <- matrix(covariates.z[c(88,71,69),],nrow=3)
PPT_All4sep <- matrix(covariates.z[c(72,71,70,69),],nrow=4)
Temp_all4tog <- matrix(covariates.z[83,],nrow=1) #annual average temperature area-weighted
Temp_OT_alone <- matrix(covariates.z[c(86,79),],nrow=2)
Temp_OT_ARK_alone <- matrix(covariates.z[c(89,79,77),],nrow=3)
Temp_All4sep <- matrix(covariates.z[c(80,79,78,77),],nrow=4)
  #climate variables (one year lag)
Max_MEI_lag1 <- matrix(covariates.z[14,],nrow=1)
NAO_lag1 <- matrix(covariates.z[12,],nrow=1)
  
  
#List of all covariate matrices....UPDATE THIS!!
c.list <- list(nocovar, 
               All4_tog,OT_alone,OT_Ark_alone, All4_sep,#crop acres (don't include as a covariate for flow models)
               temp_anomaly, 
               NAO, 
               Max_MEI, 
               PDSI_all4tog, PDSI_OT_alone, PDSI_OT_ARK_alone, PDSI_All4sep,
               PPT_all4tog, PPT_OT_alone, PPT_OT_ARK_alone, PPT_All4sep,
               Temp_all4tog, Temp_OT_alone,Temp_OT_ARK_alone,Temp_All4sep,
               Max_MEI_lag1, 
               NAO_lag1)
c.list[[10]]

##C (matrix that maps covariates to the state variables)
#one state
C_1.1 <- matrix(0)
C_1.2 <- matrix("SAT",1,1)
C_1.3 <- matrix("NAO",1,1)
C_1.4 <- matrix("MEI",1,1)
C_1.5 <- matrix("crop_acres",1,1)
C_1.6 <- matrix("PDSI",1,1)
C_1.7 <- matrix("precip",1,1)
C_1.8 <- matrix("temp",1,1)
C_1.9 <- matrix("MEI_lag1",1,1)
C_1.10 <- matrix("NAO_lag1",1,1)

# 2 states
C_2.1 <- matrix(0,nrow=2)
C_2.2_sh <-  matrix(c("SAT"),2,1)
C_2.2_sep <- matrix(c("SAT_UpperMiss","SAT_OhioTenn"),2,1)
C_2.3_sh <- matrix(c("NAO"),2,1)
C_2.3_sep <- matrix(c("NAO_UpperMiss","NAO_OhioTenn"),2,1)
C_2.4_sh <- matrix(c("MEI"),2,1)
C_2.4_sep <- matrix(c("MEI_UpperMiss","MEI_OhioTenn"),2,1)
 C_2.5_sh <- matrix(list(0),2,2)
        diag(C_2.5_sh) <- "crop_acres"
 C_2.5_sep <- matrix(list(0),2,2)
        diag(C_2.5_sep) <- c("crops_UpperMiss","crops_OhioTenn")
C_2.6_sh <- matrix(list(0),2,2)
diag(C_2.6_sh) <- "PDSI"
C_2.7_sh <- matrix(list(0),2,2)
diag(C_2.7_sh) <- "precip"
C_2.8_sh <-matrix(list(0),2,2)
diag(C_2.8_sh) <- "temp"
C_2.9_sh <- matrix(c("MEI_lag1"),2,1)
C_2.9_sep <- matrix(c("MEI_lag1_UpperMiss","MEI_lag1_OhioTenn"),2,1)
C_2.10_sh <- matrix(c("NAO_lag1"),2,1)
C_2.10_sep <- matrix(c("NAO_lag1_UpperMiss","NAO_lag1_OhioTenn"),2,1)


# 3 states
C_3.1 <- matrix(0,nrow=3)
C_3.2_sh <- matrix("SAT",3,1)
C_3.2_sep <- matrix(c("SAT_UpperMiss","SAT_OhioTenn","SAT_ArkRed"),3,1)
C_3.3_sh <- matrix("NAO",3,1)
C_3.3_sep <- matrix(c("NAO_UpperMiss","NAO_OhioTenn","NAO_ArkRed"),3,1)
C_3.4_sh <- matrix("MEI",3,1)
C_3.4_sep <- matrix(c("MEI_UpperMiss","MEI_OhioTenn","MEI_ArkRed"),3,1)
 C_3.5_sh <- matrix(list(0),3,3)
       diag(C_3.5_sh) <- "crop_acres"
C_3.5_sep <- matrix(list(0),3,3)
       diag(C_3.5_sep) <- c("crops_UpperMiss","crops_OhioTenn","crops_ArkRed")
C_3.6_sh <- matrix(list(0),3,3)
diag(C_3.6_sh) <- "PDSI"
C_3.7_sh <- matrix(list(0),3,3)
diag(C_3.7_sh) <- "precip"
C_3.8_sh <- matrix(list(0),3,3)
diag(C_3.8_sh) <- "temp"
C_3.9_sh <- matrix("MEI_lag1",3,1)
C_3.9_sep <- matrix(c("MEI_lag1_UpperMiss","MEI_lag1_OhioTenn","MEI_lag1_ArkRed"),3,1)
C_3.10_sh <- matrix("NAO_lag1",3,1)
C_3.10_sep <- matrix(c("NAO_lag1_UpperMiss","NAO_lag1_OhioTenn","NAO_lag1_ArkRed"),3,1)


# 4 states
C_4.1 <- matrix(0,nrow=4)
C_4.2_sh <- matrix("SAT",4,1)
C_4.2_sep <- matrix(c("SAT_UpperMiss","SAT_OhioTenn","SAT_MO","SAT_ArkRed"),4,1)
C_4.3_sh <- matrix("NAO",4,1)
C_4.3_sep <- matrix(c("NAO_UpperMiss","NAO_OhioTenn","NAO_MO","NAO_ArkRed"),4,1)
C_4.4_sh <- matrix("MEI",4,1)
C_4.4_sep <- matrix(c("MEI_UpperMiss","MEI_OhioTenn","MEI_MO","MEI_ArkRed"),4,1)
 C_4.5_sh <- matrix(list(0),4,4)
        diag(C_4.5_sh) <- "crop_acres"
 C_4.5_sep <- matrix(list(0),4,4)
       diag(C_4.5_sep) <- c("crops_UpperMiss","crops_OhioTenn","crops_MO","crops_ArkRed")
C_4.6_sh <- matrix(list(0),4,4)
diag(C_4.6_sh) <- "PDSI"
C_4.7_sh <-  matrix(list(0),4,4)
diag(C_4.7_sh) <- "precip"
C_4.8_sh <- matrix(list(0),4,4)
diag(C_4.8_sh) <- "temp"
C_4.9_sh <- matrix("MEI_lag1",4,1)
C_4.9_sep <- matrix(c("MEI_lag1_UpperMiss","MEI_lag1_OhioTenn","MEI_lag1_MO","MEI_lag1_ArkRed"),4,1)
C_4.10_sh <- matrix("NAO_lag1",4,1)
C_4.10_sep <- matrix(c("NAO_lag1_UpperMiss","NAO_lag1_OhioTenn","NAO_lag1_MO","NAO_lag1_ArkRed"),4,1)



#list of all C matrices ##UPDATE THIS!!!!!
C.list <- list(C_1.1, C_1.2, C_1.3, C_1.4, 
               C_1.5, 
               C_1.6, C_1.7, C_1.8, C_1.9, C_1.10,
               C_2.1, C_2.2_sh, C_2.2_sep ,C_2.3_sh, C_2.3_sep, C_2.4_sh, C_2.4_sep,
              C_2.5_sh, C_2.5_sep, 
               C_2.6_sh,  C_2.7_sh,  C_2.8_sh, C_2.9_sh,C_2.9_sep, C_2.10_sh, C_2.10_sep,
               C_3.1, C_3.2_sh, C_3.2_sep, C_3.3_sh, C_3.3_sep, C_3.4_sh, C_3.4_sep, 
              C_3.5_sh, C_3.5_sep,
               C_3.6_sh,  C_3.7_sh,  C_3.8_sh, C_3.9_sh, C_3.9_sep, C_3.10_sh, C_3.10_sep,
               C_4.1, C_4.2_sh, C_4.2_sep, C_4.3_sh, C_4.3_sep, C_4.4_sh, C_4.4_sep, 
              C_4.5_sh, C_4.5_sep,
               C_4.6_sh, C_4.7_sh, C_4.8_sh, C_4.9_sh, C_4.9_sep, C_4.10_sh, C_4.10_sep)         
C.list[[15]]

########################################################################################################################################################################################################################################
## Create a matrix of parameter combinations: (combos),
#nrow = number of models tested, 
#ncol= number of parameter matrices in MARSS equation
#WARNING: UPDATE THIS EVERY TIME YOU TEST NEW SETS OF MODELS OR SWITCH BETWEEN RUNNING
#A SOLUTE MODEL VERSUS A FLOW MODEL!!!!!

#Empty matrix
combos <- matrix(0,nmods, length(mat.names), dimnames=list(mod.names,mat.names))
#B
combos[,1] <- rep(1,nmods)
#Q
combos[,2] <-c(rep(1,10), rep(1:3,48)) 
#Z
combos[,3] <- c(rep(1,10),rep(2,48),rep(3,48),rep(4,48))
#R
combos[,4] <- rep(3,nmods)
#c:  which covariates to use
combos[,7] <- c(c(1,6,7,8,
                  2,
                  9,13,17,21,22), 
                c(rep(1,3),rep(6,6),rep(7,6),rep(8,6),
                  rep(3,6), 
                  rep(10,3), rep(14,3),rep(18,3),rep(21,6),rep(22,6)),
                c(rep(1,3),rep(6,6),rep(7,6),rep(8,6), 
                  rep(4,6),
                  rep(11,3), rep(15,3),rep(19,3),rep(21,6),rep(22,6)),
                c(rep(1,3),rep(6,6),rep(7,6),rep(8,6), 
                 rep(5,6),
                  rep(12,3), rep(16,3),rep(20,3),rep(21,6), rep(22,6)))
#C
combos[,8] <- c(1:4,
                5,
                6:10,
              rep(11,3), rep(12,3), rep(13,3),rep(14,3), rep(15,3),rep(16,3),rep(17,3),
              rep(18,3),rep(19,3),
              rep(20,3),rep(21,3),rep(22,3),rep(23,3), rep(24,3),rep(25,3),rep(26,3),
              rep(27,3),rep(28,3),rep(29,3), rep(30,3),rep(31,3),rep(32,3),rep(33,3),
              rep(34,3),rep(35,3),
              rep(36,3),rep(37,3),rep(38,3),rep(39,3), rep(40,3),rep(41,3),rep(42,3),
              rep(43,3),rep(44,3),rep(45,3),rep(46,3),rep(47,3),rep(48,3),rep(49,3),
              rep(50,3),rep(51,3),
              rep(52,3),rep(53,3),rep(54,3),rep(55,3),rep(56,3),rep(57,3),rep(58,3))

#check to make sure it's correct!
combos

####################################################################################################################################################################################################################################################
# This loop runs the MARSS function for all the combinations of model parameters contained in the matrix 'combos'
#and stores the output in the list 'mod.output'
for(i in 1:nrow(combos)){
  #select model structure and parameters
  #change this so that it uses combo matrix to index correct parameters
  mod.list <- list(B=B.list[combos[i,1]], Q=Q.list[combos[i,2]], Z=Z.list[[combos[i,3]]], R=R.list[[combos[i,4]]] ,
                   A=A.list, U=U.list, c=c.list[[combos[i,7]]], C=C.list[[combos[i,8]]], x0=x0.model, V0=V0.model, tinitx=0)
  #MARSS function call, can change number of iterations if desired
  mod <- MARSS(dat.z, model=mod.list, control=list(maxit=1000))
  #put MARSS output into a list
  mod.output[[i]] <- mod 
  #give the model output the right name
  names(mod.output)[i] <- mod.names[i]
}

##############################################################################################################################################################################################################################################################
#Check what's in mod.output
mod.output[[4]]
names(mod.output)

#############################################################################################################################################################################################################################
#Save model output as Rdata file, to use in other scripts:
save(mod.output, file="03_TN_MARSS_output_FixedR.Rdata")#CHANGE FILE NAME DEPENDING ON WHICH RESPONSE VARIABLE YOU'RE LOOKING AT
