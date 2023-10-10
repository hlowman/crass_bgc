{\rtf1\ansi\ansicpg1252\cocoartf2709
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fmodern\fcharset0 Courier;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\margl1440\margr1440\vieww14620\viewh18780\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs26 \cf0 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 ### Simulation demonstrating why covariances (Q matrix) among states are not indicative of coherence among states ###\
\
## Author: Mark Scheurell\
## Comments: TK Harms\
\
### Simulation:\
## 2 autoregressive states from a common discrete-time sine wave\
## Add negatively correlated process errors\
## States remain negatively correlated\
## Interpretation: A single large-scale synchronizing factor could be combined with small-scale perturbations. The covariances in Q represent marginal covariance in the states. A non-diagnonal covariance matrix with positive or negative values in the off-diagonals cannot be interpreted as states increasing/decreasing synchronously.\
\
set.seed(550)\
\
## number of time steps\
tt <- 48\
\
## period (months)\
pp <- 12\
\
## stationary autoregressive model (coefs = 0.6 & 0.8)\
BB <- matrix(c(0.6, 0, 0, 0.8), 2, 2)\
\
## common environmental driver (discrete sine wave)\
cc <- sin(2 * pi * seq(tt) / pp)\
\
## plot the sine wave\
plot.ts(cc)\
\
## common positive effect of the environment (= 5)\
CC <- matrix(c(5, 5), 2, 1)\
\
## mean of the process errors (= 0)\
mu <- c(0, 0)\
\
## covariance matrix for process errors\
sigma <- matrix(c(3, -2, -2, 3), 2, 2)\
\
## simulate process errors\
ww <- MASS::mvrnorm(n = tt,\
                    mu = mu,\
                    Sigma = sigma) |>\
  t()\
\
## initialize states\
xx <- ww\
\
## first state is function of sine wave and error\
xx[,1] <- CC %*% cc[1] + ww[,1]\
\
## calculate states 2:T\
for(t in 2:tt) \{\
  xx[,t] <- BB %*% xx[,t-1] + CC %*% cc[t] + ww[,t]\
\}\
\
## plot both states over time\
xx |>\
  t() |>\
  plot.ts()\
\
## plot both states against each other\
plot(xx[1,], xx[2,], pch = 16)\
\
## correlation of states (positive)\
xx |>\
  t() |>\
  cor()\
\
## plot both process errors against each other\
plot(ww[1,], ww[2,], pch = 16)\
\
## correlation of process errors (negative)\
ww |>\
  t() |>\
  cor()\
}