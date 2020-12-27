# Trial of general ET programs - step 1 is DLM
# Step 2 is a backtransformation procedure that is not used here
# Copyright Stephen R. Carpenter 2020-01-27

rm(list = ls())
graphics.off()

source('ODLMAR_NoBoot_2018-10-20.R')

# Choose data from Load_All_MARSSed_squeal2_data_V2.R
# Each data frame is Squeal II, 1 lake, all years
# columns are:  select=c('Year','Lake','DoY','BGA_HYLB','BGA_logged_HYLB')
load(file='Peter_Enrichment_2013-15.Rdata')

# show variates in PeteBGA
print(Pete_BGA[1,])

# Example: choose variates, and title-------------------------------------
lakesub = subset(Pete_BGA,select=c('Year','DoY','BGA_logged_HYLB'))
title = c('log phyco, Peter, 2013-2015')
Fname = c('Peter1315_DLMresult.Rdata')

# set time index and time series for analysis
pig = lakesub$BGA_logged_HYLB
doy = lakesub$DoY
yr = lakesub$Year
mindoy = min(doy)
maxdoy = max(doy)
Tscore = yr + ( (doy - mindoy)/(maxdoy-mindoy+1) )

windows(width=12,height=6)
par(mfrow=c(1,1),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Tscore,pig,type='l',
     xlab='Day of Year',ylab='log10 BGA',ylim=c(1.5,4.5),main=title)
grid()

# Start DLM
# z score function
zscore=function(x) {
  xbar = mean(x,na.rm=T)
  xsd = sd(x,na.rm=T)
  z = (x - xbar)/xsd
  return(z)
}

# optionally convert to logs and take zscore
Tstep = Tscore
X.raw = pig
X.z = zscore(X.raw)
X.rawmean = mean(X.raw,na.rm=T)
X.rawsd = sd(X.raw,na.rm=T)
print(c('mean & s.d. of X.raw ',X.rawmean,X.rawsd),quote=F)

# Select input to dlm
X.dlm = X.z

windows(width=12,height=6)
plot(Tstep,X.dlm,type='l',col='forestgreen',xlab='DoY index',ylab='x.dlm',
     main='time series for DLM')
grid()

# Set up DLM
nobs = length(X.dlm)
nl = 1 # number of lags
delta = 0.95 # 0<delta<1; see advice in functions

# Run DLM
ODL.out = ODLMAR(nl,delta,X.dlm,Tstep,title)

# Output matrices are stored sideways, like MARSS
Yyhat = ODL.out[[1]]
EigenVals = ODL.out[[2]]
B.ests = ODL.out[[3]]  
B.sd = ODL.out[[4]]
errvar = ODL.out[[5]] # updated error variance

# Post process DLM -----------------------------------------------

# Calculate moving equilibrium
X.eq = B.ests[1,]/(1 - B.ests[2,])
# Calculate its variance
SDterm1 = X.eq*X.eq
SDterm2 = (B.sd[1,]*B.sd[1,] + errvar)/(B.ests[1,]*B.ests[1,])
SDterm3 = (B.sd[2,]*B.sd[2,])/((1 - B.ests[2,])*(1 - B.ests[2,]))
SD.eq = sqrt(SDterm1*(SDterm2 + SDterm3))
# Z score
Z.eq = X.eq/SD.eq

# Time steps start at 2
Nstep = length(Tstep)

# Plot components of steady-state estimate
windows(width=6,height=12)
par(mfrow=c(3,1),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Tstep[2:Nstep],X.eq,type='l',col='blue',ylim=c(-10,10),
     ylab='Steady State',xlab='DoY',
     main='Local Steady-State estimate, sd, and ratio')
grid()
plot(Tstep[2:Nstep],SD.eq,type='l',col='red',ylim=c(0,10),
     ylab='S.D.',xlab='DoY')
grid()
plot(Tstep[2:Nstep],Z.eq,type='l',col='purple',
     ylab='Z score',xlab='DoY')
grid()

# Calculate level estimates
level = B.ests[1,]
stdlevel = B.ests[1,]/B.sd[1,]

# Plot components of level estimate
windows(width=12,height=9)
par(mfrow=c(2,1),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Tstep[2:Nstep],level,type='l',col='blue',#ylim=c(-10,10),
     ylab='level',xlab='DoY',
     main='Level and Std level estimate')
grid()
plot(Tstep[2:Nstep],stdlevel,type='l',col='red',#ylim=c(0,10),
     ylab='Std Level',xlab='DoY')
grid()

# Save post-processed DLM data
# X.raw is original series, X.z is z score, X.rawmean and X.rawsd are mean & sd for z score
# nl is DLM lags, delta is DLM delta, X.dlm is input to DLM, Tstep is time,
# local equilibria are X.eq and SD.eq, with ratio z.eq
# level (intercept) measures are level or stdlevel (for level/sd)
# Yyhat, B.ests, B.sd, errvar are DLM outputs
save(nl,delta,X.raw,X.z,X.rawmean,X.rawsd,X.dlm,Tstep,X.eq,SD.eq,level,stdlevel,Z.eq,
     Yyhat,B.ests,B.sd,errvar,file=Fname)

