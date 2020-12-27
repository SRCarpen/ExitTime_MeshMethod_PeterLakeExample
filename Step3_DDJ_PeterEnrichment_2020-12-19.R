# Step 3 is DDJ and optional use of backtransform
# NOTE: Step 2 is a backtransformation procedure that is not used here
# Copyright Stephen R. Carpenter 2020-01-30

rm(list = ls())
graphics.off()

source('DriftDiffJumpFunction.r')

# save results
# Save post-processed DLM data
# X.raw is original series, X.z is z score, X.rawmean and X.rawsd are mean & sd for z score
# nl is DLM lags, delta is DLM delta, X.dlm is input to DLM, Tstep is time,
# local equilibria are X.eq and SD.eq, with ratio z.eq
# level (intercept) measures are level or stdlevel (for level/sd)
# Yyhat, B.ests, B.sd, errvar are DLM outputs
#save(nl,delta,X.raw,X.z,X.rawmean,X.rawsd,X.dlm,Tstep,X.eq,SD.eq,level,stdlevel,Z.eq,
#     Yyhat,B.ests,B.sd,errvar,file=Fname)
load(file='Peter1315_DLMresult.Rdata')
fname=c('Peter1315_DDJresult.Rdata')

Xvar = stdlevel #[10:length(stdlevel)]
Tstep = Tstep[2:length(Tstep)]  # if Xvar = DLM output then trim Tstep

windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
plot(Tstep,Xvar,type='l',lwd=2,xlab='day',
     ylab='log10(phycocyanin)')
grid()

# Set up for DDJ
title = c('Pigment variate')
nx = length(Xvar)
#Tstep = c(1:nx)
X0 = Xvar[1:(nx-1)]
X1 = Xvar[2:nx]
dx = X1-X0
DT <- 1 #Tstep[2]-Tstep[1]
bw <- 0.2*sd(Xvar)  #nominal 0.5*sd
alow <- min(X0)+bw
ahigh <- max(X0)-bw
na <- 200
avec <- seq(alow,ahigh,length.out=na)

windows()
par(mfrow=c(2,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
plot(Tstep,Xvar,type='l',lwd=2,xlab='day',
     ylab='log10(phycocyanin)')
grid()
plot(Tstep[2:nx],dx,type='l',xlab='day',
     ylab='daily change')

# Drift-Diff-Jump analysis ==================================================

# Bandi estimates

ParEst = Bandi4d(X0,dx,(nx-1),DT,bw,na,avec)
#
# Bandi outputs
Drift.vec <- ParEst[[2]]
TotSig.vec <- ParEst[[3]]
DiffSig.vec <- ParEst[[4]]
JumpSig <- ParEst[[5]]
Lamda.vec <- ParEst[[6]]

# Display results
print(' ',quote=FALSE)
jsig = round(JumpSig,digits=5)
print(c('Jump sigma  ',round(JumpSig,digits=5)),quote=FALSE)
j.title = paste('Jump magnitude = ',jsig)

windows()
par(mfrow=c(2,2),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)

plot(avec,Drift.vec,type='l',lwd=2,col='black',xlab='State Value',ylab='drift',main=title)
abline(h=0,col='blue')
grid()

plot(avec,TotSig.vec,type='l',lwd=2,col='black',xlab='State Value',
     ylab='total sigma') #,main='check sensitivity to x-range') 
grid()

plot(avec,Lamda.vec,type='l',lwd=2,col='black',xlab='State Value',
     ylab='jump frequency',main=j.title) 
grid()

plot(avec,DiffSig.vec,type='l',lwd=2,col='black',#ylim=yrange,
     xlab='State Value',ylab='diffusion sigma') 
grid()

# Summary graphic
jump.component = jsig*Lamda.vec

windows(width=8,height=8)
par(mfrow=c(2,1),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
#xrange=c(0.5,2.3)
#yrange=c(-0.02,0.02) #range(Drift.vec)
plot(avec,Drift.vec,type='l',lwd=2,col='black',#xlim=xrange,ylim=yrange,
     xlab='Pigment',
     ylab='drift')
abline(h=0,col='blue')
grid()
#yrange = range(c(TotSig.vec,DiffSig.vec,jump.component),na.rm=T)
plot(avec,TotSig.vec,type='l',lwd=2,col='black',#ylim=yrange,#xlim=xrange,
     xlab='Pigment',ylab='SD of environmental noise')
grid()

# Arts and crafts
potential = -1*(cumsum(Drift.vec))

windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 5, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
#xrange=c(0.5,2.5)
#yrange=c(-0.5,1.5)
plot(avec,potential,type='l',lwd=2,col='blue',#xlim=xrange,ylim=yrange,
     xlab='Pigment',
     ylab='Negative Potential')
grid()


# Find equilibria
sdrift = sign(Drift.vec)
dsdrift = c(0,-diff(sdrift))
xeq = avec[which(!dsdrift == 0)]
ixeq = which(!dsdrift == 0)  # indices of the equilibria

print('Equilibria and their time-series indices',quote=F)
print(xeq)
print(ixeq)

# Save drift and diffusion curves and backtransform parameters
# lr1 is a linear regression, lo1 is a quadratic loess
# nlpar is the nonlinear fitted parameters L, k, x0, and C where 
# deltax = Xmeans[,1] - x0
# yhat = ( L/(1 + exp(-k*deltax)) ) - C
#save(avec,Drift.vec,TotSig.vec,xeq,lr1,lo1,nlpar,file=fname)
#
# save line if there is no BT
save(Tstep,Xvar,avec,Drift.vec,TotSig.vec,xeq,potential,file=fname)
