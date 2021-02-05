# Step 3 is DDJ and optional use of backtransform
# SRC 2020-01-30

rm(list = ls())
graphics.off()

source('DriftDiff_tau_weights_v2.r')

# save results
# Save post-processed DLM data
# X.raw is original series, X.z is z score, X.rawmean and X.rawsd are mean & sd for z score
# nl is DLM lags, delta is DLM delta, X.dlm is input to DLM, Tstep is time,
# local equilibria are X.eq and SD.eq, with ratio z.eq
# level (intercept) measures are level or stdlevel (for level/sd)
# Yyhat, B.ests, B.sd, errvar are DLM outputs
#save(nl,delta,X.raw,X.z,X.rawmean,X.rawsd,X.dlm,Tstep,X.eq,SD.eq,level,stdlevel,Z.eq,
#     Yyhat,B.ests,B.sd,errvar,file=Fname)
#
#load(file='Paul0811_nominalDLMresult.Rdata')
#fname=c('Paul0811_DDtauweightsv2.Rdata')
#load(file='Paul1315_nominalDLMresult.Rdata')
#fname=c('Paul1315_DDtauweightsv2.Rdata')
#
#load(file='Peter0811_nominalDLMresult.Rdata')
#fname=c('Peter0811_DDtauweightsv2.Rdata')
load(file='Peter1315_nominalDLMresult.Rdata')
fname=c('Peter1315_DDtauweightsv2.Rdata')

Xvar = stdlevel #[10:length(stdlevel)]
Tstep = Tstep[2:length(Tstep)]  # if Xvar = DLM output then trim Tstep

windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
plot(Tstep,Xvar,type='l',lwd=2,xlab='day',
     ylab='log10(phycocyanin)')
grid()

# Set up for DDJ_tau
title = c('Pigment variate')
nx = length(Xvar)
taus = c(1:5)
ntau = length(taus)
taumax = max(taus)

# Make dx matrix for further analysis
X0all = Xvar[1:(nx-taumax)]
nxs = length(X0all)
dxmat = matrix(0,nr=nxs,nc=5)
for(i in 1:ntau)  {
  tau = taus[i]
  X0 = Xvar[1:(nx-tau)]
  X1 = Xvar[(tau+1):nx]
  dx = X1-X0
  dxmat[,i] = t(dx[1:nxs])
}
 
# Make avec
bw <- 0.1*sd(X0all)  # try between 0.1*sd and 0.5*sd 
alow <- min(X0all)+bw
ahigh <- max(X0all)-bw
na <- 100
avec <- seq(alow,ahigh,length.out=na)

# Drift-Diff-Jump analysis ==================================================

# Bandi estimates

ParEst = Bandi_tau(taus,nxs,X0all,dxmat,bw,na,avec)

# smooth D1 and D2
D1 = ParEst[[1]]
D2 = ParEst[[2]]
sigma = sqrt(2*D2)

windows()
par(mfrow=c(3,1),mar=c(4,4.5,2,2)+0.1,cex.axis=2,cex.lab=2)
plot(avec,D1,type='l',lwd=2,col='blue',xlab='state',ylab='D1')
abline(h=0,lty=2)
plot(avec,sigma,type='l',lwd=2,col='red',xlab='state',ylab='sigma')
plot(avec,D2,type='l',lwd=2,col='red',xlab='state',ylab='D2')

# Save avec for downstream steps
avec.tau = avec

# Find equilibria
sdrift = sign(D1)
dsdrift = c(0,-diff(sdrift))
xeq = avec[which(!dsdrift == 0)]
ixeq = which(!dsdrift == 0)  # indices of the equilibria

print('equilibria and their indices',quote=F)
print(xeq,quote=F)
print(ixeq,quote=F)

# Save results using tau
# save line
save(Tstep,Xvar,avec,sigma,D1,D2,xeq,file=fname)

