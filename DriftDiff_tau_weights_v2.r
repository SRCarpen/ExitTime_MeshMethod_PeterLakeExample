# Function to compute nonparametric drift, diffusion estimators
#
#    Based the paper "Early warnings of unknown nonlinear shifts: a nonparametric
#    approach" by S.R. Carpenter and W.A. Brock
#    
#    with tau correction from Rinn et al., 'The Langevin Approach: An R package for
#    modeling Markov processes', 2016, arXiv:1603.02036v1
#    See equations 11-13
# Programmer: Stephen R. Carpenter, January 2021
#
# Inputs:
#
#  x0 is the predictor vector, which can be any covariate measured at the same
#     interval as the response variable. In the paper x0 is x-subscript-i*delta.
#  dx is the response variable vector. In the paper it is
#     x-subscript-(i+1)*delta minus x-subscript-i*delta
#  nx is length of the vector dx
#  DT is time step
#  bw is the bandwidth. See Hardle 1990 and Johannes 2004 (cited in the paper)
#     for advice about choice of bandwidth. In our experience, more stable,
#     accurate and precise estimates are obtained by using bandwidths larger
#     than the "optimal" bandwidths computed using 'density' {kernsmooth} and
#     similar algorithms.
#  na is the length of the avec vector
#  avec is a vector of values of the predictor variable for which nonparametric
#      estimates will be computed
#
# Output: The function returns a list of 6 variables, as follows:
#
# avec: same as the input avec
# mu.x: vector of nonparametric drift estimates for each element of avec
# sigma.x: vector of nonparametric total sigma (not sigma^2) estimates for
#     each element of avec
# sigma.diff: vector of nonparametric diffusion sigma (not sigma^2) estimates for
#     each element of avec
# sigma.z: a scalar; nonparametric jump sigma (not sigma^2)
# lamda.z: vector of nonparametric jump frequency (or jump intensity) lamda
#     estimates for each element of avec

Bandi_tau <- function(taus,nx,X0all,dxmat,bw,na,avec)  { 
#
#nx = length(dxmat[,1])
ntau = length(taus)
x0 = X0all
SF <- 1/(bw*sqrt(2*pi))  # scale factor for kernel calculation
# square dx matrix (each column is a tau)
dx2 = dxmat*dxmat
dx4 = dx2*dx2

# Compute matrix of kernel values
Kmat <- matrix(0,nrow=na,ncol=nx)
for(i in 1:(nx)) {  # loop over columns (x0 values)
  Kmat[,i] <- SF*exp(-0.5*(x0[i]-avec)*(x0[i]-avec)/(bw*bw))
}
#
# Moments by tau
M1mat = matrix(0,nr=na,nc=ntau)
M2mat = matrix(0,nr=na,nc=ntau)
M4mat = matrix(0,nr=na,nc=ntau)
#
# Kernel-weighted moments versus a
for(i in 1:na) {  # loop over rows (a values)
  Ksum <- sum(Kmat[i,])  # sum of weights
  for(j in 1:ntau) {   # loop over columns(taus)
    # original version divides by taus[j]
    #M1mat[i,j] = (1/taus[j])*sum(Kmat[i,]*dxmat[,j])/Ksum
    #M2mat[i,j] = (1/taus[j])*sum(Kmat[i,]*dx2[,j])/Ksum
    #M4mat[i,j] = (1/taus[j])*sum(Kmat[i,]*dx4[,j])/Ksum
    # This version does not divide by taus[j]; per time comes in regression slope
    M1mat[i,j] = sum(Kmat[i,]*dxmat[,j])/Ksum
    M2mat[i,j] = sum(Kmat[i,]*dx2[,j])/Ksum
    M4mat[i,j] = sum(Kmat[i,]*dx4[,j])/Ksum
  } 
}
D1 = rep(0,na)
D2 = rep(0,na)
# Compute regressions by row (na)
for(i in 1:na) {
  yvec = M1mat[i,]
  xvec = as.vector(taus)
  # regression for D1
  inv.erD1 = nx/(M2mat[i,] -(M1mat[i,]^2) )
  WD1 = diag(inv.erD1)
  D1[i] = ( 1/(t(xvec)%*%WD1%*%xvec) )*( t(xvec)%*%WD1%*%yvec )
  # regression for D2
  inv.erD2 = nx/(M4mat[i,] -(M2mat[i,]^2) )
  WD2 = diag(inv.erD2)
  yvec = M2mat[i,] - (D1[i]*taus)*(D1[i]*taus)
  xvec = as.vector(2*taus)
  D2[i] = 0.5*( 1/(t(xvec)%*%WD2%*%xvec) )*(t(xvec)%*%WD2%*%yvec)
}

outlist = list(D1,D2)
return(outlist) 

} 

