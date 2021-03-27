# ref: https://www4.stat.ncsu.edu/~staicu/FDAtutorial/module2-4.html
library(fda)
data(CanadianWeather)
attach(CanadianWeather)

y.precip=dailyAv[,,2]
l = which(place=="Vancouver") 
t.day = 1:365  
y=y.precip[,l]

# define domain, #knots, and #order to construct b-spline basis

ybasis  <- create.bspline.basis(rangeval = c(1,365), nbasis = 365, norder=4)

bvals = eval.basis(t.day, ybasis)
Xbasis =bvals; 
lm.fit = lm(y ~ 0 + Xbasis)   
y.fit = lm.fit$fitted.values

plot(t.day, y, type="n",lwd=4, col="black",
     xlab="day", ylab="precipitation", 
     main=paste(365, "Fourier fns"), cex=1)
points(t.day, y, pch=1, cex=.5, col="blue", lwd=1)
lines(t.day, lm.fit$fitted.values, lwd=1, col="red")
# 2 ###########################################################################
lambda <- 10^4

# int2Lfd(m)  : use this to define the m-th order derivative penalty term
# fdPar() : defines functional parameters; in this case the 2nd order derivative penalty term and the smoothing parameter.

# ybasis  <- create.bspline.basis(rangeval = c(1,365), nbasis = 365, norder=4)
tD2fdPar = fdPar(ybasis, Lfdobj=int2Lfd(2), lambda=lambda)

# smooth.basis() : smoothes the data using the roughness penalty and smoothing parameter specified in 'tD2fdPar' 
tyfd = smooth.basis(t.day,y,tD2fdPar) 

#names(tyfd)
#[1] "fd"      "df"      "gcv"     "beta"    "SSE"     "penmat"  "y2cMap"     
#    "argvals" "y"    

# fd   a functional data object containing a smooth of the data.
# df     a degrees of freedom measure of the smooth
# gcv  the value of the generalized cross-validation or GCV criterion. 
# beta the regression coefficients associated with covariate variables. 
# SSE    the error sums of squares. 
# penmat:the penalty matrix.
# y2cMap     the matrix mapping the data to the coefficients: 
#          (Phi^T Phi + R)^(-1) \Phi^T

main.label = paste("Vancouver (lambda =", round(lambda,2), ")", sep="")
plot(t.day, y, type="n", ylim=range(y), 
     ylab="Precipitation", xlab="day", main=main.label)
points(t.day, y, pch=1, cex=.5, col="blue", lwd=1)
lines(tyfd$fd,col="red",lwd=4)
# 3 ############################################################################
logl=seq(-5, 12, len=71)  
range(exp(logl))
gcv = rep(0,71)

for(i in c(1:length(logl))){
  lambda=exp(logl[i])
  
  tD2fdPar = fdPar(ybasis,Lfdobj=int2Lfd(2),lambda=lambda)
  tyfd = smooth.basis(t.day,y,tD2fdPar)
  
  gcv[i] = tyfd$gcv
}

# PLOT GCV of FIT versus log lambda
plot(logl,gcv[1:71],type='l',cex.lab=1.5, lwd=4, 
     xlab='log lambda',ylab='GCV', main="GCV(log.lambda)")


