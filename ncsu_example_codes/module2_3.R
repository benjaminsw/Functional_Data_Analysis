# ref: https://www4.stat.ncsu.edu/~staicu/FDAtutorial/module2-3.html
# 1 ###########################################################################
library(fda)
data(CanadianWeather)
attach(CanadianWeather)

y.precip = dailyAv[,,2]
l = which(place=="Vancouver") 
t.day = 1:365  # define function arguments
y=y.precip[,l] # define response of interest

plot(t.day, y, pch = 1, cex = 0.5, col='royalblue2',
     xlab="day", ylab="Mean Precipitation",
     main=place[l])
# 2 ###########################################################################
# install.packages("mgcv")
library(mgcv)
fit <- gam(y~s(t.day, k = 10, bs = "cr"), method="REML")
# names(fit)
yhat <- fit$fitted

plot(t.day, y, type="n", 
     xlab="day", ylab="mean precipitation", main=place[l])
points(t.day, y, pch=1, col="blue", cex=.5)
lines(t.day,  yhat, lwd=2, col="black")
# 3 ###########################################################################
Xtrue = yhat   # "true" underlying curve
Eps = y - yhat      # residual
m=length(Eps)

hist(Eps)
# 4 ###########################################################################
set.seed(1)
Eps.star = sample(Eps, m, replace=T)
Y.star = Xtrue + Eps.star 
plot(t.day, Y.star, pch = 1, cex = 0.5, col='red',
     xlab="day", ylab="Y.star",
     main="Simulated Dataset")
points(t.day, y, pch=1, col="light grey", cex=.5)       # original data
lines(t.day,  yhat, lwd=2, col="black")   # truth
# 5 ###########################################################################
K.vec = 2*c(2:8)+1;    # we consider 7 different values of k.
# K.vec

t.day <- 1:365
fbasis=create.fourier.basis(rangeval = c(1, 365), nbasis=max(K.vec), period=365)
bvals = eval.basis(t.day,fbasis)

Xfit = array(0, c(m, length(K.vec)))
index=0
for (K in K.vec){
  index=index+1
  Xbasis = bvals[, 1:K]
  lm.fit = lm(Y.star~0+Xbasis)
  Xfit[,index] = as.vector(lm.fit$fitted.values)
}

plot(t.day, Y.star, pch = 1, cex = 0.5, col='grey',
     xlab="day", ylab="Y.star",
     main="Simulated Dataset")
lines(t.day,  yhat, lwd=2, col="black")   # truth
lines(t.day, Xfit[,1], lwd=2, col="red")
lines(t.day, Xfit[,7], lwd=2, col="blue")
# 6 ###########################################################################
K.vec = 2*c(2:8)+1; 

## pretend 
Xtrue = yhat
Eps = y - Xtrue ; m=length(Eps)
B=100 #B=10000

Xfit = array(0, c(B, m, length(K.vec)))

fbasis=create.fourier.basis(rangeval = c(1,365), nbasis=max(K.vec), period=365)
bvals = eval.basis(t.day,fbasis)

set.seed0=1234
for(b in 1:B){#b=1
  set.seed(set.seed0+b)
  Eps.star = sample(Eps, m, replace=T)
  Y.star = Xtrue + Eps.star 
  
  # fit using Fourier basis and K basis functions
  index=0
  for (k in K.vec){
    index=index+1
    Xbasis = bvals[, 1:k]
    lm.fit = lm(Y.star~0+Xbasis)
    Xfit[b,,index] = as.vector(lm.fit$fitted.values)
  }
}

Mean.Est = apply(Xfit, c(2,3), mean)
Mean.Est2 = apply(Xfit, c(2,3), function(x) mean(x^2))

Bias = apply(Mean.Est, 2, function(x) Xtrue-x)
Var = Mean.Est2 - (Mean.Est)^2
Mse= Bias^2+Var

Mean_Bias2_L2 = apply(Bias^2, 2, mean) 
Var_L2 = apply(Var, 2, mean) 
MSE_L2 = apply(Mse, 2, mean)

plot(K.vec, Mean_Bias2_L2, type="n", cex=2, 
     xlab="Number of basis functions", ylab="Total squared error", 
     ylim=range(cbind(Mean_Bias2_L2, MSE_L2)))

points(K.vec, Mean_Bias2_L2, type='b', col="darkgreen", lwd=2)
points(K.vec, Var_L2, type='b', col="blue", lwd=2)
points(K.vec, MSE_L2, type='b', col="orange", lwd=2)
# 7 ###########################################################################
y.precip=dailyAv[,,2]
l = which(place=="Vancouver") 
t.day = 1:365  
y=y.precip[,l]
K.vec = 2*c(2:8)+1; 
CVfit = matrix(0,  nrow=m, ncol=length(K.vec))
for(j in 1:m){
  
  Y.star = y[-j]
  
  # fit using Fourier basis and K basis functions
  index=0
  for (K in K.vec){
    index=index+1
    Xbasis=bvals[, 1:K];
    Xbasis.j =  Xbasis[-j, ]; 
    lm.fit = lm(Y.star~0+Xbasis.j); Xbasis.coeff = lm.fit$coefficients
    y.fit = Xbasis%*%Xbasis.coeff
    CVfit[j,index] = (y[j] - y.fit[j])^2
  }
}

CV_L2 = apply(CVfit, 2, sum)
plot(K.vec, CV_L2, type="n",
     xlab="Number of basis functions", ylab="Total cross-validation error")
points(K.vec, CV_L2, type='b', col="royalblue2", lwd=2)
title(paste0("K = ", K.vec[which(CV_L2==min(CV_L2))], " with the smallest CV score!"))
