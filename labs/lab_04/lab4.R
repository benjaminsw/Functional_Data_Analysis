
#' ---
#' title: "Lab 4: Functional Linear Models "
#' author: "Surajit Ray"
#' output:
#'    pdf_document:
#'      highlight: zenburn
#' ---

#' '

#'  This lab will focus on functional linear models using the fRegress function.
#'  We will again make use of the Canadian Weather Data. We'll briefly re-do
#'  the smooth here using Fourier bases:

#'
#'  ## Loading R library
#' Firstly, let's load up the library and the data used in Lab 1.
par(ask=FALSE)
library('fda')
load("all.RData")
data(CanadianWeather)

temp = CanadianWeather$dailyAv[,,1]
precip = CanadianWeather$dailyAv[,,2]

daytime = (1:365)-0.5
day5 = seq(0,365,5)

dayrng = c(0,365)

#'  Basis and parameter objects

fbasis =  create.fourier.basis(dayrng,65)
harmLfd = vec2Lfd(c(0,(2*pi/365)^2,0),rangeval=dayrng)
temp.fdPar = fdPar(fbasis,harmLfd,1e-2)

#'  And smooth

tempSmooth = smooth.basis(daytime,temp,temp.fdPar)

#'  And extract the functional data objects

tempfd = tempSmooth$fd



#'###############################################################################
#' #                         SCALAR RESPONSE MODELS                          ###
#'###############################################################################


#'  We examine predicting a scalar response from a functional covariate.

#'  In this case, we'll predict log annual precipitation from the temperature
#'  profile.  This is one of the most over-used examples in FDA; so get it out
#'  of your system here.

#'
#' ### 1. Setup Data

#'  First we'll obtain log annual precipitation

annualprec = log10( apply(precip,2,mean))

#'  Now we need to set up a list of covariates.

xlist = list(len=2)

#'  First co-variate is just the intercept: a vector of ones

xlist[[1]] = rep(1,35)

#'  Second covariate is temperature

xlist[[2]] = tempfd


#' ### 2. fdPar objects for coeffients

#'  We also need a list of functional parameter objects to define the coefficient
#'  functions.

bwtlist = list(len=2)

#'  First is a constant basis and add it to the fdPar object

cbasis = create.constant.basis(dayrng)
bwtlist[[1]] = fdPar(cbasis)


#'  Now we need the coefficient of temperature, we'll use the same basis for it
#'  and add it as the second element in the list.

beta.fdPar = fdPar(fbasis,harmLfd,10^12.5)
bwtlist[[2]] = beta.fdPar

#' ### 3. fRegress and outcome
annualprec=as.numeric(annualprec)
prec.model = fRegress(annualprec,xlist,bwtlist)

#'  Let's look at the outcome

names(prec.model)

#'  We can see the intercept as

prec.model$betaestlist[[1]]$fd$coef

#'  We can also plot the estimated regression function

plot(prec.model$betaestlist[[2]])


#' ### 4. Cross-validation

#'  We should look at selecting lambda. We'll do this with OCV

lambdas = 10^(seq(5,15,0.5))

ocvs = rep(0,length(lambdas))

for(ilam in 1:length(lambdas)){
  bwtlisti = bwtlist             # define temporary beta.fdPar and bwtlist
  beta.fdPari = beta.fdPar
  
  beta.fdPari$lambda = lambdas[ilam]   # update lambda
  bwtlisti[[2]] = beta.fdPari
  
  prec.modeli = fRegress(y=annualprec,xlist,bwtlisti)
  
  ocvs[ilam] = prec.modeli$OCV        # record ocv
}

plot(lambdas,ocvs,type='b',log='x')

#'  It looks like our original choice of 12.5 was about right


#' ### 4. Statistics, Standard Errors and Tests

#'  Degrees of freedom

prec.model$df

#'  We'll plot y-by-yhat (in this case yhatfdobj is just a scalar despite
#'  its name).

yhat = prec.model$yhatfdobj

plot(yhat,annualprec)
abline(c(0,1))

#'  And we can caculate a residual variance

sigma = sum( (annualprec-yhat)^2 )/(35-prec.model$df)
sigma


#'  To obtain standard error estiamtes we can now call

sigmaE = sigma*diag(35)
prec.stderr = fRegress.stderr(prec.model,NULL,sigmaE)

#'  And we can obtain plots for beta from the estimated and standarderror

betahat = prec.model$betaestlist[[2]]
betastd = prec.stderr$betastderrlist[[2]]

plotbeta(betahat,betastd)

#'
#'    * **EXERCISE**: obtain a confidence interval for the intercept.
alphahat = prec.model$betaestlist[[1]]$fd$coef
alphastd = prec.stderr$betastderrlist[[1]]$coefs

conf1=alphahat-2*alphastd
conf2=alphahat+2*alphastd
c(conf1,conf2)

#'    * **THOUGHT EXERCISE**: $\beta_1(t)$ integrates to zero; why?


#'  Finally, we can run a permutation test on this model; this will take some
#'  time.

par(mfrow=c(1,1),ask=FALSE)
Fresult = Fperm.fd(annualprec,xlist,bwtlist)

#'  The histogram gives the permutation distribution. Dashed lines are the 95
#'  quantile and the solid line is the observed F value.

#'
#'    * **EXERCISE**: plot residuals against predicted values. This may not indicate poor fit if there is a non-linear relationship only with one part of the covariate function. Try a 3-dimensional plot putting time on the 'x' axis,
#' the covariate value on the 'y' axis and residuals on the 'z' axis. The #
#' library 'rgl' is particularly good for this.

res = annualprec-yhat
plot(res,yhat)
abline(c(0,0))

library(rgl)
plot3d(daytime, yhat, res)


#'    * **EXERCISE**: how sensitive are these results to the amount of smoothing of the
#'  temperature process? Try lambda at 1e-2 and 1e2.

beta.fdPar$lambda = 1e-2 # update lambda
bwtlist[[2]] = beta.fdPar

par(mfrow=c(1,1),ask=FALSE)
Fresult = Fperm.fd(annualprec,xlist,bwtlist)


beta.fdPar$lambda = 1e2  # update lambda
bwtlist[[2]] = beta.fdPar

par(mfrow=c(1,1),ask=FALSE)
Fresult = Fperm.fd(annualprec,xlist,bwtlist)


#'
#' #              FUNCTIONAL PRINCIPAL COMPONENTS REGRESSION                 ##

#'  Here we will continue the problem above, but we will tackle it from the
#'  perspective of functional Principal Components Regression.

#'
#' ## 1. Obtaining an Estimate

#'  First we need to re-obtain fPCAs

tempPCA = pca.fd(tempfd,nharm=6)

#'  We'll continue to use the first three PCAs and extract their scores into
#'  a predictor matrix

Xmat = tempPCA$scores[,1:3]

#'  Now perform linear regression

prec.lm = lm(annualprec~Xmat)

#'  We can already obtain some summary statistics

summary(prec.lm)

#'  and try the same trick

plot(prec.lm$fitted,annualprec)
abline(c(0,1))

#'  Now we want to reconstruct beta. First we'll make the code easy to read by
#'  obtaining the PCs and the coefficients

Mcoefs = prec.lm$coef
PCs = tempPCA$harmonics[1:3]

#'  and put them together to get beta; remember, we still leave out the intercept.

beta = Mcoefs[2]*PCs[1] + Mcoefs[3]*PCs[2] + Mcoefs[4]*PCs[3]


#'  Now we can plot the result

plot(beta)

#'    * **EXERCISE**: this is pretty rough -- try cross-validating the amount of
#' ## smoothing in the PCA analysis based on its ability to predict log
#' ## annual precipitation.
#' 


lambdas = 10^(seq(2,6,0.5))

Rsqs = rep(0,length(lambdas))

curv.Lfd = int2Lfd(2)
bbasis = create.bspline.basis(dayrng,nbasis=70,norder=4) 

for(ilam in 1:length(lambdas)){
  pca.fdPari = fdPar(bbasis,curv.Lfd,lambdas[ilam])
  
  tempPCAsmoothi = pca.fd(tempfd,nharm=6,harmfdPar=pca.fdPari)
  
  Xmat.smoothi = tempPCAsmoothi$scores[,1:3]
  
  prec.lm.smoothi = lm(annualprec~Xmat.smoothi)
  
  Rsqs[ilam] = summary(prec.lm.smoothi)$r.squared  
  #Note that Rsquared is not the only answer for assessing prediction ability. You can choose whatever is reasonable. 
}

plot(lambdas,Rsqs,type='b',log='x')

bestlam <- lambdas[which.max(Rsqs)] 


pca.fdPari = fdPar(bbasis,curv.Lfd,bestlam)

tempPCAsmoothi = pca.fd(tempfd,nharm=6,harmfdPar=pca.fdPari)

Xmat.smoothi = tempPCAsmoothi$scores[,1:3]

prec.lm.smoothi = lm(annualprec~Xmat.smoothi)


plot(prec.lm.smoothi$fitted,annualprec)
abline(c(0,1))

Mcoefs.smooth = prec.lm.smoothi$coef
PCs.smooth = tempPCAsmoothi$harmonics[1:3]

beta.smooth = Mcoefs.smooth[2]*PCs.smooth[1] + Mcoefs.smooth[3]*PCs.smooth[2] + Mcoefs.smooth[4]*PCs.smooth[3]

plot(beta.smooth)




#'
#' ##2. Standard Errors

#'  We will do this manually. First we will use the usual variance

varbeta = sigma*solve(t(Xmat)%*%Xmat)

#'  Now we'll obtain the values of the PCs at a set of plotting points

PCvals = eval.fd(day5,PCs)

#'  The covariance for beta is then

PCbetacov = PCvals%*%varbeta%*%t(PCvals)

#'  We can actually have a look at the whole covariance surface

contour(day5,day5,PCbetacov)

#'  But largely we just want to extract the diagonal and then plot it.

#'  First we'll get values for beta

PCbetavals = eval.fd(day5,beta)

#'  Then standard errors

PCbetastd = sqrt(diag(PCbetacov))

#'  And we can form a plot

plot(day5,PCbetavals,type='l',lwd=2,ylim=c(-6e-4,1e-3))
lines(day5,PCbetavals+2*PCbetastd,type='l',lwd=2,lty=2)
lines(day5,PCbetavals-2*PCbetastd,type='l',lwd=2,lty=2)
abline(h=0)

#'  This is actually pretty similar to the fRegress version and will improve
#'  with smoothing the PCA.
