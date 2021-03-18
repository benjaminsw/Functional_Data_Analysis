#' ---
#' title: "Lab 5: Functional Response Models "
#' author: "Surajit Ray"
#' output:
#'    pdf_document:
#'      highlight: zenburn
#' ---

#'
#' ### Loading R library
###' Firstly, let's load up the library and the data used in Lab 1.
par(ask=FALSE)
library('fda')
load("all.RData")
data(CanadianWeather)

temp = CanadianWeather$dailyAv[,,1]
precip = CanadianWeather$dailyAv[,,2]

daytime = (1:365)-0.5
day5 = seq(0,365,5)

dayrng = c(0,365)





#'###############################################################################
#' #                      FUNCTIONAL RESPONSE MODELS                         ####
#'###############################################################################

#'  We could predict total precipitation from temperature. How about the
#'  annual precipitation profile?

#'  We can also look at constant predictors -- see the weather demo for an ANOVA
#'  between different weather regions.

#'  First we'll create a smooth of the log precipitation

prec.fdPar = fdPar(fbasis,harmLfd,1e6)
precSmooth = smooth.basis(daytime,log(precip+0.5),prec.fdPar)
precfd = precSmooth$fd

#'  We can retain xlist from the scalar response model.

#'
#' ## 1. fdPar objects and estimation

bwtlist2 = list(len=2)

#'  The intercept is now a functional parameter as well as beta 1.   Since this
#'  is an identifiable model without smoothing, we'll set the smoothing parameter
#'  very low.

beta.fdPar2 = fdPar(fbasis,harmLfd,1e-5)

bwtlist2[[1]] = beta.fdPar2
bwtlist2[[2]] = beta.fdPar2

#'  We can also call fRegress with this

prec.conc = fRegress(precfd,xlist,bwtlist2)

#'  Let's have a look at what we've got

par(mfrow=c(2,1))
plot(prec.conc$betaestlist[[1]])
plot(prec.conc$betaestlist[[2]])

#'  We can also look at a comparison between predicted and observed

yhatfd  = prec.conc$yhatfdobj$fd  # Fitted smooths.

#plot(yhatdf)
plot(yhat)
plot(precfd)

#'  And compare observed with residuals

plot(precfd-yhatfd)
plot(precfd)

#'  Doesn't look like we really got very much.

#'
#' ## 2. Confidence Intervals

#'  In order to obtain confidence intervals, we can include the smoothing of
#'  precipitation as a source of error. In order to do this, we need two things

y2cmap = precSmooth$y2cMap

#'  This is the matrix that goes from the observations to the coefficients.

#'  We now need a covariance matrix for the errors of the original observed
#'  precipitation from the functional linear model

Errmat = precip - eval.fd(daytime,yhatfd)

SigmaE2 = cov(t(Errmat))

#'  We can now run fRegress.stderr

conc.std = fRegress.stderr(prec.conc,y2cmap,SigmaE2)

#'  And plot the results

plotbeta(prec.conc$betaestlist,conc.std$betastderrlist)

#'  There really doesn't appear to be much going on.

#'
#'    * **EXERCISE**: try predicting precipitation instead of log precipitation -- does
#'  this make a difference?

#'
#'    * **EXERCISE**: what diagnostics could be done to check goodness of fit? Try
#'  plotting residuals. Try plotting residuals against predicted values
#'  (this should give you a series of lines, you'll need to evaluate and use
#'  'matplot').    How could you check for dependence of precipitation on
#'  non-concurrent times?

#'
#' ## 3. Permutation Tests and Cross Validation

#'  The next two functions can take a very long while to run

#'  Permutation test for fRegress

par(mfrow=c(1,1),ask=FALSE)
Fresult = Fperm.fd(precfd,xlist,bwtlist2)

#'  Here the dotted line gives the 95th percentile of the permutation distribution
#'  at each time t, the dashed line gives the 95th percentile of the permutation
#'  distribution of the maximum F value, and the solid line gives the observed
#'  F value at each time t.

#'
#' ### Cross validated integrated squared error.

#'  This is a particularly long undertaking since the cross-validation is
#'  done manually. There are some matrix identities that can speed this up
#'  (still to be implemented).

#'  In this case we will only look at the same lambda for each coefficient
#'  function

lambdas = 10^(c(-5,-2,1,3))
SISEs = rep(0,length(lambdas))

ilam=1
# for(ilam  in 1:length(lambdas)){
beta.fdPari = fdPar(fbasis,harmLfd,lambdas[ilam])   # Update lambda
bwtlisti = list(beta.fdPari,beta.fdPari)

CVres = fRegress.CV(precfd,xlist,bwtlisti)

SISEs[ilam] = CVres$SSE.CV

print(c(ilam,SISEs[ilam]))      # Just so we know where we are.
# }

plot(lambdas,SISEs,type='b',log='x')

#'
#' ### Further data sets to try playing with:

#'   **gait** -- provides the gait data demo'd in class. Predict knee angle from
#'   hip angle, or try predicting a derivative. The data are in 'gait'

#'   Try predicting temperature from region in the weather data. You can find
#'   a region indicator in "CanadianWeather$Region".  BUT, you will need to set
#'   up an appropriate design matrix.
#'
#'   This is a functional version of an ANOVA.