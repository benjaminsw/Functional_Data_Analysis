
#' ---
#' title: "Lab 3: PCA "
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



################################################################################
#' #                   FUNCTIONAL PRINCIPAL COMPONENTS                       ####
################################################################################

#'
#' ##1. pca.fd

#'  We can conduct a fPCA through

curv.Lfd = int2Lfd(2)
bbasis = create.bspline.basis(dayrng,nbasis=70,norder=4) 
tempSmooth = smooth.basis(daytime,temp,bbasis)
tempfd = tempSmooth$fd
tempPCA = pca.fd(tempfd,nharm=6)

#'  Here we can look at proportion of variance explained:

plot(tempPCA$varprop,type='b')

#'  Looks like we could have stopped at 3.

## Looking at the principal components:

plot(tempPCA$harmonics[1:3])

#'  1 Looks like over-all temperature.
#'  2 Looks like Summer vs Winter
#'  3 Is Spring vs Fall.

## But let's plot these
par(mfrow=c(1,3))
plot(tempPCA,harm=1:3)
# Which gives a much better picture of what's going on.

#' ##2. PCA and Smoothing

#'  The PCs above are fairly rough, we can also add a smoothing penalty (see
#'  the special topics slides).

pca.fdPar = fdPar(bbasis,curv.Lfd,1e4)

tempPCAsmooth = pca.fd(tempfd,nharm=6,harmfdPar=pca.fdPar)

#'  Let's plot the PCs

plot(tempPCAsmooth$harmonics[1:3])

#'  We can see that these are considerably smoother but still have pretty much
#'  the same interpretation.

plot(tempPCAsmooth)

#'
#' ##3. PCA and Reconstructing Data

#'  We can ask how well the PCs reconstruct the observed functions. We'll focus
#'  on the first observation and reconstructing using PC score.

#'  Let's extract the scores and PCs just to make the code easier to read
mtempfd = mean.fd(tempfd)
scores = tempPCAsmooth$scores
PCs = tempPCAsmooth$harmonics

?pca.fd
plot(scores[,1:2])
#'  Firstly, just the mean + PC1

ex.temp.r1 = mtempfd + scores[1,1]*PCs[1]

#'  and plot these

plot(tempfd[1],ylim=c(-20,20))
lines(mtempfd,col=2)
lines(ex.temp.r1,col=3)

# Try adding the second PC

ex.temp.r2 = mtempfd + scores[1,1]*PCs[1] + scores[1,2]*PCs[2]
lines(ex.temp.r2,col=4)

# And the third

ex.temp.r3 = ex.temp.r2  + scores[1,3]*PCs[3]
lines(ex.temp.r3,col=6)

#'    * **THOUGHT EXERCISE**: how would you use this to choose the level of smoothing
#' in pca.fd by leave-one-curve-out cross validation?


#'
#' ##4. PCA of Multivariate Functions

#'  To look at two dimensions, we'll re-smooth the temperature data and look at
#'  it along with its derivative. To do that, let's consider a fourier basis
#'  and harmonic acceleration.

#'  First an fdPar object; we'll use 65 Fourier basis functions; more than enough
#'  but this will cut down the computational time. We'll also over-smooth so our
#'  derivatives don't look so wild.

fbasis =  create.fourier.basis(dayrng,65)
harmLfd = vec2Lfd(c(0,(2*pi/365)^2,0),rangeval=dayrng)
harm.fdPar = fdPar(fbasis,harmLfd,1e6)

#'  Do the smooth
tempSmooth2 = smooth.basis(daytime,temp,harm.fdPar)

#'  and take a derivative

temp.deriv = deriv.fd(tempSmooth2$fd)

#'  Now we need to duck under the hood to create a joint fd object for both
#'  temperature and precipitation at once.

#'  This basically means I need an fd object that stacks the coefficients for
#'  temperature and the coefficients for D temperature along a third dimension
#'  of a coefficient array.

Dtempcoefs = array(0,c(65,35,2))

Dtempcoefs[,,1] = tempSmooth2$fd$coefs
Dtempcoefs[,,2] = temp.deriv$coefs

#'  Let's also deal with the dimension names

Dtemp.fdnames = tempfd$fdnames
Dtemp.fdnames[[3]] = c('temperature','D temperature')

#'  Put it all together

Dtempfd = fd(Dtempcoefs,fbasis,Dtemp.fdnames)

#'  Now we can plot

par(mfrow=c(2,1))
plot(Dtempfd[,1])
plot(Dtempfd[,2])

#'  We can also look at covariance

Dtemp.varbifd = var.fd(Dtempfd)

#'  If we look at

Dtemp.var = eval.bifd(day5,day5,Dtemp.varbifd)

#'  We have a 74 x 74 x 1 x 3 array.

par(mfrow=c(1,1))

#'  First temperature
contour(day5,day5,Dtemp.var[,,1,1])

#'  Then d temperature
contour(day5,day5,Dtemp.var[,,1,3])

#'  Then their cross-product
contour(day5,day5,Dtemp.var[,,1,2])

#'    * **EXERCISE**: experiment with cor.fd on this multidimensional fd object.

Dtemp.cor = cor.fd(day5,Dtempfd)
#Dtemp.cor = eval.bifd(day5,day5,Dtemp.varbifd)
par(mfrow=c(1,1))
contour(day5,day5,Dtemp.cor[,,1,1])
contour(day5,day5,Dtemp.cor[,,1,3])
contour(day5,day5,Dtemp.cor[,,1,2])


## Now let's look at the PCA

Dtemp.pca = pca.fd(Dtempfd,nharm=4)

#'  The PCs are now two dimensional

par(mfrow=c(2,1))
plot(Dtemp.pca$harmonics[,1])
plot(Dtemp.pca$harmonics[,2])

#'  We can plot these in the usual manner

par(mfrow=c(1,2))
plot(Dtemp.pca,harm=1:3)

#'  But we can also plot the whole cycle

par(mfrow=c(1,1))

plot(Dtemp.pca,harm=1:3,cycle=TRUE,xlab='temperature',ylab='D temperature')

#'  **INTERPERATION**

#'  PC1 = over-all temperature and little variation in derivatives.

#'  PC2 = high summer-winter variation and also large variation in derivatives




