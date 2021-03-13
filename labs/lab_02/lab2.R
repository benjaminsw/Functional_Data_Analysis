
#' ---
#' title: "Lab 2: Smoothing "
#' author: "Surajit Ray"
#' 
#' output:
#'    html_document:
#'      toc: true
#'      toc_float: true
#'      highlight: zenburn
#' ---

#'
#' ### Loading R library
###' Firstly, let's load up the library and the data used in Lab 1.
par(ask=FALSE)
library('fda')

data(CanadianWeather)

temp = CanadianWeather$dailyAv[,,1]
precip = CanadianWeather$dailyAv[,,2]

daytime = (1:365)-0.5
day5 = seq(0,365,5)

dayrng = c(0,365)


################################################################################
#' #                           SIMPLE FIT without smoothing
################################################################################

bbasis = create.bspline.basis(dayrng,nbasis=21,norder=4) # or from Lab 1

tempSmooth1 = smooth.basis(daytime,temp,bbasis)

plot(tempSmooth1)
################################################################################
#' #                           SMOOTHING FUNCTIONS
################################################################################

#' ## 1. Lfd Objects

#'  Two common means of generating Lfd objects

#'
#' ### 1. int2Lfd -- just penalize some derivatives.

curv.Lfd = int2Lfd(2)

#'
#' ### 2. vec2Lfd -- a (constant) linear combination of derivatives; for technical
#'  reasons this also requires the range of the basis.

harmLfd = vec2Lfd(c(0,(2*pi/365)^2,0),rangeval=dayrng)

#'  looking inside these objects is not terribly enlightening.

#'
#' ##  2. fdPar objects

#'  We'll concentrate on B-splines and second-derivative penalties.

#'  First, a value of lambda  (purposefully large so that we can distinguish a fit
#'  from data below).

lambda = 1e6

#'  Now we can define the fdPar object

curv.fdPar = fdPar(bbasis,curv.Lfd,lambda)

#' ##  3. Smoothing functions

#'  We're now in a position to smooth

tempSmooth1 = smooth.basis(daytime,temp,curv.fdPar)

#'  Let's look at the result

names(tempSmooth1)

#'  First of all, let's plot it

plot(tempSmooth1$fd)

#'  There is also a neat utility to go through each curve in turn and look at its
#'  fit to the data: You can plot it without spliting the screen
par(mfrow=c(2,2))
plotfit.fd(temp,daytime,tempSmooth1$fd,index=1:4)



#'  Let's examine some fit statistics

#'  degrees of freedom

tempSmooth1$df

#'  Just about equivalent to fitting 5 parameters

#'  We'll also look at GCV, this is given for each observation

tempSmooth1$gcv

#'  Let's change to a more realistic value of lambda

lambda = 1e1
curv.fdPar$lambda = lambda

tempSmooth = smooth.basis(daytime,temp,curv.fdPar)

#'  and repeat the previous steps
par(mfrow=c(2,2))
plotfit.fd(temp,daytime,tempSmooth$fd,index=1:4)
tempSmooth$df
tempSmooth$gcv

#'  Here the fit looks a lot better and the gcv values are much smaller.

#'
#' ##  4. Choosing smoothing parameters


#'  We can search through a collection of smoothing parameters to try and find
#'  an optimal parameter.

#'  We will record the average gcv and choose lambda to be the minimum of these.

lambdas = 10^seq(-4,4,by=0.5)    #'  lambdas to look over

mean.gcv = rep(0,length(lambdas)) #'  store mean gcv


for(ilam in 1:length(lambdas)){
  #'  Set lambda
  curv.fdPari = curv.fdPar
  curv.fdPari$lambda = lambdas[ilam]

  #'  Smooth
  tempSmoothi = smooth.basis(daytime,temp,curv.fdPari)

  #'  Record average gcv
  mean.gcv[ilam] = mean(tempSmoothi$gcv)
}

#'  We can plot what we have

plot(lambdas,mean.gcv,type='b',log='x')

#'  Lets select the lowest of these and smooth

best = which.min(mean.gcv)
lambdabest = lambdas[best]



curv.fdPar$lambda = lambdabest
tempSmooth = smooth.basis(daytime,temp,curv.fdPar)

#'  And look at the same statistics
par(mfrow=c(2,2))
plotfit.fd(temp,daytime,tempSmooth$fd,index=1:4)
tempSmooth$df

#'  We'll also plot these

plot(tempSmooth)
#'
#' ##   EXERCISE

#'    *   **EXERCISE**: try obtaining a smooth of the precipitation data

#'    *   **EXERCISE**: how much does the result change if the basis has a knot every day
###'  instead of every 5 days?


################################################################################
#' #        FUNCTIONAL DATA OBJECTS: MANIPULATION AND STATISTICS              ####
################################################################################


##'  Now that we have a functional data object we can manipulate them in various
#'  ways.  First let's extract the fd object

tempfd = tempSmooth$fd

#'  if we look at what's in this we see

names(tempfd)

#'  We see a basis, plus coefficient matrix

dim(tempfd$coefs)

#'  and an array giving names

tempfd$fdnames

#'  With lists giving names for time points, replicates and dimensions. Each list
#'  also has a name that can be used in plotting.  Apart from plotting functions,
#'  fdnames isn't used and you can generally ignore it.

#'  We can also create fd objects by adding a basis and a coefficient array. Let's
#'  make a random one, say
nbasis = 21
newcoefs = matrix(rgamma(nbasis*10,5,2),nbasis,10)
newfd = fd(newcoefs,bbasis)

#'  Notice that we haven't specified fdnames.

#'  The plotting command nicely draws these.

plot(newfd)

#'  Not that this looks very nice; we'll stick with the Canadian weather data.

#'
#' ##  1. Manipulation

#'  We can do a number of things with these functions, treating them as data.
#'  These operations all result in new functional data objects, but we will plot
#'  them directly as an illustration.

#'  Subset

plot(tempfd[1:10])

#'  We can add them together; the 'lines' function also works with them

newfd = tempfd[1] + tempfd[2]
plot(newfd)
lines(tempfd[1],col=2)
lines(tempfd[2],col=4)

#'  We can also multiply

plot(tempfd[1]*tempfd[2])

#'  And obtain powers

plot(tempfd[1]^2)

#'  We can also obtain derivatives

plot(deriv.fd(tempfd))

#'  These are pretty wild because of the roughness of the resulting curves
#'  instead let's have a look at the over-smoothed data:

plot(deriv.fd(tempSmooth1$fd))

#'  We can also look at second derivatives

plot(deriv.fd(tempSmooth1$fd,2))

#'  Note that it is a property of B-splines of order m that the (m-2)th derivative
#'  is zero at the end of the interval.
#'
#' ##  2. Summary statistics

#'  The obvious thing to look at is the mean

mtempfd = mean(tempfd)
plot(tempfd,col=4)
lines(mtempfd,lwd=2,col=2)

#'  We can also examine a variance

temp.varbifd = var.fd(tempfd)

#'  temp.varbifd is a bivariate functional data object -- meaning it takes values
#'  on a rectangle.

#'  To plot this, we need to evaluate it; here we'll use day5 -- 365 points is
#'  a bit overkill.

temp.var = eval.bifd(day5,day5,temp.varbifd)
contour(day5,day5,temp.var)

#'  Mostly high variance in the winter, low in summer. Let's have a look at
#'  correlation. In this case, evaluation arguments go in with the function call

temp.cor = cor.fd(day5,tempfd)
filled.contour(day5,day5,temp.cor)

#'  Here we see high correlation between Summer and Winter temperatures, but
#'  much less in the spring and fall (although spring to fall correlation is still
#'  high).

#'
#' ##   EXERCISE


#'    *  **EXERCISE**: obtain these for the precipitation data and look at the covariance
###'  and correlation between temperature and precipitation.

#'    *   **EXERCISE**: try repeating the above with a Fourier basis and the harmonic
