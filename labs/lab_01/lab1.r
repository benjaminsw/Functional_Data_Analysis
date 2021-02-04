
#' ---
#' title: "Lab 1: Basis Functions and fitting"
#' author: "Surajit Ray"
#' output:
#'    pdf_document:
#'      highlight: zenburn
#' ---

#'
#' ### Loading R library
###' Firstly, let's load up the library and some data.

library('fda')


#' ### Loading dataset 'CanadianWeather'

#' All of our examples will come from the Canadian Weather Data. Not all of these
#' examples will therefore be particularly realistic. But they are there in order
#' to demonstrate aspects of the code. Realism is retained as much as possible.


#' We can load these data by

data(CanadianWeather)

#' Temperature and precipitation are contained in the dailyAV element, and we'll
#' extract them.

temp = CanadianWeather$dailyAv[,,1]
precip = CanadianWeather$dailyAv[,,2]

#' We need corresponding time points. I'll put them half-way through a day. This
#' is because the period is over 0:365 and we'd like the 365 data points to be
#' about symmetric in that period.

daytime = (1:365)-0.5

#' This is a bit fine for plotting purposes, so we'll also create a vector of
#' points to use for plotting every 5 days

day5 = seq(0,365,5)


#' Plot these data

#' Note: I am assuming you will use only one plotting window and won't close it
#' so I haven't insisted on re-setting the window as being two-pane below.

matplot(daytime,temp,type='l')
matplot(daytime,precip,type='l')

#' We can also plot by region; Atlantic, Pacific, Central and North.

matplot(daytime,temp,type='l',col=as.factor(CanadianWeather$region))
matplot(daytime,precip,type='l',col=as.factor(CanadianWeather$region))


#' But let's get on to fda.


#'
#' #  DEFINING A BASIS SYSTEM

#'  We'll always need the range

dayrng = c(0,365)

#' ## 1. Fourier Basis with 365 basis functions

fbasis = create.fourier.basis(dayrng,365)

#' Plot a basis with just 5 components

plot(create.fourier.basis(dayrng,5))

#' etc


#' etc

#' Let's try a simple linear regression of the first temperature record on the
#' first, say, 20 basis functions

fb.values = eval.basis(day5,fbasis)

#' the 74 by 365 matrix that results has rows as days, columns as bases

dim(fb.values)
plot(day5,fb.values[,1])
plot(day5,fb.values[,2])

#' Extract the first temperature record

ex.temp = temp[,1]
plot(ex.temp,main='Daily Temperature of St. Johns')
#' Run a linear regression on temperature with the first 5 basis functions. In
#' this case I need to evaluate the basis at the observation times

Xmat = eval.basis(daytime,fbasis)
Xmat = Xmat[,1:5]   #' First 5 basis functions

matplot(Xmat)

ex.temp.mod = lm(ex.temp~Xmat)

#' Let's plot this; the fitted values returned by lm will represent the smooth
#' well enough to plot.

plot(daytime,ex.temp)
lines(daytime,ex.temp.mod$fitted,col=2,lwd=2)
title(main='Daily Temperature of St. Johns with smoothing')

#' We can also look at residuals

plot(daytime,ex.temp.mod$resid)

#' There's some clear autocorrelation; more basis functions may be warranted.

#'
#' * EXERCISE: Repeat the above with different numbers of basis functions (say
#'  20, 50, 100, 200, 365. How many look like they give a reasonable smooth?




#'
#' ## 2. B-spline bases with knots every 5 days

#' First of all define a knot sequence; this will be the same as day5

knots = day5

#' We'll use fourth-order B-splines

norder = 4

#' this implies the number of basis functions

nbasis = length(knots) + norder - 2

#' Now we can define the basis

bbasis = create.bspline.basis(dayrng,nbasis,norder,knots)

#' If in doubt, we can obtain

bbasis$nbasis    #' number of basis functions
bbasis$rangeval   #' basis range

#' Plotting these is a little crazy

plot(bbasis)

#' but we can look at a smaller number

plot(create.bspline.basis(dayrng,nbasis=12,norder))

#' We can also look at the inner product of these

in.mat = inprod(bbasis,bbasis)

par(mfrow=c(1,1))
image(in.mat)

#' and see that it is zero outside of a diagonal band; this can help computation
#' a great deal.

#'
#' * EXERCISE: try changing the order of the basis and observe how the width
#' of the support of the basis changes and how its smoothness properties change.
#' * EXERCISE: obtain a least squares smooth of these data with the Bspline basis
#' how does this compare with a least squares smooth using a Fourier basis with
#' the same number of basis functions?







