# import library
if (!require("pacman")) install.packages("pacman")
pacman::p_load(fda, data.table)
# read in data
df_ini <- read.csv('assignment/data_version5.csv')

''' (a) 
Provide one composite visualisation of the raw data as connected lines over 
the year for every single observation, using the matplot function. 
Make sure the x-axis reveals the actual time scale.
'''

# transpose
df <- transpose(df_ini)
# get row and colnames in order
colnames(df) <- rownames(df_ini)
rownames(df) <- colnames(df_ini)
# get number of time
ticks = seq(1,length(rownames(df)),1)
# make plot
matplot(df, type='l', xaxt='n')
# apply label and title
axis(1,at=ticks,labels=rownames(df))
title('Greenness of Trees of 200 Locations Line Plot')

''' (b) 
Use a saturated B-spline basis, spanning over the whole year to fit the data 
using a standard roughness penalty and choose the penalty paramater 
using cross validation. Comment on the smoothness of the curves and 
provide graphs of GCV and the final smoothed data.

'''

n <- nrow(df)
daytime = (1:n)-0.5
knots = day5 = seq(0,n,5)
dayrng = c(0,n)
t_mat <- as.matrix(df)
norder <- 4
n_basis = length(knots) + norder - 2
bbasis <- create.bspline.basis(dayrng,nbasis=n_basis,norder=4) 
tempSmooth1 <- smooth.basis(daytime,t_mat,bbasis)
plot(tempSmooth1)
curv.Lfd = int2Lfd(2)
harmLfd = vec2Lfd(c(0,(2*pi/n)^2,0),rangeval=dayrng)
#lambda = 1e6
#curv.fdPar = fdPar(bbasis,curv.Lfd,lambda)
#tempSmooth1 = smooth.basis(daytime,t_mat,curv.fdPar)
lambdas = 10^seq(-10,10,by=0.5)    #'  lambdas to look over
mean.gcv = rep(0,length(lambdas)) #'  store mean gcv
curv.Lfd = int2Lfd(2)

for(ilam in 1:length(lambdas)){
  curv.fdPar = fdPar(bbasis,curv.Lfd,ilam)
  # Set lambda
  curv.fdPari = curv.fdPar
  curv.fdPari$lambda = lambdas[ilam]
  # Smooth
  tempSmoothi = smooth.basis(daytime,t_mat,curv.fdPari)
  # Record average gcv
  mean.gcv[ilam] = mean(tempSmoothi$gcv)
}
# We can plot what we have
plot(lambdas,mean.gcv,type='b',log='x')

# Lets select the lowest of these and smooth
best = which.min(mean.gcv)
lambdabest = lambdas[best]
curv.fdPar = fdPar(bbasis,curv.Lfd,lambdabest)
# Set lambda
curv.fdPari = curv.fdPar
curv.fdPari$lambda = lambdas[lambdabest]
# fit smooth function with best lamda
tempSmoothi_best = smooth.basis(daytime,t_mat,curv.fdPari)
plot(tempSmoothi_best)

'''(c) 
Adjust your code to use a harmonic acceleration penalty with a period of 1 year. 
Choose an appropriate penalty paramater using cross validation. 
Comment on the smoothness of the curves and provide graphs of GCV and 
the final smoothed data.

'''
# ref: http://faculty.bscb.cornell.edu/~hooker/FDA2008/Lecture7_refinery.R
dayrng = c(0,365)
bbasis <- create.bspline.basis(dayrng,nbasis=n_basis,norder=4) 
harmLfd = vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))
for(ilam in 1:length(lambdas)){
  curv.fdPar = fdPar(bbasis,curv.Lfd,ilam)
  # Set lambda
  curv.fdPari = curv.fdPar
  curv.fdPari$lambda = lambdas[ilam]
  # Smooth
  tempSmoothi = smooth.basis(daytime,t_mat,curv.fdPari)
  # Record average gcv
  mean.gcv[ilam] = mean(tempSmoothi$gcv)
}
# We can plot what we have
plot(lambdas,mean.gcv,type='b',log='x')

# Lets select the lowest of these and smooth
best = which.min(mean.gcv)
lambdabest = lambdas[best]
curv.fdPar = fdPar(bbasis,curv.Lfd,lambdabest)
# Set lambda
curv.fdPari = curv.fdPar
curv.fdPari$lambda = lambdas[lambdabest]
# Smooth
tempSmoothi_best = smooth.basis(daytime,t_mat,curv.fdPari)
plot(tempSmoothi_best)


'''(d) 
Based on the appropriate harmonic acceleration penalty fit calculate and 
plot the graphs of the first and second and derivatives of the curves.
'''
# first derivatives
plot(deriv.fd(tempSmoothi_best$fd))

# second derivatives
plot(deriv.fd(tempSmoothi_best$fd,2))

'''(e) 
Conduct a un-penalized principal components analysis of these data. 
How many components do you need to recover 80% of the variation? 
Do the components appear satisfactory?
'''
bbasis<-create.bspline.basis(dayrng,norder=4,nbasis=4)
harmfdPar <- fdPar(bbasis, harmLfd, lambda=0)
pinch.smooth<-smooth.basis(daytime,t_mat,bbasis)
pinch.pca<-pca.fd(pinch.smooth$fd,nharm = 4)
plot(pinch.pca$harmonics[1:4],lty=1)

par(mfrow=c(1,1))
plot(cumsum(pinch.pca$varprop),ylab="cumulative sum explained")
abline(h=0.9)

'''(f) 
Try a smoothed PCA analysis from the raw data. 
Choose the smoothing parameter by cross-validation. 
Plot the cross-validation curve. Plot the new smoothed principal components. 
Does this appear to be more satisfactory than the unsmoothed version?
'''
#' ##1. pca.fd
curv.Lfd = int2Lfd(2)
bbasis = create.bspline.basis(dayrng,nbasis=4,norder=4) 
tempSmooth = smooth.basis(daytime,t_mat,bbasis)
tempfd = tempSmooth$fd
pca.fdPar = fdPar(bbasis,curv.Lfd,1e4)
#tempPCAsmooth = pca.fd(tempfd,nharm=6,harmfdPar=pca.fdPar)

ilam = 1
mean.gcv = NULL 
lambdas = 10^(seq(2,6,0.5))
for(ilam in 1:length(lambdas)){
  #curv.fdPar = fdPar(bbasis,curv.Lfd,ilam)
  # Set lambda
  #curv.fdPari = curv.fdPar
  #curv.fdPari$lambda = lambdas[ilam]
  # Smooth
  #tempSmoothi = smooth.basis(daytime,t_mat,curv.fdPari)
  # Record average gcv
  #mean.gcv[ilam] = mean(tempSmoothi$gcv)
  
  
  
  pca.fdPari = fdPar(bbasis,curv.Lfd,lambdas[ilam])
  tempPCAsmoothi = pca.fd(tempfd,nharm=4,harmfdPar=pca.fdPari)
  # compute mean of score?
  Xmat.smoothi = tempPCAsmoothi$scores[,1:3]
  mean.gcv[ilam] = mean(Xmat.smoothi)
}
# We can plot what we have
plot(lambdas,mean.gcv,type='b',log='x')


'''(g) 
Provide a interpretation for the smoothed principal components.
'''


'''(h) 
Using the inprod function from the fda library numerically show that 
the un-penalized principal components are orthogonal, 
but the penalized principal components might not be. 
You can just use the first 4 principal components in each case. 
Ignore any "convergence" warning while executing the inprod function.
'''
# ref: https://mathworld.wolfram.com/OrthogonalMatrix.html
# check: A%*%t(A) = I