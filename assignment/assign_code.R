# import library
if (!require("pacman")) install.packages("pacman")
pacman::p_load(fda, data.table)
# read in data
df_ini <- read.csv('assignment/data_version5.csv')

# (a) Provide one composite visualisation of the raw data as connected lines over the year for
# every single observation, using the matplot function. Make sure the x-axis reveals the
# actual time scale.


# transpose dataframe
df <- as.data.frame(t(as.matrix(df_ini)))
# get row name
rownames(df) <- colnames(df_ini)
# get number of time
ticks <- seq(1,length(rownames(df)),1)
# plot data 
matplot(df,type='l', xaxt='n')
# apply label and title
axis(1,at=ticks,labels=rownames(df))
title('Greenness of Trees of 200 Locations Plot')

# (b) Use a saturated B-spline basis, spanning over the whole year to fit the data using a
# standard roughness penalty and choose the penalty para mater using cross validation.
# Comment on the smoothness of the curves and provide graphs of GCV and the final
# smoothed data.

# set up params
# number of rows/days
days_no <- nrow(df)
# days sequence
day_seq <- 1:days_no - 0.5
# day range
day_range <- c(0,days_no)
# observed sampling
obs_ticks <- seq(0,days_no,1)
# create saturated knots 
knots <- c(seq(0,days_no,1))
norder <- 4
nbasis <- length(knots) + norder - 2
# create saturated B-spline basis
saturated_basis <- create.bspline.basis(day_range,nbasis,norder,knots)
##### sanity check #####
# number of basis functions
saturated_basis$nbasis
# basis range
saturated_basis$rangeval   
# plot 
plot(saturated_basis)

# eval basis values at each sampling
basismat <- eval.basis(obs_ticks, saturated_basis);
dim(basismat)
# plot basis at each sampling
plot(obs_ticks,basismat[,1],type = "l",col=1,lwd=3)
for(i in 1:49) lines(obs_ticks,basismat[,i],type = "l",col=i+1,lwd=3)

# range of lambdas to search for optimal
lambdas <- 10^seq(-10,10,by=0.5)    
# keep track of generalized cross-validation (GCV) mean
gcv_mean <- rep(0,length(lambdas)) 
# convert integer to linear differential operator and define the m-th order derivative penalty term
int_to_linear_diff <- int2Lfd(2)

for(ilam in 1:length(lambdas)){
  # define a functional parameter object
  func_param_obj <- fdPar(saturated_basis, Lfdobj=int_to_linear_diff, lambda=lambdas[ilam])
  # smooth the data using the roughness penalty and smoothing parameter specified in 'func_param_obj' 
  smooth_func <- smooth.basis(day_seq, as.matrix(df), func_param_obj)
  # record average gcv
  gcv_mean[ilam] <- mean(smooth_func$gcv)
}
# We can plot what we have
plot(lambdas, gcv_mean, type='b', log='x')

# Lets select the lowest of these and smooth
best <- which.min(gcv_mean)
#best_lambda <- lambdas[best]
best_func_param_obj <- fdPar(saturated_basis, Lfdobj=int_to_linear_diff, lambda=lambdas[best])
# fit smooth functionn with the best lambda
best_smooth <- smooth.basis(day_seq, as.matrix(df), best_func_param_obj)
plot(best_smooth)

##########################
logl <- seq(-5, 12, len=71)  
range(exp(logl))
gcv <- rep(0,71)

for(i in c(1:length(logl))){
  lambda <- exp(logl[i])
  
  tD2fdPar <- fdPar(saturated_basis,Lfdobj=int2Lfd(2),lambda=lambda)
  tyfd <- smooth.basis(day_seq, as.matrix(df),tD2fdPar)
  
  gcv[i] = tyfd$gcv
}

# PLOT GCV of FIT versus log lambda
plot(logl,gcv[1:71],type='l',cex.lab=1.5, lwd=4, 
     xlab='log lambda',ylab='GCV', main="GCV(log.lambda)")
#####################


# (c) Adjust your code to use a harmonic acceleration penalty with a period of 1 year. Choose
# an appropriate penalty para mater using cross validation. Comment on the smoothness
#of the curves and provide graphs of GCV and the final smoothed data.


# create saturated fourier basis functions
fourier_basis <- create.fourier.basis(day_range,45)
# make a linear differential operator object from a vector
fourier_diff_operator_obj <- vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))
# range of lambdas to search for optimal
lambdas <- 10^seq(-10,10,by=0.5)    
# keep track of generalized cross-validation (GCV) mean
gcv_mean <- rep(0,length(lambdas)) 


for(ilam in 1:length(lambdas)){
  # define a functional parameter object
  func_param_obj <- fdPar(fourier_basis, Lfdobj=fourier_diff_operator_obj, lambda=lambdas[ilam]) #**bspline or fourior???
  # smooth the data using the roughness penalty and smoothing parameter specified in 'func_param_obj' 
  smooth_func <- smooth.basis(day_seq, as.matrix(df), func_param_obj)
  # record average gcv
  gcv_mean[ilam] <- mean(smooth_func$gcv)
}
# We can plot what we have
plot(lambdas, gcv_mean, type='b', log='x')

# Lets select the lowest of these and smooth
best <- which.min(gcv_mean)
#best_lambda <- lambdas[best]
best_func_param_obj <- fdPar(saturated_basis, Lfdobj=int_to_linear_diff, lambda=lambdas[best])
# fit smooth functionn with the best lambda
best_smooth <- smooth.basis(day_seq, as.matrix(df), best_func_param_obj)
plot(best_smooth)

# (d) Based on the appropriate harmonic acceleration penalty fit calculate and plot the graphs
# of the first and second and derivatives of the curves.

# first derivatives
plot(deriv.fd(best_smooth$fd))

# second derivatives
plot(deriv.fd(best_smooth$fd,2))

# (e) Conduct a un-penalized principal components analysis of these data. 
# How many components do you need to recover 80% of the variation? 
# Do the components appear satisfactory?

#int_to_linear_diff <- int2Lfd(2)
#ticks <- seq(1,length(rownames(df)),1)
nbasis <- 45 #length(ticks) + norder - 2 
saturated_basis <- create.bspline.basis(day_range,nbasis,norder)
#func_param_obj <- fdPar(saturated_basis, Lfdobj=int_to_linear_diff)
# smooth the data using the roughness penalty and smoothing parameter specified in 'func_param_obj' 
smooth_func <- smooth.basis(day_seq, as.matrix(df), saturated_basis)
# average gcv
mean(smooth_func$gcv)
# fit functional principal components analysis
func_pca <- pca.fd(smooth_func$fd,nharm=10)
non_smoothed_func_pca <- func_pca
names(func_pca)
func_pca$varprop
#func_pca$values are the eigenvalues
plot(func_pca$values[1:8],xlab='component',ylab='variance',col="red",
     cex.lab=1.5,cex.axis=1.5,cex=2)

# plot the cumulative percentage explained total variations
plot(cumsum(func_pca$values[1:10])/sum(func_pca$values),xlab='Number of Components',
     ylab='cumulative variance explained',col=2,cex.lab=2,
     cex.axis=2,cex=2)
abline(h=0.80)
# plot mean function
plot(func_pca$meanfd)

# functional principal components
pca_fd <- func_pca$harmonics
pca_vals <- eval.fd(day_seq, pca_fd)
# the top 4 FPCs
dim(pca_vals) 
# plot 4 pca
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(day_seq, pca_vals[,1:4], xlab='day', ylab='PCs',
        lwd=4,lty=1,cex.lab=2.5,cex.axis=2.5,type='l')
legend(0,-0.07,c('PC1','PC2','PC3','PC4'),col=1:4,lty=1,lwd=5)
title('Principle Component Functions')

# (f) Try a smoothed PCA analysis from the raw data. Choose the smoothing parameter
# by cross-validation. Plot the cross-validation curve. Plot the new smoothed principal
# components. Does this appear to be more satisfactory than the unsmoothed version?


# create saturated fourier basis functions
# fourier_basis <- create.fourier.basis(day_range,45)
# make a linear differential operator object from a vector
fourier_diff_operator_obj <- vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))
# range of lambdas to search for optimal
lambdas <- 10^seq(-10,10,by=0.5)    
# keep track of generalized cross-validation (GCV) mean
gcv_mean <- rep(0,length(lambdas)) 
for(ilam in 1:length(lambdas)){
  # define a functional parameter object
  func_param_obj <- fdPar(fourier_basis, Lfdobj=fourier_diff_operator_obj, lambda=lambdas[ilam]) #**bspline or fourior???
  # smooth the data using the roughness penalty and smoothing parameter specified in 'func_param_obj' 
  smooth_func <- smooth.basis(day_seq, as.matrix(df), func_param_obj)
  # record average gcv
  gcv_mean[ilam] <- mean(smooth_func$gcv)
}
# We can plot what we have
plot(lambdas, gcv_mean, type='b', log='x')

# Lets select the lowest of these and smooth
best <- which.min(gcv_mean)
#best_lambda <- lambdas[best]
best_func_param_obj <- fdPar(saturated_basis, Lfdobj=int_to_linear_diff, lambda=lambdas[best])
# fit smooth functionn with the best lambda
best_smooth <- smooth.basis(day_seq, as.matrix(df), best_func_param_obj)
plot(best_smooth)

# fit functional principal components analysis
func_pca <- pca.fd(best_smooth$fd,nharm=10)
smoothed_func_pca <- func_pca
names(func_pca)
func_pca$varprop

# plot the cumulative percentage explained total variations
plot(cumsum(func_pca$values[1:10])/sum(func_pca$values),xlab='Number of Components',
     ylab='cumulative variance explained',col=2,cex.lab=2,
     cex.axis=2,cex=2)
abline(h=0.80)
# plot mean function
plot(func_pca$meanfd)

# functional principal components
pca_fd <- func_pca$harmonics
pca_vals <- eval.fd(day_seq, pca_fd)
# the top 4 FPCs
dim(pca_vals) 
# plot 4 pca
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(day_seq, pca_vals[,1:4], xlab='day', ylab='PCs',
        lwd=4,lty=1,cex.lab=2.5,cex.axis=2.5,type='l')
legend(0,-0.07,c('PC1','PC2','PC3','PC4'),col=1:4,lty=1,lwd=5)
title('Principle Component Functions')

# (g) Provide a interpretation for the smoothed principal components

# (h) Using the inprod function from the fda library numerically show that the un-penalized
# principal components are orthogonal, but the penalized principal components might not be. You can
# just use the first 4 principal components in each case. Ignore any "convergence" warning while
# executing the inprod function. 
non_smoothed_inprod <- inprod(non_smoothed_func_pca$harmonics, non_smoothed_func_pca$harmonics)
non_smoothed_inprod <- 1*(non_smoothed_inprod > 0.99)
smoothed_inprod <- inprod(smoothed_func_pca$harmonics, smoothed_func_pca$harmonics)
smoothed_inprod <- 1*(smoothed_inprod > 0.99)
# proof: inner product is 0
