rm(list=ls())

library(RcppArmadillo)
library(ggplot2)
library(signal)
library(fields)
library(viridis)

setwd('/Update/Path/To/Directory/Here/GitHub_repo')

#output from msboot function is a list of the following 7 objects:
#1) A list of matrices containing the test statistics over the Fourier frequencies for each W
#2) A list of matrices containing the bootstrap test statistics over the Fourier frequencies for each W
#3) A list of matrices containing the p-values for testing each Fourier frequency as a partition point for each W
#4) A matrix indicating whether or not each frequency was identified as a partition point
#5) A list of matrices containing the test statistics for each component over the Fourier frequencies for each W
#6) A list of arrays containing the bootstrap test statistics for each component over the Fourier frequencies for each W
#7) A list of matrices containing the p-values for each component testing each Fourier frequency as a partition point for each W

set.seed(489)

########################
## 0. Import functions and settings
########################
Rcpp::sourceCpp('mEBA_CPPfunctions.cpp')
source('mEBA_Rfunctions.R')


t <- 200; #length of time series
R <- 10; #number of components

########################
## Sinusoidal 3 Band (S3B)
########################
#simulate data
X.s3b <- matrix(NA,nrow=t,ncol=R);
df <- meba.simdata(t+200);
ll <- seq(from=0,by=-1,length.out=R); #ll as c(0,-1,-2,-3,..)
cf <- rep(1,R); #same for all cf as 1

for (m in 1:R){
  X.s3b[,m] <- cf[m]*df$bS[(101+ll[m]):(t+100+ll[m])]
}

#plot series 
plot.ts(X.s3b,main="Sinusoidal 3 Bands")

#compute and plot local periodogram and demeaned local periodogram
N <- 2*floor(t^0.7)-floor(t^0.7/2)*2; #neighborhood for local periodogram
freq <- seq(0,floor(N/2),by=1)/N
pse <- fhat(X.s3b,N,stdz=FALSE);
gpse <- ghat(pse);

idx=1 #component for which you want to view local periodogram (can be 1,2,..,R)
par(mfrow=c(1,1))
image.plot(x=(1:t)/(t+1),y=freq[-1],z=t(Re(pse[-1,idx+(idx-1)*R,])), 
           axes = TRUE, col = inferno(256),
           main = 'Local Periodogram',xlab='Time',ylab='Frequency',xaxs="i");

par(mfrow=c(1,1))
image.plot(x=(1:t)/(t+1),y=freq[-1],z=t(Re(gpse[-1,idx+(idx-1)*R,])), 
           axes = TRUE, col = inferno(256),
           main = 'Demeaned Local Periodogram',xlab='Time',ylab='Frequency',xaxs="i");

#run algorithm and store results
s3b.out=msboot(nrep=1000, X.s3b, Wsel=3, stdz=FALSE, ncore=1)

#plot output
#plot of observed test statistics across frequencies (red) and the first 50 bootstrap test statistics (black) for each W
par(mfrow=c(3,1))
plot(s3b.out[[2]][[1]][,1:2],type='l',xlab='Frequency',ylab='D(omega)',main='W1',ylim=c(0,200000))
apply(s3b.out[[2]][[1]][,3:50],2,function(x) lines(cbind(s3b.out[[2]][[1]][,1],x)))
lines(s3b.out[[1]][[1]],col='red',lwd=2)

plot(s3b.out[[2]][[2]][,1:2],type='l',xlab='Frequency',ylab='D(omega)',main='W2',ylim=c(0,200000))
apply(s3b.out[[2]][[2]][,3:50],2,function(x) lines(cbind(s3b.out[[2]][[2]][,1],x)))
lines(s3b.out[[1]][[2]],col='red',lwd=2)

plot(s3b.out[[2]][[3]][,1:2],type='l',xlab='Frequency',ylab='D(omega)',main='W3',ylim=c(0,200000))
apply(s3b.out[[2]][[3]][,3:50],2,function(x) lines(cbind(s3b.out[[2]][[3]][,1],x)))
lines(s3b.out[[1]][[3]],col='red',lwd=2)

#plot of p-values for testing each frequency as a partition point (black) and 0.05 threshold (red) for each W
par(mfrow=c(3,1))
plot(s3b.out[[3]][[1]],type='l',xlab='Frequency',ylab='p-value',main='W1',ylim=c(0,1));
abline(h=0.05,col='red')
plot(s3b.out[[3]][[2]],type='l',xlab='Frequency',ylab='p-value',main='W1',ylim=c(0,1));
abline(h=0.05,col='red')
plot(s3b.out[[3]][[3]],type='l',xlab='Frequency',ylab='p-value',main='W1',ylim=c(0,1));
abline(h=0.05,col='red')

#table of significant frequency partition points (none in the white noise case)
s3b.out[[4]][which(s3b.out[[4]][,2]==1),1]


#local periodogram with estimated cutpoints
idx=1 #component for which you want to view local periodogram (can be 1,2,..,R)
par(mfrow=c(1,1))
image.plot(x=(1:t)/(t+1),y=freq[-1],z=t(Re(pse[-1,idx+(idx-1)*R,])), 
           axes = TRUE, col = inferno(256),
           main = 'Local Periodogram',xlab='Time',ylab='Frequency',xaxs="i");
abline(h=s3b.out[[4]][which(s3b.out[[4]][,2]==1),1],col='green',lwd=2);

