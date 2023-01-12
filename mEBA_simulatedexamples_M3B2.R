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

set.seed(567)

########################
## 0. Import functions and settings
########################
Rcpp::sourceCpp('mEBA_CPPfunctions.cpp')
source('mEBA_Rfunctions.R')


t <- 200; #length of time series
R <- 10; #number of components

########################
## Differing Proportions 3 Band (M3B-2)
########################
#simulate data
X.m3b2 <- matrix(NA,nrow=t,ncol=R);
df <- meba.simdata(t+200);
ll <- seq(from=0,by=-1,length.out=R); #ll as c(0,-1,-2,-3,..)
cf <- rep(1,R); #same for all cf as 1

for (m in 1:floor(R*0.2)){
  X.m3b2[,m] <- cf[m]*df$bL2f15[(101+ll[m]):(t+100+ll[m])]
}      

for (m in (floor(R*0.2)+1):R){
  X.m3b2[,m] <- cf[m]*df$bS2f35[(101+ll[m]):(t+100+ll[m])]
}

#plot series 
plot.ts(X.m3b2,main="Differing Proportions 3 Bands")

#compute and plot local periodogram and demeaned local periodogram
N <- 2*floor(t^0.7)-floor(t^0.7/2)*2; #neighborhood for local periodogram
freq <- seq(0,floor(N/2),by=1)/N
pse <- fhat(X.m3b2,N,stdz=FALSE);
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
m3b2.out=msboot(nrep=1000, X.m3b2, Wsel=3, stdz=FALSE, ncore=1)

#plot output
#plot of observed test statistics across frequencies (red) and the first 50 bootstrap test statistics (black) for each W
par(mfrow=c(3,1))
plot(m3b2.out[[2]][[1]][,1:2],type='l',xlab='Frequency',ylab='D(omega)',main='W1',ylim=c(0,90000))
apply(m3b2.out[[2]][[1]][,3:50],2,function(x) lines(cbind(m3b2.out[[2]][[1]][,1],x)))
lines(m3b2.out[[1]][[1]],col='red',lwd=2)

plot(m3b2.out[[2]][[2]][,1:2],type='l',xlab='Frequency',ylab='D(omega)',main='W2',ylim=c(0,90000))
apply(m3b2.out[[2]][[2]][,3:50],2,function(x) lines(cbind(m3b2.out[[2]][[2]][,1],x)))
lines(m3b2.out[[1]][[2]],col='red',lwd=2)

plot(m3b2.out[[2]][[3]][,1:2],type='l',xlab='Frequency',ylab='D(omega)',main='W3',ylim=c(0,90000))
apply(m3b2.out[[2]][[3]][,3:50],2,function(x) lines(cbind(m3b2.out[[2]][[3]][,1],x)))
lines(m3b2.out[[1]][[3]],col='red',lwd=2)

#plot of p-values for testing each frequency as a partition point (black) and 0.05 threshold (red) for each W
par(mfrow=c(3,1))
plot(m3b2.out[[3]][[1]],type='l',xlab='Frequency',ylab='p-value',main='W1',ylim=c(0,1));
abline(h=0.05,col='red')
plot(m3b2.out[[3]][[2]],type='l',xlab='Frequency',ylab='p-value',main='W1',ylim=c(0,1));
abline(h=0.05,col='red')
plot(m3b2.out[[3]][[3]],type='l',xlab='Frequency',ylab='p-value',main='W1',ylim=c(0,1));
abline(h=0.05,col='red')

#table of significant frequency partition points (none in the white noise case)
m3b2.out[[4]][which(m3b2.out[[4]][,2]==1),1]


#local periodogram with estimated cutpoints
idx=1 #component for which you want to view local periodogram (can be 1,2,..,R)
par(mfrow=c(1,1))
image.plot(x=(1:t)/(t+1),y=freq[-1],z=t(Re(pse[-1,idx+(idx-1)*R,])), 
           axes = TRUE, col = inferno(256),
           main = 'Local Periodogram',xlab='Time',ylab='Frequency',xaxs="i");
abline(h=m3b2.out[[4]][which(m3b2.out[[4]][,2]==1),1],col='green',lwd=2);

