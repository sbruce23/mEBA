rm(list=ls())

library(RcppArmadillo)
library(ggplot2)
library(R.matlab)
library(signal)
library(fields)
library(viridis)
library(parallel)
library(ggpubr)

#setwd('/Update/Path/to/Local/Directory/Here/GitHub_repo')

ncores <- ceiling(detectCores()/2) #set number of cores appropriately

########################
## 0. Import functions
########################
Rcpp::sourceCpp('mEBA_CPPfunctions.cpp')
source('mEBA_Rfunctions.R')


########################
## 1. Import EEG signals for participants 2 and 13
########################

eeg_p2 <- readMat(url("https://zenodo.org/record/2348892/files/subject_03.mat"))
eeg_p13 <- readMat(url("https://zenodo.org/record/2348892/files/subject_17.mat"))
# NOTE: Six subjects (2, 7, 9, 10, 19, and 20) were excluded due to data quality
# concerns.  Remaining participants are renumbered 1 to 14.

#EO: eyes open, EC: eyes closed
#subset series to five blocks: EC->EO->EC->EO->EC
ec_p2 <- which(eeg_p2$SIGNAL[,18]==1)[1] + 512 #eyes closed first instance
eo_p2 <- which(eeg_p2$SIGNAL[,19]==1)[3] + 512 #last instance of eyes open
eeg_p2 <- eeg_p2$SIGNAL[ec_p2:eo_p2,2:17]

ec_p13 <- which(eeg_p13$SIGNAL[,18]==1)[1] + 512 #eyes closed first instance
eo_p13 <- which(eeg_p13$SIGNAL[,19]==1)[3] + 512 #last instance of eyes open
eeg_p13 <- eeg_p13$SIGNAL[ec_p13:eo_p13,2:17]

#remove linear trends, standardized variance, downsample to 64Hz, 
#and filter out frequencies near 0 (below 1 Hz)
dsfrq <- 64 # down sample frequency
#EEG channels
channels <- c("FP1","FP2","FC5","FC6","FZ","T7","CZ","T8","P7","P3","PZ","P4",
              "P8","O1","OZ","O2") 
eeg_p2 <- preprocess(eeg_p2,dsfrq,channels)
eeg_p13 <- preprocess(eeg_p13,dsfrq,channels)

########################
## 2. Visualize EEG signals for participants 2 and 13 and spectral densities
########################

#plot time series
plot.ts(eeg_p2$signal[,1:8],main="Participant 2")
plot.ts(eeg_p2$signal[,9:16],main="Participant 2")
plot.ts(eeg_p13$signal[,1:8],main="Participant 13")
plot.ts(eeg_p13$signal[,9:16],main="Participant 13")

#compute local periodogram estimate of spectral matrix and demeaned estimate
pse_p2 <- fhat(eeg_p2$signal,eeg_p2$N,stdz=FALSE)
gpse_p2 <- ghat(pse_p2)
pse_p13 <- fhat(eeg_p13$signal,eeg_p13$N,stdz=FALSE)
gpse_p13 <- ghat(pse_p13)


#plot time series and spectral density for Participant 2 Oz channel (Figure 1)
idxc=which(channels=="OZ")[1]
par(mfrow=c(1,1))
ggplot(data=as.data.frame(eeg_p2$signal[,idxc]),
       mapping=aes(x=seq(0,1,length.out=eeg_p2$Ts),y=`eeg_p2$signal[, idxc]`))+
  geom_line()+
  labs(title='Participant 2 (Oz) Standardized EEG',x='Time',y='Participant 2 (Oz)')+
  hw+scale_x_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

par(mfrow=c(1,1))
image.plot(x=seq(0,1,length.out=eeg_p2$Ts),
           y=eeg_p2$ff[-1]*dsfrq,
           z=t(Re(pse_p2[-1,idxc+(idxc-1)*length(channels),][-1,])),
           axes = TRUE, col = inferno(256),
           main = 'Participant 2 (Oz)',xlab='Time',ylab='Frequency (Hz)',xaxs="i");

par(mfrow=c(1,1))
ggplot(data=as.data.frame(eeg_p13$signal[,idxc]),
       mapping=aes(x=seq(0,1,length.out=eeg_p13$Ts),y=`eeg_p13$signal[, idxc]`))+
  geom_line()+
  labs(title='Participant 13 (Oz) Standardized EEG',x='Time',y='Participant 13 (Oz)')+
  hw+scale_x_continuous(expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

par(mfrow=c(1,1))
image.plot(x=seq(0,1,length.out=eeg_p13$Ts),
           y=eeg_p13$ff[-1]*dsfrq,
           z=t(Re(pse_p13[-1,idxc+(idxc-1)*length(channels),][-1,])),
           axes = TRUE, col = inferno(256),
           main = 'Participant 13 (Oz)',xlab='Time',ylab='Frequency (Hz)',xaxs="i");

########################
## 2. Apply proposed method to signals for participants 2 and 13
########################
#WARNING: Will take hours run the following code on a personal computer 
#due to the length of the time series and number of components.  Instead, 
#you can load the results from the RData file below

# meba_p2 <- msboot(nrep=1000, eeg_p2$signal, Wsel=3, stdz=FALSE,ncore=ncores)
# meba_p13 <- msboot(nrep=1000, eeg_p13$signal, Wsel=3, stdz=FALSE,ncore=ncores)
# save(meba_p2,meba_p13,file='meba_p2p13.RData')

#If you choose to run it, messages are provided to let you
#know code is still running.  First, observed data test
#statistics are calculated. Then, bootstrap draws of test statistics are drawn.
#Then the search for partition points takes place.  This happens for each W (Wsel in total).

load(url('https://www.dropbox.com/s/8b0zuvqb6476iak/meba_p2p13.RData?dl=1'))

########################
## 3. Visualizing results for participants 2 and 13
########################

#Figure 3: estimated spectral densities with frequency band structure for participants 7 and 13
idxc=which(channels=="OZ")[1]
par(mfrow=c(1,1))
image.plot(x=seq(0,1,length.out=eeg_p13$Ts),
           y=eeg_p13$ff[-1]*dsfrq,
           z=t(Re(pse_p13[-1,idxc+(idxc-1)*length(channels),][-1,])),
           axes = TRUE, col = inferno(256),
           main = 'Participant 13 (Oz)',xlab='Time',ylab='Frequency (Hz)',xaxs="i");
abline(h=meba_p13[[4]][which(meba_p13[[4]][,2]==1),1]*64,col='green',lwd=2)

idxc=which(channels=="P7")[1]
par(mfrow=c(1,1))
image.plot(x=seq(0,1,length.out=eeg_p13$Ts),
           y=eeg_p13$ff[-1]*dsfrq,
           z=t(Re(pse_p13[-1,idxc+(idxc-1)*length(channels),][-1,])),
           axes = TRUE, col = inferno(256),
           main = 'Participant 13 (P7)',xlab='Time',ylab='Frequency (Hz)',xaxs="i");
abline(h=meba_p13[[4]][which(meba_p13[[4]][,2]==1),1]*64,col='green',lwd=2)


#Figure 4: components significantly contributed to estimated frequency partition points
#for participants 2 and 13
bfthres=.05/choose(length(channels),2)
xy = expand.grid(channels,channels)
subj.idx=c(2,13)
plts=list()
idx=0

#get indices and frequencies for frequencies identified as partition points
freq.idx=which(meba_p2[[3]][[1]][,1] %in% 
                 meba_p2[[4]][which(meba_p2[[4]][,2]==1),1])
freqs=round(meba_p2[[3]][[1]][freq.idx,1]*64,1) 

for (i in 1:length(freq.idx)){
  idx=idx+1;
  rm(tmpdf)
  rm(p1)
  tmpdf=cbind(xy,meba_p2[[7]][[1]][freq.idx[i],-1]<bfthres)
  p1=ggplot(data=tmpdf)+aes_string(x=tmpdf[,1],y=tmpdf[,2],fill=tmpdf[,3])+
    geom_tile()+hw+ggtitle(paste("Participant 2","Frequency",freqs[i],"Hz",sep=" "))+
    labs(x="",y="")+theme(axis.text.x = element_text(angle = 45,vjust=0.5))+
    scale_fill_manual(values=c("grey90","red"))+theme(legend.position="none")
  # print(p1)
  plts[[idx]]=p1
}

#get indices and frequencies for frequencies identified as partition points
freq.idx=which(meba_p13[[3]][[1]][,1] %in% 
                 meba_p13[[4]][which(meba_p13[[4]][,2]==1),1])
freqs=round(meba_p13[[3]][[1]][freq.idx,1]*64,1) 

for (i in 1:length(freq.idx)){
  idx=idx+1;
  rm(tmpdf)
  rm(p1)
  tmpdf=cbind(xy,meba_p13[[7]][[1]][freq.idx[i],-1]<bfthres)
  p1=ggplot(data=tmpdf)+aes_string(x=tmpdf[,1],y=tmpdf[,2],fill=tmpdf[,3])+
    geom_tile()+hw+ggtitle(paste("Participant 13","Frequency",freqs[i],"Hz",sep=" "))+
    labs(x="",y="")+theme(axis.text.x = element_text(angle = 45,vjust=0.5))+
    scale_fill_manual(values=c("grey90","red"))+theme(legend.position="none")
  # print(p1)
  plts[[idx]]=p1
}
ggarrange(plotlist=plts,ncol=2,nrow=2)

