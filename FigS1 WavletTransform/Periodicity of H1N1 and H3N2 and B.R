library(Rwave)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(gridExtra)
library(boot)
library(moments)
library(mgcv)
library(tidyverse)

read.csv('USA-inc-month-2019.csv') %>%
  select(month,H1N,H3N,B) %>%
  filter(month>=2012.75,month<2019.75) -> inc

inc$H1N <- inc$H1N/(1e+7)
inc$H3N <- inc$H3N/(1e+7)
inc$B <- inc$B/(1e+7)

inc$H1N_H3N_B <- inc$H1N+inc$H3N+inc$B

no <- 8
nv <- 16
a <- 2^seq(1, no+1-1/nv, by=1/nv)

tiff('periodicity of H1 H3 B.tiff', height = 2000, width = 2500, res=250)
par(mfrow=c(7,1), mar=c(0.5,7,0.5,0.75), oma=c(1.5,0.5,1.5,0.5))

xlim <- c(2012.75,2019.75)
v <- seq(2013,2019,by=1)
at <- 2013:2019
labels <- at

##H1N1--------------------------------------------------------------------------------
plot(inc$month,inc$H1N, type="b", xaxt='n', xlab="Time(monthly)", 
     ylab=expression(paste('Incidence'%*%10^7)), 
     col="blue", xaxs="i", xlim=xlim, ylim=c(0,1.5))
axis(side=1, at=at, labels=F)
axis(side=3, at=at, labels=labels)
abline(v=v, col="grey80", lty=2)
text(2019,1.25,labels = 'H1N1', col='blue')

wfit <- cwt(sqrt(inc$H1N), no, nv, plot=F)
wspec <- Mod(wfit)
image(inc$month, wspec, y=a/12, ylim=c(0,5), xaxt="n", xlab="", ylab="Period(years)", 
      main="", xlim=xlim)
contour(inc$month, wspec, y=a/12, ylim=c(0,5), add=T, drawlabels=F)
axis(side=1, at=at, labels=F)
abline(v=v, col="grey80", lty=2)

##H3N2--------------------------------------------------------------------------------
plot(inc$month,inc$H3N, type="b", xaxt='n', xlab="", 
     ylab=expression(paste('Incidence'%*%10^7)), 
     col="red", xaxs="i", xlim=xlim, ylim=c(0,2))
axis(side=1, at=at, labels=F)
abline(v=v, col="grey80", lty=2)
text(2019,1.75,labels = 'H3N2', col='red')

wfit <- cwt(inc$H3N, no, nv, plot=F)
wspec <- Mod(wfit)
image(inc$month, wspec, y=a/12, ylim=c(0,5), xaxt="n", xlab="", ylab="Period(years)", 
      main="", xlim=xlim)
contour(inc$month, wspec, y=a/12, ylim=c(0,5), add=T, drawlabels=F)
axis(side=1, at=at, labels=F)
abline(v=v, col="grey80", lty=2)


##B--------------------------------------------------------------------------------
plot(inc$month,inc$B, type="b", xaxt='n', xlab="", 
     ylab=expression(paste('Incidence'%*%10^7)), 
     col="green", xaxs="i", xlim=xlim, ylim=c(0,2))
axis(side=1, at=at, labels=F)
abline(v=v, col="grey80", lty=2)
text(2019,1.75,labels = 'B', col='green')

wfit <- cwt(inc$B, no, nv, plot=F)
wspec <- Mod(wfit)
image(inc$month, wspec, y=a/12, ylim=c(0,5), xaxt="n", xlab="", ylab="Period(years)", 
      main="", xlim=xlim)
contour(inc$month, wspec, y=a/12, ylim=c(0,5), add=T, drawlabels=F)
axis(side=1, at=at, labels=F)
abline(v=v, col="grey80", lty=2)


##H1N1+H3N2+B-----------------------------------------------------------------------------
plot(inc$month,inc$H1N_H3N_B, type="b", xaxt='n', xlab="", 
     ylab=expression(paste('Incidence'%*%10^7)), 
     col="black", xaxs="i", xlim=xlim, ylim=c(0,3))
axis(side=1, at=at, labels=F)
abline(v=v, col="grey80", lty=2)
text(2019,2,labels = 'H1N1+H3N2+B', col='black')


dev.off()