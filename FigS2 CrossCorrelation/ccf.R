remove(list = ls())
library(tidyverse)

read.csv('USA-inc-month-2019.csv') %>%
  filter(month>=2012.75, month<2019.75) -> inc

at <- -16:16

tiff('ccf.tiff', height = 3000, width = 2000, res = 250)
par(mfrow=c(3,1))
ccfH1NH3N <- ccf(inc$H1N, inc$H3N, main='H1N1 and H3N2', ylab = 'cross-correlation',xlab='Lag (months)')
ccfH1NB <- ccf(inc$H1N, inc$B, main='H1N1 and B', ylab = 'cross-correlation', xlab='Lag (months)')
ccfH3NB <- ccf(inc$H3N, inc$B, main='H3N2 and B', ylab = 'cross-correlation', xlab='Lag (months)')
dev.off()