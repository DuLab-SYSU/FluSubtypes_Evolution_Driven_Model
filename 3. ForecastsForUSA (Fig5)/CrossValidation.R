library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)

##obs Vs pre--------------------------------------------------------------------------------------------------
cv2016 <- read.csv('CrossValidation 2016.csv')
cv2017 <- read.csv('CrossValidation 2017.csv')
cv2018 <- read.csv('CrossValidation 2018.csv')
cv2019 <- read.csv('CrossValidation 2019.csv')
cv2020 <- read.csv('CrossValidation 2020.csv')

##cv2016----------------------------------------------------------------------------------
cv2016H1N <- ggplot(cv2016,aes(x = month))+
  geom_ribbon(aes(y = sim.history.H1N.mean,
                  ymin = sim.history.H1N.lower,
                  ymax = sim.history.H1N.upper),alpha = 0.2,fill = "#0000FF") + 
  geom_line(aes(y = H1N/1e+07),size = 0.8,color = "black",alpha =0.8) +
  geom_line(aes(y = sim.history.H1N.mean),color = "blue",size = 0.8,alpha =0.9) +
  geom_ribbon(alpha = 0.2,map=aes(y = Pre.2016.H1N.mean,
                                  ymin = Pre.2016.H1N.lower,
                                  ymax = Pre.2016.H1N.upper),fill="DarkGreen") + 
  geom_line(aes(y = Pre.2016.H1N.mean),color = "DarkGreen",size = 0.8,alpha =0.9) +
  ylab(expression(paste('Incidence'%*%10^7))) + xlab('Time(monthly)') +
  scale_x_continuous(breaks = seq(2013,2020,2)) +
  scale_y_continuous(limits = c(0,6.5))+
  annotate('text', x=2019, y=5, label='H1N1', color='blue', alpha=0.8, size=1) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))


cv2016H3N <- ggplot(cv2016,aes(x = month))+
  geom_ribbon(aes(y = sim.history.H3N.mean,
                  ymin = sim.history.H3N.lower,
                  ymax = sim.history.H3N.upper),alpha = 0.2,fill = "#EE6363") + 
  geom_line(aes(y = H3N/1e+07),size = 0.8,color = "black",alpha =0.8) +
  geom_line(aes(y = sim.history.H3N.mean),color = "red",size = 0.8,alpha =0.9) +
  geom_ribbon(alpha = 0.2,map=aes(y = Pre.2016.H3N.mean,
                                  ymin = Pre.2016.H3N.lower,
                                  ymax = Pre.2016.H3N.upper),fill="DarkGreen") + 
  geom_line(aes(y = Pre.2016.H3N.mean),color = "DarkGreen",size = 0.8,alpha =0.9) +
  ylab('') + xlab('Time(monthly)') +
  scale_x_continuous(breaks = seq(2013,2020,2)) +
  scale_y_continuous(limits = c(0,6.5))+
  annotate('text', x=2019, y=5, label='H3N2', color='red', alpha=0.8, size=1) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))

g2016 <- grid.arrange(cv2016H1N, cv2016H3N, nrow=1)


##cv2017----------------------------------------------------------------------------------
cv2017H1N <- ggplot(cv2017,aes(x = month))+
  geom_ribbon(aes(y = sim.history.H1N.mean,
                  ymin = sim.history.H1N.lower,
                  ymax = sim.history.H1N.upper),alpha = 0.2,fill = "#0000FF") + 
  geom_line(aes(y = H1N/1e+07),size = 0.5,color = "black",alpha =0.8) +
  geom_line(aes(y = sim.history.H1N.mean),color = "blue",size = 0.5,alpha =0.9) +
  geom_ribbon(alpha = 0.2,map=aes(y = Pre.2017.H1N.mean,
                                  ymin = Pre.2017.H1N.lower,
                                  ymax = Pre.2017.H1N.upper),fill="DarkGreen") + 
  geom_line(aes(y = Pre.2017.H1N.mean),color = "DarkGreen",size = 0.5,alpha =0.9) +
  ylab(expression(paste('Incidence'%*%10^7))) + xlab('Time(monthly)') +
  scale_x_continuous(breaks = seq(2013,2020,2)) +
  scale_y_continuous(limits = c(0,6.5))+
  annotate('text', x=2019, y=5, label='H1N1', color='blue', alpha=0.8, size=2) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))


cv2017H3N <- ggplot(cv2017,aes(x = month))+
  geom_ribbon(aes(y = sim.history.H3N.mean,
                  ymin = sim.history.H3N.lower,
                  ymax = sim.history.H3N.upper),alpha = 0.2,fill = "#EE6363") + 
  geom_line(aes(y = H3N/1e+07),size = 0.5,color = "black",alpha =0.8) +
  geom_line(aes(y = sim.history.H3N.mean),color = "red",size = 0.5,alpha =0.9) +
  geom_ribbon(alpha = 0.2,map=aes(y = Pre.2017.H3N.mean,
                                  ymin = Pre.2017.H3N.lower,
                                  ymax = Pre.2017.H3N.upper),fill="DarkGreen") + 
  geom_line(aes(y = Pre.2017.H3N.mean),color = "DarkGreen",size = 0.5,alpha =0.9) +
  ylab('') + xlab('Time(monthly)') +
  scale_x_continuous(breaks = seq(2013,2020,2)) +
  scale_y_continuous(limits = c(0,6.5))+
  annotate('text', x=2019, y=5, label='H3N2', color='red', alpha=0.8, size=2) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))

g2017 <- grid.arrange(cv2017H1N, cv2017H3N, nrow=1)

##cv2018----------------------------------------------------------------------------------
cv2018H1N <- ggplot(cv2018,aes(x = month))+
  geom_ribbon(aes(y = sim.history.H1N.mean,
                  ymin = sim.history.H1N.lower,
                  ymax = sim.history.H1N.upper),alpha = 0.2,fill = "#0000FF") + 
  geom_line(aes(y = H1N/1e+07),size = 0.5,color = "black",alpha =0.8) +
  geom_line(aes(y = sim.history.H1N.mean),color = "blue",size = 0.5,alpha =0.9) +
  geom_ribbon(alpha = 0.2,map=aes(y = Pre.2018.H1N.mean,
                                  ymin = Pre.2018.H1N.lower,
                                  ymax = Pre.2018.H1N.upper),fill="DarkGreen") + 
  geom_line(aes(y = Pre.2018.H1N.mean),color = "DarkGreen",size = 0.5,alpha =0.9) +
  ylab(expression(paste('Incidence'%*%10^7))) + xlab('Time(monthly)') +
  scale_x_continuous(breaks = seq(2013,2020,2)) +
  scale_y_continuous(limits = c(0,6.5))+
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))


cv2018H3N <- ggplot(cv2018,aes(x = month))+
  geom_ribbon(aes(y = sim.history.H3N.mean,
                  ymin = sim.history.H3N.lower,
                  ymax = sim.history.H3N.upper),alpha = 0.2,fill = "#EE6363") + 
  geom_line(aes(y = H3N/1e+07),size = 0.5,color = "black",alpha =0.8) +
  geom_line(aes(y = sim.history.H3N.mean),color = "red",size = 0.5,alpha =0.9) +
  geom_ribbon(alpha = 0.2,map=aes(y = Pre.2018.H3N.mean,
                                  ymin = Pre.2018.H3N.lower,
                                  ymax = Pre.2018.H3N.upper),fill="DarkGreen") + 
  geom_line(aes(y = Pre.2018.H3N.mean),color = "DarkGreen",size = 0.5,alpha =0.9) +
  ylab('') + xlab('Time(monthly)') +
  scale_x_continuous(breaks = seq(2013,2020,2)) +
  scale_y_continuous(limits = c(0,6.5))+
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))

g2018 <- grid.arrange(cv2018H1N, cv2018H3N, nrow=1)
##cv2019--------------------------------------------------------------------------------------------
cv2019H1N <- ggplot(cv2019,aes(x = month))+
  geom_ribbon(aes(y = sim.history.H1N.mean,
                  ymin = sim.history.H1N.lower,
                  ymax = sim.history.H1N.upper),alpha = 0.2,fill = "#0000FF") + 
  geom_line(aes(y = H1N/1e+07),size = 0.5,color = "black",alpha =0.8) +
  geom_line(aes(y = sim.history.H1N.mean),color = "blue",size = 0.5,alpha =0.9) +
  geom_ribbon(alpha = 0.2,map=aes(y = Pre.2019.H1N.mean,
                                  ymin = Pre.2019.H1N.lower,
                                  ymax = Pre.2019.H1N.upper),fill="DarkGreen") + 
  geom_line(aes(y = Pre.2019.H1N.mean),color = "DarkGreen",size = 0.5,alpha =0.9) +
  ylab(expression(paste('Incidence'%*%10^7))) + xlab('Time(monthly)') +
  scale_x_continuous(breaks = seq(2013,2020,2)) +
  scale_y_continuous(limits = c(0,6.5))+
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))

cv2019H3N <- ggplot(cv2019,aes(x = month))+
  geom_ribbon(aes(y = sim.history.H3N.mean,
                  ymin = sim.history.H3N.lower,
                  ymax = sim.history.H3N.upper),alpha = 0.2,fill = "#EE6363") + 
  geom_line(aes(y = H3N/1e+07),size = 0.5,color = "black",alpha =0.8) +
  geom_line(aes(y = sim.history.H3N.mean),color = "red",size = 0.5,alpha =0.9) +
  geom_ribbon(alpha = 0.2,map=aes(y = Pre.2019.H3N.mean,
                                  ymin = Pre.2019.H3N.lower,
                                  ymax = Pre.2019.H3N.upper),fill="DarkGreen") + 
  geom_line(aes(y = Pre.2019.H3N.mean),color = "DarkGreen",size = 0.5,alpha =0.9) +
  ylab('') + xlab('Time(monthly)') +
  scale_x_continuous(breaks = seq(2013,2020,2)) +
  scale_y_continuous(limits = c(0,6.5))+
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))

g2019 <- grid.arrange(cv2019H1N, cv2019H3N, nrow=1)

##cv2020----------------------------------------------------------------------------------------
cv2020H1N <- ggplot(cv2020,aes(x = month))+
  geom_ribbon(aes(y = sim.history.H1N.mean,
                  ymin = sim.history.H1N.lower,
                  ymax = sim.history.H1N.upper),alpha = 0.2,fill = "#0000FF") + 
  geom_line(aes(y = H1N/1e+07),size = 0.5,color = "black",alpha =0.8) +
  geom_line(aes(y = sim.history.H1N.mean),color = "blue",size = 0.5,alpha =0.9) +
  geom_ribbon(alpha = 0.2,map=aes(y = Pre.2020.H1N.mean,
                                  ymin = Pre.2020.H1N.lower,
                                  ymax = Pre.2020.H1N.upper),fill="orange") + 
  geom_line(aes(y = Pre.2020.H1N.mean),color = "DarkOrange",size = 0.5,alpha =0.9) +
  ylab(expression(paste('Incidence'%*%10^7))) + xlab('Time(monthly)') +
  scale_x_continuous(breaks = seq(2013,2020,2)) +
  scale_y_continuous(limits = c(0,6.5))+
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))

cv2020H3N <- ggplot(cv2020,aes(x = month))+
  geom_ribbon(aes(y = sim.history.H3N.mean,
                  ymin = sim.history.H3N.lower,
                  ymax = sim.history.H3N.upper),alpha = 0.2,fill = "#EE6363") + 
  geom_line(aes(y = H3N/1e+07),size = 0.5,color = "black",alpha =0.8) +
  geom_line(aes(y = sim.history.H3N.mean),color = "red",size = 0.5,alpha =0.9) +
  geom_ribbon(alpha = 0.2,map=aes(y = Pre.2020.H3N.mean,
                                  ymin = Pre.2020.H3N.lower,
                                  ymax = Pre.2020.H3N.upper),fill="orange") + 
  geom_line(aes(y = Pre.2020.H3N.mean),color = "DarkOrange",size = 0.5,alpha =0.9) +
  ylab('') + xlab('Time(monthly)') +
  scale_x_continuous(breaks = seq(2013,2020,2)) +
  scale_y_continuous(limits = c(0,6.5))+
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))

g2020 <- grid.arrange(cv2020H1N, cv2020H3N, nrow=1)

library(ggpubr)
tiff('CV.tiff', height = 2000, width = 1500, res = 250)
ggarrange(g2017, g2018, g2019, g2020, nrow=4, labels = c('A','B','C','D'))
dev.off()
g.CV <- ggarrange(g2017, g2018, g2019, g2020, nrow=4, ncol=1)

##corr obs and pre-------------------------------------------------------------------------------------------
dat <- read.csv('CrossValidation 2017 2018 2019 2020.csv')

cor.dat <- read.csv('CrossValidation 2017 2018 2019 2020 corr.csv')
cor.dat <- cor.dat[which(cor.dat$H1N != 'NA'),]
cor.dat <- cor.dat[which(cor.dat$sim.H1N.mean != 'NA'),]
cor.dat <- cor.dat[which(cor.dat$sim.H1N.mean != 0),]

cor.H1N <- cor(cor.dat$H1N/(1e+07),cor.dat$sim.H1N.mean)
cor.H3N <- cor(cor.dat$H3N/(1e+07),cor.dat$sim.H3N.mean)
cor.test(cor.dat$H1N/(1e+07),cor.dat$sim.H1N.mean)
cor.test(cor.dat$H3N/(1e+07),cor.dat$sim.H3N.mean)

#label1 <- sprintf('r=%.2f, P<0.01', cor.H1N)
label1 <- expression(paste('r=0.85,',italic('P'), '<0.01'))

cv.H1N <- ggplot(dat, aes(x=H1N/(1e+07))) + 
  #geom_point(aes(y=sim.history.H1N.mean),color='blue') +
  geom_point(aes(y=Pre.H1N.mean),color='DarkGreen') +
  geom_errorbar(aes(ymin=Pre.H1N.lower, ymax=Pre.H1N.upper)) +
  geom_point(aes(y=Pre.H1N.mean19_20),color='DarkOrange') +
  geom_errorbar(aes(ymin=Pre.H1N.lower19_20, ymax=Pre.H1N.upper19_20)) +
  geom_abline(intercept=0,slope=1,linetype="dotted")+
  scale_x_continuous(limits=c(0, 1.5),breaks=seq(0, 1.5, 0.3))+
  scale_y_continuous(limits=c(0, 2.5),breaks=seq(0, 2.5, 0.3))+
  labs(x=expression(paste('Observed incidence'%*%10^7)),
       y=expression(paste('Predicted incidence'%*%10^7)))+
  annotate("text", x = 0.3, y = 2.3, label = label1, size=2)+
  annotate("text", x = 0.3, y = 2.45, label = 'H1N1', color='blue', size=2)+
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))

dat.tem.H1N <- dat[which(dat$H1N/(1e+07)<=0.6),]
cv.H1N.vp <- ggplot(dat.tem.H1N, aes(x=H1N/(1e+07))) + 
  geom_point(aes(y=sim.history.H1N.mean),color='blue') +
  geom_point(aes(y=Pre.H1N.mean),color='purple') +
  geom_abline(intercept=0,slope=1,linetype="dotted")+
  scale_x_continuous(limits=c(0, 0.6),breaks=seq(0, 0.6, 0.2))+
  scale_y_continuous(limits=c(0, 0.6),breaks=seq(0, 0.6, 0.2))+
  ylab('') + xlab('') +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

vp <- viewport(width = 0.4, height = 0.4, x = 0.125,y = 0.9,just=c("left","top"))
print(cv.H1N)
print(cv.H1N.vp,vp=vp)

#label2 <- sprintf('r=%.2f, P<0.01', cor.H3N)
label2 <- expression(paste('r=0.54,',italic('P'), '<0.01'))

cv.H3N <- ggplot(dat, aes(x=H3N/(1e+07))) + 
  #geom_point(aes(y=sim.history.H3N.mean),color='red') +
  geom_point(aes(y=Pre.H3N.mean),color='DarkGreen') +
  geom_errorbar(aes(ymin=Pre.H3N.lower, ymax=Pre.H3N.upper)) +
  geom_point(aes(y=Pre.H3N.mean19_20),color='DarkOrange') +
  geom_errorbar(aes(ymin=Pre.H3N.lower19_20, ymax=Pre.H3N.upper19_20)) +
  geom_abline(intercept=0,slope=1,linetype="dotted")+
  scale_x_continuous(limits=c(0, 1.5),breaks=seq(0, 1.5, 0.3))+
  scale_y_continuous(limits=c(0, 5),breaks=seq(0, 5, 0.5))+
  labs(x=expression(paste('Observed incidence'%*%10^7)),
       y=expression(paste('Predicted incidence'%*%10^7)))+
  annotate("text", x = 0.3, y = 4.65, label = label2, size=2)+
  annotate("text", x = 0.3, y = 4.95, label = 'H3N2', color='red', size=2)+
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))

dat.tem.H3N <- dat[which(dat$H3N/(1e+07)<=5),]
cv.H3N.vp <- ggplot(dat.tem.H3N, aes(x=H1N/(1e+07))) + 
  #geom_point(aes(y=sim.history.H3N.mean),color='red') +
  geom_point(aes(y=Pre.H3N.mean),color='DarkGreen') +
  geom_errorbar(aes(ymin=Pre.H3N.lower, ymax=Pre.H3N.upper)) +
  geom_point(aes(y=Pre.H3N.mean19_20),color='DarkOrange') +
  geom_errorbar(aes(ymin=Pre.H3N.lower19_20, ymax=Pre.H3N.upper19_20)) +
  geom_abline(intercept=0,slope=1,linetype="dotted")+
  scale_x_continuous(limits=c(0, 0.4),breaks=seq(0, 0.4, 0.05))+
  scale_y_continuous(limits=c(0, 5),breaks=seq(0, 5, 1))+
  ylab('') + xlab('') +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

vp <- viewport(width = 0.4, height = 0.4, x = 0.9,y = 0.9,just=c("right","top"))
print(cv.H3N)
print(cv.H3N.vp,vp=vp)

tiff('CV corr.tiff', height = 1000, width = 2000, res = 250)
ggarrange(cv.H1N, cv.H3N, nrow=1, labels=c('A','B'))
dev.off()

g.CV.cor <- ggarrange(cv.H1N, cv.H3N, nrow=2)

tiff('Fig5.tiff', width = 2200, height = 2000, res = 350, units = 'px',pointsize = 12, compression = 'jpeg')
ggdraw() +
  draw_plot(g.CV, 0, 0, width = 3/5, height = 1) +
  draw_plot(g.CV.cor, 3/5, 0, width = 2/5, height = 1) +
  draw_plot_label(label = c('A','B'), x=c(0,3/5), y=c(1,1),size = 10)
dev.off()
