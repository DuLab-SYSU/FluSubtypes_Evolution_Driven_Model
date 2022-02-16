library(ggplot2)
library(tidyverse)
library(gridExtra)
##Germany---------------------------------------------------------------------------------------------------------------------
sim1_H1N_summary <- read.csv('Germany sim1_H1N_summary.csv')
sim1_H3N_summary <- read.csv('Germany sim1_H3N_summary.csv')

figure_H1N <- ggplot(sim1_H1N_summary,aes(x = month,y = mean,ymin = lower,ymax = upper))+
  geom_ribbon(alpha = 0.2,fill = "#0000FF") + 
  geom_line(aes(y = H1N),size = 0.8,color = "black",alpha =0.8) +
  geom_line(aes(y = mean),color = "blue",size = 0.8,alpha =0.9) +
  labs(y=expression(paste('Incidence'%*%10^6)), x='', title='Germany')+
  annotate("text", x = 2019, y = 4, label = 'H1N1', color='blue', size=2)+
  scale_y_continuous(limits = c(0,4.5))+
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5), title=element_text(size=7.5))

figure_H3N <- ggplot(sim1_H3N_summary,aes(x = month,y = mean,ymin = lower,ymax = upper))+
  geom_ribbon(alpha = 0.2,fill = "#EE6363") + 
  geom_line(aes(y = H3N),size = 0.8,color = "black",alpha =0.8) +
  geom_line(aes(y = mean),color = "red",size = 0.8,alpha =0.9) +
  labs(y='',x='',title = '')+
  annotate("text", x = 2019, y = 4, label = 'H3N2', color='red', size=2)+
  scale_y_continuous(limits = c(0,4.5))+
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))

Et <- read.csv('GermanyEt.csv')
figure_Et_Germany <- ggplot(Et, aes(x=month)) +
  geom_line(aes(y=Et_H1N), color='blue') +
  geom_line(aes(y=Et_H3N), color='red') +
  labs(x='', y='Evolutionary index',title = '') +
  scale_x_continuous(breaks = seq(2013,2019,1)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))

g1 <- grid.arrange(figure_H1N,figure_H3N,nrow=1)

##Austria--------------------------------------------------------------------------------------------------------------------------------------
sim1_H1N_summary <- read.csv('Austria sim1_H1N_summary.csv')
sim1_H3N_summary <- read.csv('Austria sim1_H3N_summary.csv')

figure_H1N <- ggplot(sim1_H1N_summary,aes(x = month,y = mean,ymin = lower,ymax = upper))+
  geom_ribbon(alpha = 0.2,fill = "#0000FF") + 
  geom_line(aes(y = H1N),size = 0.8,color = "black",alpha =0.8) +
  geom_line(aes(y = mean),color = "blue",size = 0.8,alpha =0.9) +
  labs(y=expression(paste('Incidence'%*%10^6)),x='',  title='Austria')+
  #annotate("text", x = 2019, y = 4, label = 'H1N1', color='blue')+
  scale_y_continuous(limits = c(0,4.5))+
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),title=element_text(size=7.5))

figure_H3N <- ggplot(sim1_H3N_summary,aes(x = month,y = mean,ymin = lower,ymax = upper))+
  geom_ribbon(alpha = 0.2,fill = "#EE6363") + 
  geom_line(aes(y = H3N),size = 0.8,color = "black",alpha =0.8) +
  geom_line(aes(y = mean),color = "red",size = 0.8,alpha =0.9) +
  labs(y='',x='',title = '')+
  #annotate("text", x = 2019, y = 4, label = 'H3N2', color='red')+
  scale_y_continuous(limits = c(0,4.5))+
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))

Et <- read.csv('AustriaEt.csv')
figure_Et_Austria <- ggplot(Et, aes(x=month)) +
  geom_line(aes(y=Et_H1N), color='blue') +
  geom_line(aes(y=Et_H3N), color='red') +
  labs(x='', y='Evolutionary index',title = '') +
  scale_x_continuous(breaks = seq(2013,2019,1)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))

g2 <- grid.arrange(figure_H1N,figure_H3N,nrow=1)


##Norway--------------------------------------------------------------------------------------------------------------------------------------
sim1_H1N_summary <- read.csv('Norway sim1_H1N_summary.csv')
sim1_H3N_summary <- read.csv('Norway sim1_H3N_summary.csv')

figure_H1N <- ggplot(sim1_H1N_summary,aes(x = month,y = mean,ymin = lower,ymax = upper))+
  geom_ribbon(alpha = 0.2,fill = "#0000FF") + 
  geom_line(aes(y = H1N),size = 0.8,color = "black",alpha =0.8) +
  geom_line(aes(y = mean),color = "blue",size = 0.8,alpha =0.9) +
  labs(y=expression(paste('Incidence'%*%10^4)),x="Time(monthly)", title='Norway')+
  #annotate("text", x = 2019, y = 1, label = 'H1N1', color='blue')+
  scale_y_continuous(limits = c(0,1.5))+
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),title=element_text(size=7.5))


figure_H3N <- ggplot(sim1_H3N_summary,aes(x = month,y = mean,ymin = lower,ymax = upper))+
  geom_ribbon(alpha = 0.2,fill = "#EE6363") + 
  geom_line(aes(y = H3N),size = 0.8,color = "black",alpha =0.8) +
  geom_line(aes(y = mean),color = "red",size = 0.8,alpha =0.9) +
  labs(y='',x="Time(monthly)",title = '')+
  #annotate("text", x = 2019, y = 1, label = 'H3N2', color='red')+
  scale_y_continuous(limits = c(0,1.5))+
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))

Et <- read.csv('NorwayEt.csv')

figure_Et_Norway <- ggplot(Et, aes(x=month)) +
  geom_line(aes(y=Et_H1N), color='blue') +
  geom_line(aes(y=Et_H3N), color='red') +
  labs(x='Time(monthly)', y='Evolutionary index',title = '') +
  scale_x_continuous(breaks = seq(2013,2019,1)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))

g3 <- grid.arrange(figure_H1N,figure_H3N,nrow=1)


library(cowplot)
tiff('Fig3.tiff', height = 1750, width = 2000, res = 300, units = 'px', pointsize = 12,compression = 'jpeg')
ggdraw() +
  draw_plot(g1, 0, 2/3, width = 0.65, height = 1/3) +
  draw_plot(figure_Et_Germany, 0.65, 2/3, width = 0.35, height = 1/3) +
  draw_plot(g2, 0, 1/3, width = 0.65, height = 1/3) +
  draw_plot(figure_Et_Austria, 0.65, 1/3, width = 0.35, height = 1/3) +
  draw_plot(g3, 0, 0, width = 0.65, height = 1/3) +
  draw_plot(figure_Et_Norway, 0.65, 0, width = 0.35, height = 1/3) +
  draw_plot_label(label = c('A','B','C'), x=c(0,0,0), y=c(1,2/3,1/3), size = 10) 
dev.off()

