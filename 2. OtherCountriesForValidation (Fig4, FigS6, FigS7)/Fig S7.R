library(ggplot2)
library(tidyverse)
library(gridExtra)

sim1_H1N_summary <- read.csv('Germany sim1_H1N_summary.csv')
sim1_H3N_summary <- read.csv('Germany sim1_H3N_summary.csv')
sim1_H1N_summary <- na.omit(sim1_H1N_summary)
sim1_H3N_summary <- na.omit(sim1_H3N_summary)

cor.H1N <- cor(sim1_H1N_summary$mean, sim1_H1N_summary$H1N)
cortest.H1N <- cor.test(sim1_H1N_summary$mean, sim1_H1N_summary$H1N)
cor.H3N <- cor(sim1_H3N_summary$mean, sim1_H3N_summary$H3N)
cortest.H3N <- cor.test(sim1_H3N_summary$mean, sim1_H3N_summary$H3N)

label1 <- sprintf('r=%.2f',cor.H1N)
f1 <- ggplot(data=sim1_H1N_summary, aes(x=H1N, y=mean, ymin=lower, ymax=upper)) + 
  geom_linerange(color='LightBlue', size=0.5) + 
  geom_point(color='blue',size=0.5) +
  geom_abline(intercept = 0,slope=1, color='black', linetype="dotted") +
  labs(x=expression(paste('Observed incidence'%*%10^6)),
       y=expression(paste('Simulated incidence'%*%10^6)),
       title=label1) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 5), axis.title.y = element_text(size = 5)) +
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),title=element_text(size=5))

label2 <- sprintf('r=%.2f',cor.H3N)
f2 <- ggplot(data=sim1_H3N_summary, aes(x=H3N, y=mean, ymin=lower, ymax=upper)) + 
  geom_linerange(color='LightCoral', size=0.5) + 
  geom_point(color='red',size=0.5) +
  geom_abline(intercept = 0,slope=1, color='black', linetype="dotted") +
  labs(x=expression(paste('Observed incidence'%*%10^6)),
       y=expression(paste('Simulated incidence'%*%10^6)),
       title=label2) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 5), axis.title.y = element_text(size = 5)) +
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),title=element_text(size=5))

g1 <- grid.arrange(f1, f2, nrow=1)


sim1_H1N_summary <- read.csv('Austria sim1_H1N_summary.csv')
sim1_H3N_summary <- read.csv('Austria sim1_H3N_summary.csv')
sim1_H1N_summary <- na.omit(sim1_H1N_summary)
sim1_H3N_summary <- na.omit(sim1_H3N_summary)

cor.H1N <- cor(sim1_H1N_summary$mean, sim1_H1N_summary$H1N)
cortest.H1N <- cor.test(sim1_H1N_summary$mean, sim1_H1N_summary$H1N)
cor.H3N <- cor(sim1_H3N_summary$mean, sim1_H3N_summary$H3N)
cortest.H3N <- cor.test(sim1_H3N_summary$mean, sim1_H3N_summary$H3N)

label1 <- sprintf('r=%.2f',cor.H1N)
f1 <- ggplot(data=sim1_H1N_summary, aes(x=H1N, y=mean, ymin=lower, ymax=upper)) + 
  geom_linerange(color='LightBlue', size=0.5) + 
  geom_point(color='blue',size=0.5) +
  geom_abline(intercept = 0,slope=1, color='black', linetype="dotted") +
  labs(x=expression(paste('Observed incidence'%*%10^6)),
       y=expression(paste('Simulated incidence'%*%10^6)),
       title=label1) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 5), axis.title.y = element_text(size = 5)) +
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),title=element_text(size=5))

label2 <- sprintf('r=%.2f',cor.H3N)
f2 <- ggplot(data=sim1_H3N_summary, aes(x=H3N, y=mean, ymin=lower, ymax=upper)) + 
  geom_linerange(color='LightCoral', size=0.5) + 
  geom_point(color='red',size=0.5) +
  geom_abline(intercept = 0,slope=1, color='black', linetype="dotted") +
  labs(x=expression(paste('Observed incidence'%*%10^6)),
       y=expression(paste('Simulated incidence'%*%10^6)),
       title=label2) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 5), axis.title.y = element_text(size = 5)) +
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),title=element_text(size=5))

g2 <- grid.arrange(f1, f2, nrow=1)


sim1_H1N_summary <- read.csv('Norway sim1_H1N_summary.csv')
sim1_H3N_summary <- read.csv('Norway sim1_H3N_summary.csv')
sim1_H1N_summary <- na.omit(sim1_H1N_summary)
sim1_H3N_summary <- na.omit(sim1_H3N_summary)

cor.H1N <- cor(sim1_H1N_summary$mean, sim1_H1N_summary$H1N)
cortest.H1N <- cor.test(sim1_H1N_summary$mean, sim1_H1N_summary$H1N)
cor.H3N <- cor(sim1_H3N_summary$mean, sim1_H3N_summary$H3N)
cortest.H3N <- cor.test(sim1_H3N_summary$mean, sim1_H3N_summary$H3N)

label1 <- sprintf('r=%.2f',cor.H1N)
f1 <- ggplot(data=sim1_H1N_summary, aes(x=H1N, y=mean, ymin=lower, ymax=upper)) + 
  geom_linerange(color='LightBlue', size=0.5) + 
  geom_point(color='blue',size=0.5) +
  geom_abline(intercept = 0,slope=1, color='black', linetype="dotted") +
  labs(x=expression(paste('Observed incidence'%*%10^6)),
       y=expression(paste('Simulated incidence'%*%10^6)),
       title=label1) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 5), axis.title.y = element_text(size = 5)) +
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),title=element_text(size=5))

label2 <- sprintf('r=%.2f',cor.H3N)
f2 <- ggplot(data=sim1_H3N_summary, aes(x=H3N, y=mean, ymin=lower, ymax=upper)) + 
  geom_linerange(color='LightCoral', size=0.5) + 
  geom_point(color='red',size=0.5) +
  geom_abline(intercept = 0,slope=1, color='black', linetype="dotted") +
  labs(x=expression(paste('Observed incidence'%*%10^6)),
       y=expression(paste('Simulated incidence'%*%10^6)),
       title=label2) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 5), axis.title.y = element_text(size = 5)) +
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),title=element_text(size=5))

g3 <- grid.arrange(f1, f2, nrow=1)

tiff('FigS7.tiff', height = 1500, width = 1500, res = 350, units = 'px',pointsize = 6,compression = 'jpeg')
ggarrange(g1,g2,g3, labels = c('A','B','C'), nrow=3,font.label = list(size=10))
dev.off()

