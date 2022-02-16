remove(list=ls())

library(tidyverse)
library(pomp)
library(dplyr)
library(foreach)
library(doParallel)
library(iterators)
library(doRNG)
library(gridExtra)

registerDoParallel()

estimates <- read.csv('profile c1')
estimates <- arrange(estimates, -loglik)
tem.ests <- estimates[c(1:5),]
tem.ests %>%
  filter(loglik>=max(loglik)-5) %>%
  filter(loglik.se<2) %>%
  select(-loglik,-loglik.se) %>%
  apply(2,quantile) %>%
  as.data.frame() -> PredefinedEstimates

best <- PredefinedEstimates[1,]

expand.grid(c1=rep(seq(from=0,to=2,length=40),each=3), 
            c2=rep(seq(from=0,to=2,length=40),each=3)) -> params

read.csv('sim c1=1 c2=1.csv') %>%
  filter(month>=2012.75) %>%
  filter(month<2019.75) %>% 
  select(month,H1N,H3N)  -> inc
inc$H1N <- round(inc$H1N*(1e+7),0)
inc$H3N <- round(inc$H3N*(1e+7),0)


ggplot(inc,aes(x=month))+
  geom_line(aes(y=H1N),color='blue')+
  geom_line(aes(y=H3N),color='red')
  
cov <- read.csv('USA-covar-2019.csv')
cov <- cov[which(cov$month>=2012.75),]
cov <- cov[which(cov$month<2019.75),]

##pomp------------------------------------------------------------------------------
rproc <- Csnippet('
                  double mu = 0.015;
                  double dW = rgammawn(sigma,dt);
                  double A = 0.1*365;
                  double beta = exp(b1*month1+b2*month2+b3*month3+b4*month4+b5*month5+b6*month6)*(dW/dt);
                  double Et1 = y11*exp(-1/theta)+y12*exp(-2/theta);//H3
                  double Et2 = y21*exp(-1/theta)+y22*exp(-2/theta);//H1

                  double dASI = A*dt;
                  double dBS = (pop+dpop)*mu*dt;
                  double dSI1 = beta*S*pow((xi*I1+J1)/pop,alpha)*dt;
                  double dSI2 = beta*S*pow((xi*I2+J2)/pop,alpha)*dt;
                  
                  double dI1Rc1 = I1/D*dt;
                  double dI2Rc2 = I2/D*dt;
                  
                  double dRc1R1 = Rc1/P1*dt;
                  double dRc2R2 = Rc2/P2*dt;

                  double dR1S = (R1/L1)*exp(Et1*wL1)*dt;
                  double dR2S = (R2/L2)*exp(Et2*wL2)*dt;

                  double dRc1J2 = beta*c1*Rc1*pow((xi*I2+J2)/pop,alpha)*dt;
                  double dRc2J1 = beta*c2*Rc2*pow((xi*I1+J1)/pop,alpha)*dt;

                  double dR1J2 = beta*R1*pow((xi*I2+J2)/pop,alpha)*dt;
                  double dR2J1 = beta*R2*pow((xi*I1+J1)/pop,alpha)*dt;
                  
                  double dJ1R = J1/(omega*D)*dt;
                  double dJ2R = J2/(omega*D)*dt;
                  
                  double dRR1 = (R/L2)*exp(Et2*wL2)*dt;
                  double dRR2 = (R/L1)*exp(Et1*wL1)*dt;

                  double dSD = mu*S*dt;
                  double dI1D = mu*I1*dt;
                  double dI2D = mu*I2*dt;
                  double dRc1D = mu*Rc1*dt;
                  double dRc2D = mu*Rc2*dt;
                  double dR1D = mu*R1*dt;
                  double dR2D = mu*R2*dt;
                  double dJ1D = mu*J1*dt;
                  double dJ2D = mu*J2*dt;
                  double dRD = mu*R*dt;
                  
                  S += dBS-dSD-dSI1-dSI2+dR1S+dR2S-2*dASI;
                  I1 += dSI1-dI1Rc1-dI1D+dASI;
                  I2 += dSI2-dI2Rc2-dI2D+dASI;
                  Rc1 += dI1Rc1-dRc1J2-dRc1R1-dRc1D;
                  Rc2 += dI2Rc2-dRc2J1-dRc2R2-dRc2D;
                  R1 += dRc1R1-dR1J2-dR1S+dRR1-dR1D;
                  R2 += dRc2R2-dR2J1-dR2S+dRR2-dR2D;
                  J1 += dRc2J1+dR2J1-dJ1R-dJ1D;
                  J2 += dRc1J2+dR1J2-dJ2R-dJ2D;
                  R += dJ1R+dJ2R-dRR1-dRR2-dRD;
                  C1 += dSI1+dRc2J1+dR2J1+dASI;
                  C2 += dSI2+dRc1J2+dR1J2+dASI;
                  
                  if(S<0){err-=S;S=0;}
                  if(I1<0){err-=I1;I1=0;}
                  if(I2<0){err-=I2;I2=0;}
                  if(Rc1<0){err-=Rc1;Rc1=0;}
                  if(Rc2<0){err-=Rc2;Rc2=0;}
                  if(R1<0){err-=R1;R1=0;}
                  if(R2<0){err-=R2;R2=0;}
                  if(J1<0){err-=J1;J1=0;}
                  if(J2<0){err-=J2;J2=0;}
                  if(R<0){err-=R;R=0;}
                  if(C1<0){err-=C1;C1=0;}
                  if(C2<0){err-=C2;C2=0;}
                  
                  ')

rinit <- Csnippet('
                  double m = pop/(S_0+I1_0+I2_0+Rc1_0+Rc2_0+R1_0+R2_0+J1_0+J2_0+R_0);
                  S = nearbyint(S_0*m);
                  I1 = nearbyint(I1_0*m);
                  I2 = nearbyint(I2_0*m);
                  Rc1 = nearbyint(Rc1_0*m);
                  Rc2 = nearbyint(Rc2_0*m);
                  R1 = nearbyint(R1_0*m);
                  R2 = nearbyint(R2_0*m); 
                  J1 = nearbyint(J1_0*m);
                  J2 = nearbyint(J2_0*m);
                  R = nearbyint(R_0*m);
                  C1 = 0;
                  C2 = 0;
                  err=0;
                  ')

rmeas <- Csnippet("
                  double m1 = phi*C1;
                  double m2 = phi*C2;
                  double v1 = m1*rho;
                  double v2 = m2*rho;
                  double vS1 = m1*rhoS;
                  double vS2 = m2*rhoS;
                  double tol = 1.0e-18;
                  
                  if(((t-floor(t)) > 0.25) && ((t-floor(t)) < 0.75)){
                  H3N = rnorm(m1,vS1+tol);
                  H1N = rnorm(m2,vS2+tol);
                  }else{
                  H3N = rnorm(m1,v1+tol);
                  H1N = rnorm(m2,v2+tol);
                  }
                  
                  if(H3N > 0.0){
                  H3N = nearbyint(H3N);
                  }else{
                  H3N = 0.0;
                  }
                  
                  if(H1N > 0.0){
                  H1N = nearbyint(H1N);
                  }else{
                  H1N = 0.0;
                  }
                  ")

dmeas <- Csnippet("
                  double m1 = phi*C1;
                  double m2 = phi*C2;
                  double v1 = m1*rho;
                  double v2 = m2*rho;
                  double vS1 = m1*rhoS;
                  double vS2 = m2*rhoS;
                  double tol = 1.0e-18;
                  
                  if(m1<=0 || m2<=0){
                  lik = -10000;
                  }else if(((t-floor(t)) > 0.25) && ((t-floor(t)) < 0.75)){
                  lik = dnorm(H3N,m1,vS1+tol,1)+dnorm(H1N,m2,vS2+tol,1);
                  }else{
                  lik = dnorm(H3N,m1,v1+tol,1)+dnorm(H1N,m2,v2+tol,1);
                  }
                  if(!give_log){
                  lik=exp(lik);
                  }
                  
                  ")

##crude estimate
paratransform <- parameter_trans(log=c("sigma","b1","b2","b3","b4","b5","b6",
                                       'xi',"omega",'P1','P2',"c1","c2",
                                       "theta","D","L1","L2","wL1","wL2"),
                                 barycentric=c("S_0","I1_0","I2_0","Rc1_0","Rc2_0",
                                               "R1_0","R2_0","J1_0","J2_0","R_0"),
                                 logit=c("alpha","phi","rho","rhoS"))

covartab <- covariate_table(subset(cov,
                                   select=c("month","pop","dpop",
                                            "month1","month2","month3","month4","month5","month6",
                                            "y11","y12","y13","y21","y22","y23")),
                            times = 'month')

t0 <- inc$month[1]-1/12

inc %>% 
  select(month,H3N,H1N) %>%
  rbind(data.frame(month=t0,H3N=NA,H1N=NA)) %>%
  arrange(month) %>%
  pomp(
    times ='month',
    t0=t0,
    rinit = rinit,
    rprocess=euler(step.fun = rproc,delta.t=1/365),
    rmeasure=rmeas,
    dmeasure=dmeas,
    partrans = paratransform,
    covar=covartab,
    statenames=c("S","I1","I2","Rc1","Rc2","R1","R2","C1","C2","J1","J2","R","W","err"),
    accumvars =c("C1","C2","W","err"),
    paramnames=c("sigma","b1","b2","b3","b4","b5","b6",
                 "theta",'alpha',"D","L1","L2","wL1","wL2",
                 "c1","c2","omega",'xi','P1','P2',
                 "phi","rho","rhoS",
                 "S_0","I1_0","I2_0","Rc1_0","Rc2_0","R1_0","R2_0","J1_0","J2_0","R_0")
  ) -> tsim

##expand.grid----------------------------------------------------------------
foreach (p=iter(params,"row"), 
         .combine=rbind, .inorder=FALSE) %dopar% { 
           library(pomp) 
           best$c1 <- p$c1
           best$c2 <- p$c2
           tsim %>% pfilter(params=best,Np=2000) -> pf 
           p$loglik <- logLik(pf) 
           p } -> profile

write.csv(profile,'profile c1=1 c2=1.csv')

profile %>% 
  mutate(loglik=ifelse(loglik>max(loglik)-20,loglik,NA)) %>%
  ggplot(aes(x=c1,y=c2,z=loglik,fill=loglik))+
  geom_tile(color=NA)+
  scale_fill_gradient()+
  labs(x='c1',y='c2',title = 'c1=1 c2=1')+
  annotate("text", x = 1.75, y = 0.3, 
           label = "c1 (0.66-1.14)",parse=T,color='grey75')+
  annotate("text", x = 1.75, y = 0.5, 
           label = "c2 (0.70-1.20) ",parse=T,color='grey75')

profile %>%
  select(c1,loglik) -> dat.c1
fit.c1 <- loess(dat.c1$loglik~dat.c1$c1,span=0.3)
ranges.c1 <- seq(min(dat.c1$c1),max(dat.c1$c1),0.01)
loglik.c1 <- predict(fit.c1,ranges.c1)
df.c1 <- data.frame(c1=ranges.c1,loglik=loglik.c1)
tem.c1 <- df.c1[which(df.c1$loglik>=max(df.c1$loglik)-1.92),]
min.c1 <- min(tem.c1$c1)
max.c1 <- max(tem.c1$c1)

profile %>%
  select(c2,loglik) %>%
  group_by(c2) %>%
  filter(rank(-loglik)<=20) %>%
  ungroup() -> dat.c2
profile %>%
  select(c2,loglik) -> dat.c2
fit.c2 <- loess(dat.c2$loglik~dat.c2$c2,span=0.35)
ranges.c2 <- seq(min(dat.c2$c2),max(dat.c2$c2),0.1)
loglik.c2 <- predict(fit.c2,ranges.c2)
df.c2 <- data.frame(c2=ranges.c2,loglik=loglik.c2)
tem.c2 <- df.c2[which(df.c2$loglik>=max(df.c2$loglik)-1.92),]
min.c2 <- min(tem.c2$c2)
max.c2 <- max(tem.c2$c2)

profile %>%
  filter(loglik>max(loglik)-20) %>%
  ggplot(aes(x=c2,y=loglik))+
  geom_point()+
  #scale_x_log10()+
  geom_smooth(method="loess",span=0.35)
##Conclude------------------------------------------------------------------

inc1 <- read.csv('sim c1=0.5 c2=1.csv')
inc1 <- inc1[which(inc1$month>=2012.75),]
inc1 <- inc1[which(inc1$month<2019.75),]

inc2 <- read.csv('sim c1=1 c2=0.5.csv')
inc2 <- inc2[which(inc2$month>=2012.75),]
inc2 <- inc2[which(inc2$month<2019.75),]

inc3 <- read.csv('sim c1=1 c2=1.csv')
inc3 <- inc3[which(inc3$month>=2012.75),]
inc3 <- inc3[which(inc3$month<2019.75),]

profile1 <- read.csv('profile c1=0.5 c2=1.csv')
profile2 <- read.csv('profile c1=1 c2=0.5.csv')
profile3 <- read.csv('profile c1=1 c2=1.csv')

##95% CIs------------------------------------------------------------------------------------------------------
span <- 0.5
profile1 %>%
  select(c1,loglik) %>%
  group_by(c1) %>%
  filter(rank(-loglik)<=2) %>%
  ungroup() -> dat.c1
fit.c1 <- loess(dat.c1$loglik~dat.c1$c1,span=span)
ranges.c1 <- seq(min(dat.c1$c1),max(dat.c1$c1),0.01)
loglik.c1 <- predict(fit.c1,ranges.c1)
df.c1 <- data.frame(c1=ranges.c1,loglik=loglik.c1)
tem.c1 <- df.c1[which(df.c1$loglik>=max(df.c1$loglik)-1.92),]
profile1.min.c1 <- min(tem.c1$c1)
profile1.max.c1 <- max(tem.c1$c1)

profile1 %>%
  select(c2,loglik) %>%
  group_by(c2) %>%
  filter(rank(-loglik)<=2) %>%
  ungroup() -> dat.c2
fit.c2 <- loess(dat.c2$loglik~dat.c2$c2,span=span)
ranges.c2 <- seq(min(dat.c2$c2),max(dat.c2$c2),0.01)
loglik.c2 <- predict(fit.c2,ranges.c2)
df.c2 <- data.frame(c2=ranges.c2,loglik=loglik.c2)
tem.c2 <- df.c2[which(df.c2$loglik>=max(df.c2$loglik)-1.92),]
profile1.min.c2 <- min(tem.c2$c2)
profile1.max.c2 <- max(tem.c2$c2)

profile2 %>%
  select(c1,loglik) %>%
  group_by(c1) %>%
  filter(rank(-loglik)<=2) %>%
  ungroup() -> dat.c1
fit.c1 <- loess(dat.c1$loglik~dat.c1$c1,span=span)
ranges.c1 <- seq(min(dat.c1$c1),max(dat.c1$c1),0.01)
loglik.c1 <- predict(fit.c1,ranges.c1)
df.c1 <- data.frame(c1=ranges.c1,loglik=loglik.c1)
tem.c1 <- df.c1[which(df.c1$loglik>=max(df.c1$loglik)-1.92),]
profile2.min.c1 <- min(tem.c1$c1)
profile2.max.c1 <- max(tem.c1$c1)

profile2 %>%
  select(c2,loglik) %>%
  group_by(c2) %>%
  filter(rank(-loglik)<=2) %>%
  ungroup() -> dat.c2
fit.c2 <- loess(dat.c2$loglik~dat.c2$c2,span=span)
ranges.c2 <- seq(min(dat.c2$c2),max(dat.c2$c2),0.01)
loglik.c2 <- predict(fit.c2,ranges.c2)
df.c2 <- data.frame(c2=ranges.c2,loglik=loglik.c2)
tem.c2 <- df.c2[which(df.c2$loglik>=max(df.c2$loglik)-1.92),]
profile2.min.c2 <- min(tem.c2$c2)
profile2.max.c2 <- max(tem.c2$c2)

profile3 %>%
  select(c1,loglik) %>%
  group_by(c1) %>%
  filter(rank(-loglik)<=2) %>%
  ungroup() -> dat.c1
fit.c1 <- loess(dat.c1$loglik~dat.c1$c1,span=span)
ranges.c1 <- seq(min(dat.c1$c1),max(dat.c1$c1),0.01)
loglik.c1 <- predict(fit.c1,ranges.c1)
df.c1 <- data.frame(c1=ranges.c1,loglik=loglik.c1)
tem.c1 <- df.c1[which(df.c1$loglik>=max(df.c1$loglik)-1.92),]
profile3.min.c1 <- min(tem.c1$c1)
profile3.max.c1 <- max(tem.c1$c1)

profile3 %>%
  select(c2,loglik) %>%
  group_by(c2) %>%
  filter(rank(-loglik)<=2) %>%
  ungroup() -> dat.c2
fit.c2 <- loess(dat.c2$loglik~dat.c2$c2,span=span)
ranges.c2 <- seq(min(dat.c2$c2),max(dat.c2$c2),0.01)
loglik.c2 <- predict(fit.c2,ranges.c2)
df.c2 <- data.frame(c2=ranges.c2,loglik=loglik.c2)
tem.c2 <- df.c2[which(df.c2$loglik>=max(df.c2$loglik)-1.92),]
profile3.min.c2 <- min(tem.c2$c2)
profile3.max.c2 <- max(tem.c2$c2)

profile3 %>%
  group_by(c2) %>%
  filter(rank(-loglik)<=2) %>%
  ungroup() %>%
  ggplot(aes(x=c2,y=loglik))+
  geom_point(color='grey')+
  geom_smooth(method="loess",span=span,color='black',size=0.5)+
  geom_hline(aes(yintercept=max(df.c2$loglik)-1.92), linetype="dashed",colour="blue",size=0.5)+
  geom_vline(aes(xintercept=profile3.min.c2), linetype="dashed",colour="blue",size=0.5)+
  geom_vline(aes(xintercept=profile3.max.c2), linetype="dashed",colour="blue",size=0.5)+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

##plot---------------------------------------------------------------------------------------
f1 <- ggplot(inc1,aes(x=month))+
  geom_line(aes(y=H1N),color='blue')+
  geom_line(aes(y=H3N),color='red')+
  labs(title=expression(paste('Scenario I : ',c[H3]==0.5,',',c[H1]==1)),
       x='Time(monthly)',
       y=expression(paste('incidence '%*%10^7)))+
  scale_x_continuous(breaks=seq(2013,2019,2))+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),title=element_text(size=7.5))


f2 <- ggplot(inc2,aes(x=month))+
  geom_line(aes(y=H1N),color='blue')+
  geom_line(aes(y=H3N),color='red')+
  labs(title=expression(paste('Scenario II : ',c[H3]==1,',',c[H1]==0.5)),
       x='Time(monthly)',
       y=expression(paste('incidence '%*%10^7)))+
  scale_x_continuous(breaks=seq(2013,2019,2))+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),title=element_text(size=7.5))


f3 <- ggplot(inc3,aes(x=month))+
  geom_line(aes(y=H1N),color='blue')+
  geom_line(aes(y=H3N),color='red')+
  labs(title=expression(paste('Scenario III : ',c[H3]==1,',',c[H1]==1)),
       x='Time(monthly)',
       y=expression(paste('incidence '%*%10^7)))+
  scale_x_continuous(breaks=seq(2013,2019,2))+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),title=element_text(size=7.5))


cut.loglik <- 50
size <- 2
f4 <- profile1 %>% 
  mutate(loglik=ifelse(loglik>max(loglik)-cut.loglik,loglik,NA)) %>%
  ggplot(aes(x=c1,y=c2,z=loglik,fill=loglik))+
  geom_tile(color=NA)+
  scale_fill_gradient(low = 'LightBlue', high = 'DarkBlue')+
  labs(x=expression(c[H3]),y=expression(c[H1]))+
  annotate("text", x = 1.25, y = 1.85, 
           label = expression(paste(c[H3],' : ','0.35~0.53')),
           parse=T,color='black',size = size)+
  annotate("text", x = 1.25, y = 1.65, 
           label = expression(paste(c[H1],' : ','0.71~1.30')),
           parse=T,color='black',size = size)+
  theme(legend.position = "bottom",legend.key.height = unit(5, "pt"),
        legend.key.width = unit(20, "pt"), legend.text=element_text(size = 5))+
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),title=element_text(size=7.5))


f5 <- profile2 %>% 
  mutate(loglik=ifelse(loglik>max(loglik)-cut.loglik,loglik,NA)) %>%
  ggplot(aes(x=c1,y=c2,z=loglik,fill=loglik))+
  geom_tile(color=NA)+
  scale_fill_gradient(low = 'LightBlue', high = 'DarkBlue')+
  labs(x=expression(c[H3]),y=expression(c[H1]))+
  annotate("text", x = 1.25, y = 1.85, 
           label = expression(paste(c[H3],' : ','0.62~1.04')),           
           parse=T,color='black',size = size)+
  annotate("text", x = 1.25, y = 1.65, 
           label = expression(paste(c[H1],' : ','0.23~0.62')),           
           parse=T,color='black',size = size)+
  theme(legend.position = "bottom",legend.key.height = unit(5, "pt"),
        legend.key.width = unit(20, "pt"), legend.text=element_text(size = 5))+
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),title=element_text(size=7.5))


f6 <- profile3 %>% 
  mutate(loglik=ifelse(loglik>max(loglik)-cut.loglik,loglik,NA)) %>%
  ggplot(aes(x=c1,y=c2,z=loglik,fill=loglik))+
  geom_tile(color=NA)+
  scale_fill_gradient(low = 'LightBlue', high = 'DarkBlue')+
  labs(x=expression(c[H3]),y=expression(c[H1]))+
  annotate("text", x = 1.25, y = 1.85, 
           label = expression(paste(c[H3],' : ','0.68~1.22')),           
           parse=T,color='black',size = size)+
  annotate("text", x = 1.25, y = 1.65, 
           label = expression(paste(c[H1],' : ','0.63~1.29')),           
           parse=T,color='black',size = size)+
  theme(legend.position = "bottom",legend.key.height = unit(5, "pt"),
        legend.key.width = unit(20, "pt"), legend.text=element_text(size = 5))+
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5),title=element_text(size=7.5))

tiff('FigS3.tiff', height = 1500, width = 2000, res = 300, units = 'px', pointsize = 12,compression = 'jpeg')
grid.arrange(f1,f2,f3,f4,f5,f6,nrow=2)
dev.off()