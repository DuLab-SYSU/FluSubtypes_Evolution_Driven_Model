remove(list=ls())

library(tidyverse)
library(pomp)
library(dplyr)
library(foreach)
library(doParallel)
library(iterators)
library(doRNG)
library(gridExtra)
library(cowplot)
library(ggplot2)
library(plot3D)
library(scatterplot3d)
library(rgl)

registerDoParallel()

inc <- read.csv('USA-inc-month-2019.csv')
inc <- inc[which(inc$month>=2012.75),]
inc <- inc[which(inc$month<2019.75),]

cov <- read.csv('USA-covar-2019.csv')
cov <- cov[which(cov$month>=2012.75),]
cov <- cov[which(cov$month<2019.75),]

denominator <- 1e+07
##plot inc--------------------------------------------------------------------------------------------
figure_inc <- ggplot(inc,aes(x = month))+
  geom_line(aes(y = H1N/denominator),size = 0.5,color = "blue",alpha =0.8) +
  geom_point(aes(y = H1N/denominator),size = 0.5,color = "blue",alpha =0.8) +
  geom_line(aes(y = H3N/denominator),size = 0.5,color = "red",alpha =0.8) +
  geom_point(aes(y = H3N/denominator),size = 0.5,color = "red",alpha =0.8) +
  geom_line(aes(y = B/denominator),size = 0.5,color = "green",alpha =0.8) +
  geom_point(aes(y = B/denominator),size = 0.5,color = "green",alpha =0.8) +
  labs(y=expression(paste('Incidence'%*%10^7)),x="Time(monthly)")+
  scale_x_continuous(breaks = seq(2013,2019,2)) +
  scale_y_continuous(limits = c(0,2))+
  annotate("text", x = 2019, y = 1.75, label = 'H1N1', color='blue', size=2)+
  annotate("text", x = 2019, y = 1.65, label = 'H3N2', color='red', size=2)+
  annotate("text", x = 2019, y = 1.55, label = 'B', color='green', size=2)+
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))

##constrcut the model-------------------------------------------------------------------------------------------------------------
rproc <- Csnippet('
                  double mu = 0.015;
                  double dW = rgammawn(sigma,dt);
                  double A = 0.1*365;
                  double beta = exp(b1*month1+b2*month2+b3*month3+b4*month4+b5*month5+b6*month6)*(dW/dt);
                  double Et1 = y11*exp(-1/theta)+y12*exp(-2/theta);//where subscrip 1 denotes H3N2
                  double Et2 = y21*exp(-1/theta)+y22*exp(-2/theta);//where subscrip 2 denotes H3N2
                  
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

##initial priors--------------------------------------------------------------------------------------------------
sobolDesign(
  lower = c(sigma=0.01,b1=1,b2=1,b3=1,b4=1,b5=1,b6=1,
            theta=0.5,alpha=0.9,L1=10,L2=10,wL1=0.9,wL2=0.09,
            c1=0,c2=0,omega=0,xi=0,P1=0,P2=0,D=0,
            S_0=0.1,I1_0=0.0000001,I2_0=0.0000001,R1_0=0.1,R2_0=0.1,Rc1_0=0.1,Rc2_0=0.1,
            J1_0=0.0000001,J2_0=0.0000001,R_0=0.1,
            phi=0.1,rho=0.1,rhoS=0.1),
  upper = c(sigma=0.10,b1=10,b2=10,b3=10,b4=10,b5=10,b6=10,
            theta=3.0,alpha=1.0,L1=100,L2=100,wL1=5,wL2=5,
            c1=3,c2=3,omega=3,xi=5,P1=1,P2=1,D=0.05,
            S_0=0.9,I1_0=0.0000002,I2_0=0.0000002,R1_0=0.9,R2_0=0.9,Rc1_0=0.5,Rc2_0=0.5,
            J1_0=0.0000002,J2_0=0.0000002,R_0=0.9,
            phi=1.0,rho=0.9,rhoS=0.9),
  nseq = 5000
) -> guess

paratransform <- parameter_trans(log=c("sigma","b1","b2","b3","b4","b5","b6",
                                       'xi',"omega",'P1','P2',"c1","c2",'D',
                                       "theta","L1","L2","wL1","wL2"),
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
                 "theta",'alpha',"L1","L2","wL1","wL2",
                 "c1","c2","omega",'xi','P1','P2','D',
                 "phi","rho","rhoS",
                 "S_0","I1_0","I2_0","Rc1_0","Rc2_0","R1_0","R2_0","J1_0","J2_0","R_0")
  ) -> tsim

##First search---------------------------------------------------------------------------------------------------------
tsim %>% 
  mif2(params=guess[1,],
       Np=2000,
       Nmif=50,
       rw.sd=rw.sd(sigma=0.02,
                   b1=0.02,b2=0.02,b3=0.02,b4=0.02,b5=0.02,b6=0.02,
                   theta=0.02,alpha=0.02,L1=0.02,L2=0.02,wL1=0.02,wL2=0.02,
                   c1=0.02,c2=0.02,omega=0.02,xi=0.02,P1=0.02,P2=0.02,D=0.02,
                   phi=0.02,rho=0.02,rhoS=0.02,
                   S_0=ivp(0.02),I1_0=ivp(0.02),I2_0=ivp(0.02),Rc1_0=ivp(0.02),Rc2_0=ivp(0.02),
                   R1_0=ivp(0.02),R2_0=ivp(0.02),J1_0=ivp(0.02),J2_0=ivp(0.02),R_0=ivp(0.02)),
       cooling.fraction.50=0.5) -> mf1

foreach (p=iter(guess,"row"),
         .combine=c, .packages=c("pomp"),
         .errorhandling="remove", .inorder=FALSE) %dopar% {
           
           mf1 %>% mif2(params=p)
           
         } -> mfs 

foreach (mf=mfs,
         .combine=rbind, .packages=c("pomp"), 
         .errorhandling="remove", .inorder=FALSE) %dopar% {
           replicate(5, mf %>% pfilter() %>% logLik()) %>%
             logmeanexp(se=TRUE) -> ll
           data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
         } -> estimates1

estimates1 %>%
  filter(!is.na(loglik)) %>%
  filter(loglik > max(loglik)-50) %>%
  select(-loglik,-loglik.se) -> starts

##Second search----------------------------------------------------------------------------------------------------
foreach (start=iter(starts,"row"),
         .combine=rbind, .packages=c("pomp"), 
         .errorhandling="remove", .inorder=FALSE) %dopar% {
           mf1 %>%  
             mif2(params=start,rw.sd=rw.sd(sigma=0.01,
                                           b1=0.01,b2=0.01,b3=0.01,b4=0.01,b5=0.01,b6=0.01,
                                           theta=0.01,alpha=0.01,L1=0.01,L2=0.01,wL1=0.01,wL2=0.01,
                                           c1=0.01,c2=0.01,omega=0.01,xi=0.01,P1=0.01,P2=0.01,D=0.01,
                                           phi=0.01,rho=0.01,rhoS=0.01,
                                           S_0=ivp(0.01),I1_0=ivp(0.01),I2_0=ivp(0.01),Rc1_0=ivp(0.01),Rc2_0=ivp(0.01),
                                           R1_0=ivp(0.01),R2_0=ivp(0.01),J1_0=ivp(0.01),J2_0=ivp(0.01),R_0=ivp(0.01)),
                  cooling.fraction.50=0.3) %>% mif2() -> mf
           replicate(5, mf %>% pfilter() %>% logLik()) %>%
             logmeanexp(se=TRUE) -> ll
           data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
         } -> estimates2

estimates <- rbind(estimates1,estimates2)
estimates <- arrange(estimates,-loglik)

write.csv(estimates, 'Parameters_estimation.csv', row.names = F)
##simulation--------------------------------------------------------------------------------------------------------
estimates %>%
  filter(!is.na(loglik)) %>%
  filter(loglik == max(loglik)) %>%
  select(-loglik,-loglik.se) -> best  

Et <- data.frame(month=inc$month,
                 Et_H3N=cov$y11*exp(-1/best$theta)+cov$y12*exp(-2/best$theta),
                 Et_H1N=cov$y21*exp(-1/best$theta)+cov$y22*exp(-2/best$theta))

Figure.Et <- ggplot(Et, aes(x=month)) +
  geom_line(aes(y=Et_H1N), color='blue') +
  geom_line(aes(y=Et_H3N), color='red') +
  labs(x='Time(monthly)', y='Evolutionary index') +
  scale_x_continuous(breaks = seq(2013,2019,2)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))

tsim %>% simulate(params=best,nsim = 5000) -> sim1

sim1 %>% as.data.frame() %>%
  select(month,H1N) %>%
  group_by(month) %>%
  summarise(mean=quantile(H1N,0.5)/denominator,
            lower=quantile(H1N,0.025)/denominator,
            upper=quantile(H1N,0.975)/denominator) %>%
  as.data.frame() %>%
  filter(month>=2012.75) -> sim1_H1N_summary
sim1_H1N_summary$H1N <- inc$H1N/denominator

sim1 %>% as.data.frame() %>%
  select(month,H3N) %>%
  group_by(month) %>%
  summarise(mean=quantile(H3N,0.5)/denominator,
            lower=quantile(H3N,0.025)/denominator,
            upper=quantile(H3N,0.975)/denominator) %>%
  as.data.frame() %>%
  filter(month>=2012.75) -> sim1_H3N_summary
sim1_H3N_summary$H3N <- inc$H3N/denominator

figure_H1N <- ggplot(sim1_H1N_summary,aes(x = month,y = mean,ymin = lower,ymax = upper))+
  geom_ribbon(alpha = 0.2,fill = "#0000FF") + 
  geom_line(aes(y = H1N),size = 0.8,color = "black",alpha =0.8) +
  geom_line(aes(y = mean),color = "blue",size = 0.8,alpha =0.9) +
  labs(y=expression(paste('Incidence'%*%10^7)),x="Time(monthly)")+
  scale_x_continuous(breaks = seq(2013,2019,2)) +
  scale_y_continuous(limits = c(0,6.5))+
  annotate("text", x = 2019, y = 5.5, label = 'H1N1', color='blue', size=2)+
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))

figure_H3N <- ggplot(sim1_H3N_summary,aes(x = month,y = mean,ymin = lower,ymax = upper))+
  geom_ribbon(alpha = 0.2,fill = "#EE6363") + 
  geom_line(aes(y = H3N),size = 0.8,color = "black",alpha =0.8) +
  geom_line(aes(y = mean),color = "red",size = 0.8,alpha =0.9) +
  labs(x="Time(monthly)", y='')+
  scale_x_continuous(breaks = seq(2013,2019,2)) +
  scale_y_continuous(limits = c(0,6.5))+
  annotate("text", x = 2019, y = 5.5, label = 'H3N2', color='red',size=2)+
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5))


g.sim <- grid.arrange(figure_H1N,figure_H3N,nrow=1)

tiff('Fig2.tiff', height = 1500, width = 2000, res = 350, units = 'px', pointsize = 12, compression = 'jpeg')
ggdraw() +
  draw_plot(figure_inc, 0, 1/2, width = 1, height = 1/2) +
  draw_plot(g.sim, 0, 0, width = 2/3, height = 1/2) +
  draw_plot(Figure.Et, 2/3, 0, width = 1/3, height = 1/2) +
  draw_plot_label(label = c('A','B','C'), x=c(0,0,2/3), y=c(1,1/2,1/2), size = 10)
dev.off()

estimates <- read.csv('Parameters_estimation.csv')
estimates %>%
  filter(loglik>=max(loglik)-5) %>%
  filter(loglik.se<=2) %>%
  select(-loglik,-loglik.se) %>%
  apply(2,range) -> ranges

##profiling design-----------------------------------------------------------------------------------
profileDesign(
  c1=seq(from=0,to=3,length=60),
  lower=ranges[1,-15],
  upper=ranges[2,-15],
  nprof=100,type="runif"
) -> begins.c1

profileDesign(
  c2=seq(from=0,to=3,length=60),
  lower=ranges[1,-16],
  upper=ranges[2,-16],
  nprof=100,type="runif"
) -> begins.c2

profileDesign(
  omega=seq(from=0,to=2,length=80),
  lower=ranges[1,-17],
  upper=ranges[2,-17],
  nprof=50,type="runif"
) -> begins.omega


profileDesign(
  xi=seq(from=0,to=2,length=40),
  lower=ranges[1,-18],
  upper=ranges[2,-18],
  nprof=50,type="runif"
) -> begins.xi

profileDesign(
  L1=seq(from=0,to=240,length=40),
  lower=ranges[1,-11],
  upper=ranges[2,-11],
  nprof=100,type="runif"
) -> begins.L1


profileDesign(
  L2=seq(from=0,to=240,length=96),
  lower=ranges[1,-12],
  upper=ranges[2,-12],
  nprof=50,type="runif"
) -> begins.L2


profileDesign(
  P1=seq(from=ranges[1,19],to=ranges[2,19],length=80),
  lower=ranges[1,-19],
  upper=ranges[2,-19],
  nprof=50,type="runif"
) -> begins.P1


profileDesign(
  P2=seq(from=ranges[1,20],to=ranges[2,20],length=80),
  lower=ranges[1,-20],
  upper=ranges[2,-20],
  nprof=50,type="runif"
) -> begins.P2

profileDesign(
  D=seq(from=0,to=0.2,length=40),
  lower=ranges[1,-10],
  upper=ranges[2,-10],
  nprof=100,type="runif"
) -> begins.D


profileDesign(
  theta=seq(from=0,to=5,length=100),
  lower=ranges[1,-8],
  upper=ranges[2,-8],
  nprof=50,type="runif"
) -> begins.theta


profileDesign(
  wL1=seq(from=0,to=20,length=80),
  lower=ranges[1,-13],
  upper=ranges[2,-13],
  nprof=50,type="runif"
) -> begins.wL1


profileDesign(
  wL2=seq(from=0,to=15,length=60),
  lower=ranges[1,-14],
  upper=ranges[2,-14],
  nprof=50,type="runif"
) -> begins.wL2


profileDesign(
  alpha=seq(from=0.5,to=1,length=20),
  lower=ranges[1,-9],
  upper=ranges[2,-9],
  nprof=100,type="runif"
) -> begins.alpha

##profiling search-----------------------------------------------------------------------------------------------
foreach (begin=iter(begins.c1,"row"),
         .combine=rbind, 
         .errorhandling="remove", .inorder=FALSE) %dopar% {
           
           mf1 %>%
             mif2(
               params=begin,
               rw.sd=rw.sd(sigma=0.02,
                           b1=0.02,b2=0.02,b3=0.02,b4=0.02,b5=0.02,b6=0.02,
                           theta=0.02,alpha=0.02,D=0.02,L1=0.02,L2=0.02,wL1=0.02,wL2=0.02,
                           c2=0.02,omega=0.02,xi=0.02,P1=0.02,P2=0.02,
                           phi=0.02,rho=0.02,rhoS=0.02,
                           S_0=ivp(0.02),I1_0=ivp(0.02),I2_0=ivp(0.02),Rc1_0=ivp(0.02),Rc2_0=ivp(0.02),
                           R1_0=ivp(0.02),R2_0=ivp(0.02),J1_0=ivp(0.02),J2_0=ivp(0.02),R_0=ivp(0.02))
               
             ) %>% mif2(cooling.fraction.50=0.3) -> mf
           
           replicate(5, mf %>% pfilter() %>% logLik()) %>% logmeanexp(se=TRUE) -> ll
           data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
         } -> c1.prof


foreach (begin=iter(begins.c2,"row"),
         .combine=rbind, 
         .errorhandling="remove", .inorder=FALSE) %dopar% {
           library(pomp)
           library(tidyverse)
           
           mf1 %>%
             mif2(
               params=begin,
               rw.sd=rw.sd(sigma=0.02,
                           b1=0.02,b2=0.02,b3=0.02,b4=0.02,b5=0.02,b6=0.02,
                           theta=0.02,alpha=0.02,D=0.02,L1=0.02,L2=0.02,wL1=0.02,wL2=0.02,
                           c1=0.02,omega=0.02,xi=0.02,P1=0.02,P2=0.02,
                           phi=0.02,rho=0.02,rhoS=0.02,
                           S_0=ivp(0.02),I1_0=ivp(0.02),I2_0=ivp(0.02),Rc1_0=ivp(0.02),Rc2_0=ivp(0.02),
                           R1_0=ivp(0.02),R2_0=ivp(0.02),J1_0=ivp(0.02),J2_0=ivp(0.02),R_0=ivp(0.02))
               
             ) -> mf
           
           replicate(5, mf %>% pfilter() %>% logLik()) %>% logmeanexp(se=TRUE) -> ll
           data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
         } -> c2.prof


foreach (begin=iter(begins.omega,"row"),
         .combine=rbind, 
         .errorhandling="remove", .inorder=FALSE) %dopar% {
           
           mf1 %>%
             mif2(
               params=begin,
               rw.sd=rw.sd(sigma=0.02,
                           b1=0.02,b2=0.02,b3=0.02,b4=0.02,b5=0.02,b6=0.02,
                           theta=0.02,alpha=0.02,D=0.02,L1=0.02,L2=0.02,wL1=0.02,wL2=0.02,
                           c1=0.02,c2=0.02,xi=0.02,P1=0.02,P2=0.02,
                           phi=0.02,rho=0.02,rhoS=0.02,
                           S_0=ivp(0.02),I1_0=ivp(0.02),I2_0=ivp(0.02),Rc1_0=ivp(0.02),Rc2_0=ivp(0.02),
                           R1_0=ivp(0.02),R2_0=ivp(0.02),J1_0=ivp(0.02),J2_0=ivp(0.02),R_0=ivp(0.02))
             ) %>% mif2(cooling.fraction.50=0.3) -> mf
           
           replicate(5, mf %>% pfilter() %>% logLik()) %>% logmeanexp(se=TRUE) -> ll
           data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
         } -> omega.prof


foreach (begin=iter(begins.xi,"row"),
         .combine=rbind, 
         .errorhandling="remove", .inorder=FALSE) %dopar% {
           
           mf1 %>%
             mif2(
               params=begin,
               rw.sd=rw.sd(sigma=0.02,
                           b1=0.02,b2=0.02,b3=0.02,b4=0.02,b5=0.02,b6=0.02,
                           theta=0.02,alpha=0.02,D=0.02,L1=0.02,L2=0.02,wL1=0.02,wL2=0.02,
                           c1=0.02,omega=0.02,c2=0.02,P1=0.02,P2=0.02,
                           phi=0.02,rho=0.02,rhoS=0.02,
                           S_0=ivp(0.02),I1_0=ivp(0.02),I2_0=ivp(0.02),Rc1_0=ivp(0.02),Rc2_0=ivp(0.02),
                           R1_0=ivp(0.02),R2_0=ivp(0.02),J1_0=ivp(0.02),J2_0=ivp(0.02),R_0=ivp(0.02))
               
             ) %>% mif2(cooling.fraction.50=0.3) -> mf
           
           replicate(5, mf %>% pfilter() %>% logLik()) %>% logmeanexp(se=TRUE) -> ll
           data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
         } -> xi.prof


foreach (begin=iter(begins.L1,"row"),
         .combine=rbind, 
         .errorhandling="remove", .inorder=FALSE) %dopar% {
           
           mf1 %>%
             mif2(
               params=begin,
               rw.sd=rw.sd(sigma=0.02,
                           b1=0.02,b2=0.02,b3=0.02,b4=0.02,b5=0.02,b6=0.02,
                           theta=0.02,alpha=0.02,D=0.02,L2=0.02,wL1=0.02,wL2=0.02,
                           c1=0.02,c2=0.02,omega=0.02,xi=0.02,P1=0.02,P2=0.02,
                           phi=0.02,rho=0.02,rhoS=0.02,
                           S_0=ivp(0.02),I1_0=ivp(0.02),I2_0=ivp(0.02),Rc1_0=ivp(0.02),Rc2_0=ivp(0.02),
                           R1_0=ivp(0.02),R2_0=ivp(0.02),J1_0=ivp(0.02),J2_0=ivp(0.02),R_0=ivp(0.02))
               
             ) %>% mif2() -> mf
           
           replicate(5, mf %>% pfilter() %>% logLik()) %>% logmeanexp(se=TRUE) -> ll
           data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
         } -> L1.prof

foreach (begin=iter(begins.L2,"row"),
         .combine=rbind, 
         .errorhandling="remove", .inorder=FALSE) %dopar% {
           
           mf1 %>%
             mif2(
               params=begin,
               rw.sd=rw.sd(sigma=0.02,
                           b1=0.02,b2=0.02,b3=0.02,b4=0.02,b5=0.02,b6=0.02,
                           theta=0.02,alpha=0.02,D=0.02,L1=0.02,wL1=0.02,wL2=0.02,
                           c1=0.02,c2=0.02,omega=0.02,xi=0.02,P1=0.02,P2=0.02,
                           phi=0.02,rho=0.02,rhoS=0.02,
                           S_0=ivp(0.02),I1_0=ivp(0.02),I2_0=ivp(0.02),Rc1_0=ivp(0.02),Rc2_0=ivp(0.02),
                           R1_0=ivp(0.02),R2_0=ivp(0.02),J1_0=ivp(0.02),J2_0=ivp(0.02),R_0=ivp(0.02))
               
             ) %>% mif2() -> mf
           
           replicate(5, mf %>% pfilter() %>% logLik()) %>% logmeanexp(se=TRUE) -> ll
           data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
         } -> L2.prof


foreach (begin=iter(begins.P1,"row"),
         .combine=rbind, 
         .errorhandling="remove", .inorder=FALSE) %dopar% {
           
           mf1 %>%
             mif2(
               params=begin,
               rw.sd=rw.sd(sigma=0.02,
                           b1=0.02,b2=0.02,b3=0.02,b4=0.02,b5=0.02,b6=0.02,
                           theta=0.02,alpha=0.02,D=0.02,L1=0.02,L2=0.02,wL1=0.02,wL2=0.02,
                           c1=0.02,c2=0.02,omega=0.02,xi=0.02,P2=0.02,
                           phi=0.02,rho=0.02,rhoS=0.02,
                           S_0=ivp(0.02),I1_0=ivp(0.02),I2_0=ivp(0.02),Rc1_0=ivp(0.02),Rc2_0=ivp(0.02),
                           R1_0=ivp(0.02),R2_0=ivp(0.02),J1_0=ivp(0.02),J2_0=ivp(0.02),R_0=ivp(0.02))
               
             ) %>% mif2() -> mf
           
           replicate(5, mf %>% pfilter() %>% logLik()) %>% logmeanexp(se=TRUE) -> ll
           data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
         } -> P1.prof


foreach (begin=iter(begins.P2,"row"),
         .combine=rbind, 
         .errorhandling="remove", .inorder=FALSE) %dopar% {
           library(pomp)
           library(tidyverse)
           
           mf1 %>%
             mif2(
               params=begin,
               rw.sd=rw.sd(sigma=0.02,
                           b1=0.02,b2=0.02,b3=0.02,b4=0.02,b5=0.02,b6=0.02,
                           theta=0.02,alpha=0.02,D=0.02,L1=0.02,L2=0.02,wL1=0.02,wL2=0.02,
                           c1=0.02,c2=0.02,omega=0.02,xi=0.02,P1=0.02,
                           phi=0.02,rho=0.02,rhoS=0.02,
                           S_0=ivp(0.02),I1_0=ivp(0.02),I2_0=ivp(0.02),Rc1_0=ivp(0.02),Rc2_0=ivp(0.02),
                           R1_0=ivp(0.02),R2_0=ivp(0.02),J1_0=ivp(0.02),J2_0=ivp(0.02),R_0=ivp(0.02))
               
             ) %>% mif2() -> mf
           
           replicate(5, mf %>% pfilter() %>% logLik()) %>% logmeanexp(se=TRUE) -> ll
           data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
         } -> P2.prof

foreach (begin=iter(begins.D,"row"),
         .combine=rbind, 
         .errorhandling="remove", .inorder=FALSE) %dopar% {
           
           mf1 %>%
             mif2(
               params=begin,
               rw.sd=rw.sd(sigma=0.02,
                           b1=0.02,b2=0.02,b3=0.02,b4=0.02,b5=0.02,b6=0.02,
                           theta=0.02,alpha=0.02,L1=0.02,L2=0.02,wL1=0.02,wL2=0.02,
                           c1=0.02,c2=0.02,omega=0.02,xi=0.02,P1=0.02,P2=0.02,
                           phi=0.02,rho=0.02,rhoS=0.02,
                           S_0=ivp(0.02),I1_0=ivp(0.02),I2_0=ivp(0.02),Rc1_0=ivp(0.02),Rc2_0=ivp(0.02),
                           R1_0=ivp(0.02),R2_0=ivp(0.02),J1_0=ivp(0.02),J2_0=ivp(0.02),R_0=ivp(0.02))
               
             ) %>% mif2() -> mf
           
           replicate(5, mf %>% pfilter() %>% logLik()) %>% logmeanexp(se=TRUE) -> ll
           data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
         } -> D.prof


foreach (begin=iter(begins.theta,"row"),
         .combine=rbind, 
         .errorhandling="remove", .inorder=FALSE) %dopar% {
           
           mf1 %>%
             mif2(
               params=begin,
               rw.sd=rw.sd(sigma=0.02,
                           b1=0.02,b2=0.02,b3=0.02,b4=0.02,b5=0.02,b6=0.02,
                           alpha=0.02,D=0.02,L1=0.02,L2=0.02,wL1=0.02,wL2=0.02,
                           c1=0.02,c2=0.02,omega=0.02,xi=0.02,P1=0.02,P2=0.02,
                           phi=0.02,rho=0.02,rhoS=0.02,
                           S_0=ivp(0.02),I1_0=ivp(0.02),I2_0=ivp(0.02),Rc1_0=ivp(0.02),Rc2_0=ivp(0.02),
                           R1_0=ivp(0.02),R2_0=ivp(0.02),J1_0=ivp(0.02),J2_0=ivp(0.02),R_0=ivp(0.02))
               
             ) %>% mif2(cooling.fraction.50=0.3) -> mf
           
           replicate(5, mf %>% pfilter() %>% logLik()) %>% logmeanexp(se=TRUE) -> ll
           data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
         } -> theta.prof


foreach (begin=iter(begins.wL1,"row"),
         .combine=rbind, 
         .errorhandling="remove", .inorder=FALSE) %dopar% {
           
           mf1 %>%
             mif2(
               params=begin,
               rw.sd=rw.sd(sigma=0.02,
                           b1=0.02,b2=0.02,b3=0.02,b4=0.02,b5=0.02,b6=0.02,
                           theta=0.02,alpha=0.02,D=0.02,L1=0.02,L2=0.02,wL2=0.02,
                           c1=0.02,c2=0.02,omega=0.02,xi=0.02,P1=0.02,P2=0.02,
                           phi=0.02,rho=0.02,rhoS=0.02,
                           S_0=ivp(0.02),I1_0=ivp(0.02),I2_0=ivp(0.02),Rc1_0=ivp(0.02),Rc2_0=ivp(0.02),
                           R1_0=ivp(0.02),R2_0=ivp(0.02),J1_0=ivp(0.02),J2_0=ivp(0.02),R_0=ivp(0.02))
               
             ) %>% mif2() -> mf
           
           replicate(5, mf %>% pfilter() %>% logLik()) %>% logmeanexp(se=TRUE) -> ll
           data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
         } -> wL1.prof

foreach (begin=iter(begins.wL2,"row"),
         .combine=rbind, 
         .errorhandling="remove", .inorder=FALSE) %dopar% {
           
           mf1 %>%
             mif2(
               params=begin,
               rw.sd=rw.sd(sigma=0.02,
                           b1=0.02,b2=0.02,b3=0.02,b4=0.02,b5=0.02,b6=0.02,
                           theta=0.02,alpha=0.02,D=0.02,L1=0.02,L2=0.02,wL1=0.02,
                           c1=0.02,c2=0.02,omega=0.02,xi=0.02,P1=0.02,P2=0.02,
                           phi=0.02,rho=0.02,rhoS=0.02,
                           S_0=ivp(0.02),I1_0=ivp(0.02),I2_0=ivp(0.02),Rc1_0=ivp(0.02),Rc2_0=ivp(0.02),
                           R1_0=ivp(0.02),R2_0=ivp(0.02),J1_0=ivp(0.02),J2_0=ivp(0.02),R_0=ivp(0.02))
               
             ) %>% mif2(cooling.fraction.50=0.3) -> mf
           
           replicate(5, mf %>% pfilter() %>% logLik()) %>% logmeanexp(se=TRUE) -> ll
           data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
         } -> wL2.prof



foreach (begin=iter(begins.alpha,"row"),
         .combine=rbind, 
         .errorhandling="remove", .inorder=FALSE) %dopar% {
           
           mf1 %>%
             mif2(
               params=begin,
               rw.sd=rw.sd(sigma=0.02,
                           b1=0.02,b2=0.02,b3=0.02,b4=0.02,b5=0.02,b6=0.02,
                           theta=0.02,D=0.02,L1=0.02,L2=0.02,wL1=0.02,wL2=0.02,
                           c1=0.02,c2=0.02,omega=0.02,xi=0.02,P1=0.02,P2=0.02,
                           phi=0.02,rho=0.02,rhoS=0.02,
                           S_0=ivp(0.02),I1_0=ivp(0.02),I2_0=ivp(0.02),Rc1_0=ivp(0.02),Rc2_0=ivp(0.02),
                           R1_0=ivp(0.02),R2_0=ivp(0.02),J1_0=ivp(0.02),J2_0=ivp(0.02),R_0=ivp(0.02))
               
             ) %>% mif2(cooling.fraction.50=0.3) -> mf
           
           replicate(5, mf %>% pfilter() %>% logLik()) %>% logmeanexp(se=TRUE) -> ll
           data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
         } -> alpha.prof

rank <- 2
span <- 0.5;## span <- 0.75

##profiling plot-------------------------------------------------------------------------------------------
alpha.prof <- read.csv('profile alpha.csv')
theta.prof <- read.csv('profile theta.csv')
P1.prof <- read.csv('profile P1.csv')
P2.prof <- read.csv('profile P2.csv')
L1.prof <- read.csv('profile L1.csv')
L2.prof <- read.csv('profile L2.csv')
wL1.prof <- read.csv('profile wL1.csv')
wL2.prof <- read.csv('profile wL2.csv')
c1.prof <- read.csv('profile c1.csv')
c2.prof <- read.csv('profile c2.csv')
xi.prof <- read.csv('profile xi.csv')
omega.prof <- read.csv('profile omega.csv')


alpha.prof %>%
  filter(loglik>=max(loglik)-50) %>%
  group_by(alpha) %>%
  filter(rank(-loglik)<=rank) %>%
  ungroup() %>%
  filter(loglik.se<=2) -> dat.alpha
fit.alpha <- loess(dat.alpha$loglik~dat.alpha$alpha,span=span)
ranges.alpha <- seq(min(dat.alpha$alpha),max(dat.alpha$alpha),0.001)
loglik.alpha <- predict(fit.alpha,ranges.alpha)
df.alpha <- data.frame(alpha=ranges.alpha,loglik=loglik.alpha)
tem.alpha <- df.alpha[which(df.alpha$loglik>=max(df.alpha$loglik)-1.92),]
min.alpha <- min(tem.alpha$alpha)
max.alpha <- max(tem.alpha$alpha)

theta.prof %>%
  filter(loglik>=max(loglik)-50) %>%
  group_by(theta) %>%
  filter(rank(-loglik)<=rank) %>%
  ungroup() %>%
  filter(loglik.se<=2) -> dat.theta
fit.theta <- loess(dat.theta$loglik~dat.theta$theta,span=span)
ranges.theta <- seq(min(dat.theta$theta),max(dat.theta$theta),0.05)
loglik.theta <- predict(fit.theta,ranges.theta)
df.theta <- data.frame(theta=ranges.theta,loglik=loglik.theta)
tem.theta <- df.theta[which(df.theta$loglik>=max(df.theta$loglik)-1.92),]
min.theta <- min(tem.theta$theta)
max.theta <- max(tem.theta$theta)


P1.prof %>%
  filter(loglik>=max(loglik)-50) %>%
  filter(P1<=7) %>%
  group_by(P1) %>%
  filter(rank(-loglik)<=rank) %>%
  ungroup() %>%
  filter(loglik.se<=2) -> dat.P1
fit.P1 <- loess(dat.P1$loglik~dat.P1$P1,span=span)
ranges.P1 <- seq(min(dat.P1$P1),max(dat.P1$P1),0.1)
loglik.P1 <- predict(fit.P1,ranges.P1)
df.P1 <- data.frame(P1=ranges.P1,loglik=loglik.P1)
tem.P1 <- df.P1[which(df.P1$loglik>=max(df.P1$loglik)-1.92),]
min.P1 <- min(tem.P1$P1)
max.P1 <- max(tem.P1$P1)

P2.prof %>%
  filter(loglik>=max(loglik)-50) %>%
  group_by(P2) %>%
  filter(rank(-loglik)<=rank) %>%
  ungroup() %>%
  filter(loglik.se<=2) -> dat.P2
fit.P2 <- loess(dat.P2$loglik~dat.P2$P2,span=span)
ranges.P2 <- seq(min(dat.P2$P2),max(dat.P2$P2),0.1)
loglik.P2 <- predict(fit.P2,ranges.P2)
df.P2 <- data.frame(P2=ranges.P2,loglik=loglik.P2)
tem.P2 <- df.P2[which(df.P2$loglik>=max(df.P2$loglik)-1.92),]
min.P2 <- min(tem.P2$P2)
max.P2 <- max(tem.P2$P2)


L1.prof %>%
  filter(loglik>=max(loglik)-20) %>%
  filter(L1<=230) %>%
  group_by(L1) %>%
  filter(rank(-loglik)<=rank) %>%
  ungroup() %>%
  filter(loglik.se<=2) -> dat.L1
fit.L1 <- loess(dat.L1$loglik~dat.L1$L1,span=span)
ranges.L1 <- seq(min(dat.L1$L1),max(dat.L1$L1),0.1)
loglik.L1 <- predict(fit.L1,ranges.L1)
df.L1 <- data.frame(L1=ranges.L1,loglik=loglik.L1)
tem.L1 <- df.L1[which(df.L1$loglik>=max(df.L1$loglik)-1.92),]
min.L1 <- min(tem.L1$L1)
max.L1 <- max(tem.L1$L1)


L2.prof %>%
  filter(loglik>=max(loglik)-40) %>%
  group_by(L2) %>%
  filter(rank(-loglik)<=rank) %>%
  ungroup() %>%
  filter(loglik.se<=2) -> dat.L2
fit.L2 <- loess(dat.L2$loglik~dat.L2$L2,span=span)
ranges.L2 <- seq(min(dat.L2$L2),max(dat.L2$L2),0.1)
loglik.L2 <- predict(fit.L2,ranges.L2)
df.L2 <- data.frame(L2=ranges.L2,loglik=loglik.L2)
tem.L2 <- df.L2[which(df.L2$loglik>=max(df.L2$loglik)-1.92),]
min.L2 <- min(tem.L2$L2)
max.L2 <- max(tem.L2$L2)

wL1.prof %>%
  filter(loglik>=max(loglik)-50) %>%
  group_by(wL1) %>%
  filter(rank(-loglik)<=rank) %>%
  ungroup() %>%
  filter(loglik.se<=2) -> dat.wL1
fit.wL1 <- loess(dat.wL1$loglik~dat.wL1$wL1,span=span)
ranges.wL1 <- seq(min(dat.wL1$wL1),max(dat.wL1$wL1),0.05)
loglik.wL1 <- predict(fit.wL1,ranges.wL1)
df.wL1 <- data.frame(wL1=ranges.wL1,loglik=loglik.wL1)
tem.wL1 <- df.wL1[which(df.wL1$loglik>=max(df.wL1$loglik)-1.92),]
min.wL1 <- min(tem.wL1$wL1)
max.wL1 <- max(tem.wL1$wL1)


wL2.prof %>%
  filter(loglik>=max(loglik)-50) %>%
  group_by(wL2) %>%
  filter(rank(-loglik)<=rank) %>%
  ungroup() %>%
  filter(loglik.se<=2) -> dat.wL2
fit.wL2 <- loess(dat.wL2$loglik~dat.wL2$wL2,span=span)
ranges.wL2 <- seq(min(dat.wL2$wL2),max(dat.wL2$wL2),0.05)
loglik.wL2 <- predict(fit.wL2,ranges.wL2)
df.wL2 <- data.frame(wL2=ranges.wL2,loglik=loglik.wL2)
tem.wL2 <- df.wL2[which(df.wL2$loglik>=max(df.wL2$loglik)-1.92),]
min.wL2 <- min(tem.wL2$wL2)
max.wL2 <- max(tem.wL2$wL2)


c1.prof %>%
  filter(c1>0) %>%
  filter(c1<=3) %>%
  group_by(c1) %>%
  filter(rank(-loglik)<=rank) %>%
  ungroup() %>%
  filter(loglik.se<=2) -> dat.c1
fit.c1 <- loess(dat.c1$loglik~dat.c1$c1,span=span)
ranges.c1 <- seq(min(dat.c1$c1),max(dat.c1$c1),0.01)
loglik.c1 <- predict(fit.c1,ranges.c1)
df.c1 <- data.frame(c1=ranges.c1,loglik=loglik.c1)
tem.c1 <- df.c1[which(df.c1$loglik>=max(df.c1$loglik)-1.92),]
min.c1 <- min(tem.c1$c1)
max.c1 <- max(tem.c1$c1)
c1.label <- sprintf('(%.2f~%.2f)',round(min.c1,2),round(max.c1,2))

f.c1 <- c1.prof %>%
  filter(c1>0) %>%
  filter(c1<=3) %>%
  group_by(c1) %>%
  filter(rank(-loglik)<=rank) %>%
  ungroup() %>%
  filter(loglik.se<=2) %>%
  ggplot(aes(x=c1,y=loglik,
             ymin=loglik-2*loglik.se,ymax=loglik+2*loglik.se))+
  geom_point(color='LightGrey',size=0.5)+
  geom_errorbar(color='LightGrey')+
  geom_smooth(method="loess",span=span,color='black',size=0.5)+
  geom_hline(aes(yintercept=max(df.c1$loglik)-1.92), linetype="dashed",colour="blue",size=0.5)+
  geom_vline(aes(xintercept=min.c1), linetype="dashed",colour="blue",size=0.5)+
  geom_vline(aes(xintercept=max.c1), linetype="dashed",colour="blue",size=0.5)+
  labs(title = c1.label, size=0.5) +
  xlab(expression(c[H3])) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5)) +
  theme(plot.title = element_text(size=7.5))

theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(f.c1)


c2.prof %>%
  filter(loglik>=max(loglik)-50) %>%
  filter(c2<=3) %>%
  group_by(c2) %>%
  filter(rank(-loglik)<=rank) %>%
  ungroup() %>%
  filter(loglik.se<=2) -> dat.c2
fit.c2 <- loess(dat.c2$loglik~dat.c2$c2,span=span)
ranges.c2 <- seq(min(dat.c2$c2),max(dat.c2$c2),0.01)
loglik.c2 <- predict(fit.c2,ranges.c2)
df.c2 <- data.frame(c2=ranges.c2,loglik=loglik.c2)
tem.c2 <- df.c2[which(df.c2$loglik>=max(df.c2$loglik)-1.92),]
min.c2 <- min(tem.c2$c2)
max.c2 <- max(tem.c2$c2)
c2.label <- sprintf('(%.2f~%.2f)',round(min.c2,2),round(max.c2,2))

f.c2 <- c2.prof %>%
  filter(loglik>=max(loglik)-50) %>%
  filter(c2<=3) %>%
  group_by(c2) %>%
  filter(rank(-loglik)<=rank) %>%
  ungroup() %>%
  filter(loglik.se<=2) %>%
  ggplot(aes(x=c2,y=loglik,
             ymin=loglik-2*loglik.se,ymax=loglik+2*loglik.se))+
  geom_point(color='LightGrey', size=0.5)+
  geom_errorbar(color='LightGrey', size=0.5)+
  geom_smooth(method="loess",span=span,color='black',se=T,size=0.5)+
  geom_hline(aes(yintercept=max(df.c2$loglik)-1.92), linetype="dashed",colour="blue",size=0.5)+
  geom_vline(aes(xintercept=min.c2), linetype="dashed",colour="blue",size=0.5)+
  geom_vline(aes(xintercept=max.c2), linetype="dashed",colour="blue",size=0.5)+
  labs(title = c2.label) +
  xlab(expression(c[H1])) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 7.5), axis.title.y = element_text(size = 7.5)) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5)) +
  theme(plot.title = element_text(size=7.5))

theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
print(f.c2)

xi.prof %>%
  filter(loglik>=max(loglik)-60) %>%
  group_by(xi) %>%
  filter(rank(-loglik)<=rank) %>%
  ungroup() -> dat.xi
fit.xi <- loess(dat.xi$loglik~dat.xi$xi,span=span)
ranges.xi <- sexi(min(dat.xi$xi),max(dat.xi$xi),0.01)
loglik.xi <- predict(fit.xi,ranges.xi)
df.xi <- data.frame(xi=ranges.xi,loglik=loglik.xi)
tem.xi <- df.xi[which(df.xi$loglik>=max(df.xi$loglik)-1.92),]
min.xi <- min(tem.xi$xi)
max.xi <- max(tem.xi$xi)

omega.prof %>%
  filter(loglik>=max(loglik)-50) %>%
  group_by(omega) %>%
  filter(rank(-loglik)<=rank) %>%
  ungroup() %>%
  filter(loglik.se<=2) -> dat.omega
fit.omega <- loess(dat.omega$loglik~dat.omega$omega,span=span)
ranges.omega <- seq(min(dat.omega$omega),max(dat.omega$omega),0.01)
loglik.omega <- predict(fit.omega,ranges.omega)
df.omega <- data.frame(omega=ranges.omega,loglik=loglik.omega)
tem.omega <- df.omega[which(df.omega$loglik>=max(df.omega$loglik)-1.92),]
min.omega <- min(tem.omega$omega)
max.omega <- max(tem.omega$omega)

##Fig 2--------------------------------------------------------------------------------------------------------------------
g.c1c2 <- ggarrange(f.c1, f.c2, nrow = 1)

estimates <- read.csv('Parameters_estimation.csv')
tem.ests <- estimates[c(1:10),]
tem.ests %>%
  filter(loglik>=max(loglik)-5) %>%
  filter(loglik.se<2) %>%
  select(-loglik,-loglik.se) %>%
  apply(2,quantile)  -> PredefinedEstimates
best <- PredefinedEstimates[3,]

c1 <- seq(0,1,0.05)
c2 <- seq(0,1,0.05)
res1 <- data.frame(c1=c(NA),c2=c(NA),sum.H1N=c(NA),sum.H3N=c(NA))
write.csv(res1,'c1 c2 sum inc.csv')

best <- PredefinedEstimates[3,]
row <- 1
for (i in 1:length(c1)) {
  for (j in 1:length(c2)) {
    best$c1 <- c1[i]
    best$c2 <- c2[j]
    
    tsim %>% simulate(params=best,nsim = 2000) -> sim1
    
    sim1 %>% as.data.frame() %>%
      select(month,H1N,H3N) %>%
      group_by(month) %>%
      summarise(H1N=quantile(H1N,0.5),
                H3N=quantile(H3N,0.5)) %>%
      as.data.frame() %>%
      filter(month>=2012.75) -> sim1.inc
    
    res1[row,1] <- c1[i]
    res1[row,2] <- c2[j]
    res1[row,3] <- sum(sim1.inc$H1N/(1e+07))
    res1[row,4] <- sum(sim1.inc$H3N/(1e+07))
    row <- row + 1
    print(row)
  }
}

res1 <- read.csv('c1 c2 sum inc.csv')
res1$sum.H1N.rate <- NA
res1$sum.H3N.rate <- NA

for (i in 2:nrow(res1)) {
  res1$sum.H1N.rate[i] <- res1$sum.H1N[i]/res1$sum.H1N[1]
  res1$sum.H3N.rate[i] <- res1$sum.H3N[i]/res1$sum.H3N[1]
}


f.H1N <- res1 %>% 
  ggplot(aes(x=c1,y=c2,z=sum.H1N.rate,fill=sum.H1N.rate))+
  geom_tile(colour = "NA")+
  scale_fill_gradient()+
  labs(x=expression(c[H3]),
       y=expression(c[H1]),
       title = 'H1N1',
       fill='Fold change of incidence'
       #fill=expression(paste('H3N2 incidence'%*%10^7))
  )+
  theme_bw()+
  theme(axis.title.x =element_text(size=7.5), 
        axis.title.y=element_text(size=7.5),
        axis.text.x =element_text(size=7.5), 
        axis.text.y=element_text(size=7.5),
        title=element_text(size=7.5),
        legend.title=element_text(size=7.5),
        legend.text = element_text(size=7.5),
        legend.position = "bottom",
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(10, "pt"))+
  scale_fill_gradient2(low = "blue",mid = "grey50",midpoint = 4,high = "red")

f.H3N <- res1 %>% 
  ggplot(aes(x=c1,y=c2,z=sum.H3N.rate,fill=sum.H3N.rate))+
  geom_tile(colour = "NA")+
  scale_fill_gradient()+
  labs(x=expression(c[H3]),
       y=expression(c[H1]),
       title = 'H3N2',
       fill='Fold change of incidence'
       #fill=expression(paste('H3N2 incidence'%*%10^7))
  )+
  theme_bw()+
  theme(axis.title.x =element_text(size=7.5), 
        axis.title.y=element_text(size=7.5),
        axis.text.x =element_text(size=7.5), 
        axis.text.y=element_text(size=7.5),
        title=element_text(size=7.5),
        legend.title=element_text(size=7.5),
        legend.text = element_text(size=7.5),
        legend.position = "bottom",
        legend.key.height = unit(5, "pt"),
        legend.key.width = unit(10, "pt"))+
  scale_fill_gradient2(low = "blue",mid = "grey50",midpoint = 4,high = "red")


grid_arrange_shared_legend <- function(..., ncol = 2, nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

g.sum.inc.c1c2 <- grid_arrange_shared_legend(f.H1N,f.H3N)

tiff('Fig2.tiff', height = 1500, width = 2000, res = 350, compression = 'jpeg')
ggdraw() +
  draw_plot(g.c1c2, 0, 1/2, width = 1, height = 1/2) +
  draw_plot(g.sum.inc.c1c2, 0, 0, width = 1, height = 1/2) +
  draw_plot_label(label = c('A','B'), x=c(0,0), y=c(1,1/2), size = 10)
dev.off()

