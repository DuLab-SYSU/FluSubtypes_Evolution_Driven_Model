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

read.csv('Austria Incidence.csv') %>%
  filter(month>=2014.75) -> inc
read.csv('Austria covar YearlyDistance.csv') %>%
  filter(month>=2014.75) -> cov

denominator <- 1e+05
inc %>%
  select(month,H1N,H3N) %>%
  melt(id.vars='month',variable.name="subtype",value.name="incidence") %>%
  ggplot(aes(x=month,y=incidence/denominator,group=subtype,shape=subtype,color=subtype)) +
  scale_colour_discrete(name  ="Subtype",breaks=c("H1N", "H3N"),labels=c("H1N1", "H3N2")) +
  scale_color_manual(values = c("blue", "red")) +
  geom_line() +
  geom_point() +
  labs(x='Time(monthly)',y=expression(paste('Incidence'%*%10^5)),title='Austria') +
  theme_bw() +
  theme(legend.position = 'top') -> f.inc

cov %>%
  ggplot(aes(x=month))+
  geom_point(aes(y=y11),color='black') +
  geom_point(aes(y=y12),color='grey50') +
  labs(x='Time(monthly)',y='H3N2 dist',title='Austria') +
  theme_bw() -> f.Et.H3N

cov %>%
  ggplot(aes(x=month))+
  geom_point(aes(y=y21),color='black') +
  geom_point(aes(y=y22),color='grey50') +
  labs(x='Time(monthly)',y='H1N1 dist',title='Austria') +
  theme_bw() -> f.Et.H1N

grid.arrange(f.inc, f.Et.H3N, f.Et.H1N, ncol=1)
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
                  
                  if (ISNA(H1N) || ISNA(H3N)) {
                  lik = (give_log) ? 0 : 1;
                  }else{
                  
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

                  }
                  
                  ")

##crude estimate-------------------------------------------------------------------------------------
sobolDesign(
  lower = c(sigma=0.01,b1=1,b2=1,b3=1,b4=1,b5=1,b6=1,
            theta=0.5,alpha=0.9,D=0.01,L1=1,L2=1,wL1=0.9,wL2=0.09,
            c1=0,c2=1,omega=0,xi=0,
            P1=0.1,P2=0.1,
            S_0=0.1,I1_0=0.0000001,I2_0=0.0000001,R1_0=0.1,R2_0=0.1,Rc1_0=0.1,Rc2_0=0.1,
            J1_0=0.0000001,J2_0=0.0000001,R_0=0.1,
            phi=0.1,rho=0.1,rhoS=0.1),
  upper = c(sigma=0.10,b1=10,b2=10,b3=10,b4=10,b5=10,b6=10,
            theta=3.0,alpha=1.0,D=0.10,L1=100,L2=100,wL1=5,wL2=5,
            c1=1,c2=1,omega=1,xi=1,
            P1=10,P2=10,
            S_0=0.9,I1_0=0.0000002,I2_0=0.0000002,R1_0=0.9,R2_0=0.9,Rc1_0=0.5,Rc2_0=0.5,
            J1_0=0.0000002,J2_0=0.0000002,R_0=0.9,
            phi=1.0,rho=0.9,rhoS=0.9),
  nseq = 3000
) -> guess


paratransform <- parameter_trans(log=c("sigma","b1","b2","b3","b4","b5","b6",
                                       'P1','P2','c2',
                                       "theta","D1","D2","L1","L2","wL1","wL2"),
                                 barycentric=c("S_0","I1_0","I2_0","Rc1_0","Rc2_0",
                                               "R1_0","R2_0","J1_0","J2_0","R_0"),
                                 logit=c("alpha","phi","rho","rhoS",'xi',"omega",'c1'))

covartab <- covariate_table(subset(cov,
                                   select=c("month","pop","dpop",
                                            "month1","month2","month3","month4","month5","month6",
                                            "y11","y12","y21","y22")),
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

##First search-------------------------------------------------------------------------
tsim %>% 
  mif2(params=guess[1,],
       Np=2000,
       Nmif=50,
       rw.sd=rw.sd(sigma=0.02,
                   b1=0.02,b2=0.02,b3=0.02,b4=0.02,b5=0.02,b6=0.02,
                   theta=0.02,alpha=0.02,D1=0.02,D2=0.02,L1=0.02,L2=0.02,wL1=0.02,wL2=0.02,
                   c1=0.02,c2=0,omega=0.02,xi=0.02,P1=0.02,P2=0.02,
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

##Second search--------------------------------------------------------------------
foreach (start=iter(starts,"row"),
         .combine=rbind, .packages=c("pomp"), 
         .errorhandling="remove", .inorder=FALSE) %dopar% {
           mf1 %>%  
            mif2(params=start,rw.sd=rw.sd(sigma=0.01,
                                          b1=0.01,b2=0.01,b3=0.01,b4=0.01,b5=0.01,b6=0.01,
                                          theta=0.01,alpha=0.01,D=0.01,L1=0.01,L2=0.01,wL1=0.01,wL2=0.01,
                                          c1=0.01,c2=0,omega=0.01,xi=0.01,P1=0.02,P2=0.02,
                                          phi=0.01,rho=0.01,rhoS=0.01,
                                          S_0=ivp(0.01),I1_0=ivp(0.01),I2_0=ivp(0.01),Rc1_0=ivp(0.01),Rc2_0=ivp(0.01),
                                          R1_0=ivp(0.01),R2_0=ivp(0.01),J1_0=ivp(0.01),J2_0=ivp(0.01),R_0=ivp(0.01))) %>% 
             mif2() -> mf
            replicate(5, mf %>% pfilter() %>% logLik()) %>%
             logmeanexp(se=TRUE) -> ll
           data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
         } -> estimates2

estimates <- rbind(estimates1,estimates2)
estimates <- arrange(estimates,-loglik)

write.csv(estimates,'AustriaEstimates.csv', row.names = F)

##simulation----------------------------------------------------------------------
estimates %>%
  filter(!is.na(loglik)) %>%
  filter(loglik == max(loglik)) %>%
  select(-loglik,-loglik.se) -> best
tsim %>% simulate(params=best,nsim = 5000) -> sim1

sim1 %>% as.data.frame() %>%
  select(month,H1N) %>%
  group_by(month) %>%
  summarise(mean=quantile(H1N,0.5)/denominator,
            lower=quantile(H1N,0.025)/denominator,
            upper=quantile(H1N,0.975)/denominator) %>%
  as.data.frame() %>%
  filter(month>=2014.75) -> sim1_H1N_summary
sim1_H1N_summary$H1N <- inc$H1N/denominator

sim1 %>% as.data.frame() %>%
  select(month,H3N) %>%
  group_by(month) %>%
  summarise(mean=quantile(H3N,0.5)/denominator,
            lower=quantile(H3N,0.025)/denominator,
            upper=quantile(H3N,0.975)/denominator) %>%
  as.data.frame() %>%
  filter(month>=2014.75) -> sim1_H3N_summary
sim1_H3N_summary$H3N <- inc$H3N/denominator


write.csv(sim1_H1N_summary,'Austria sim1_H1N_summary.csv', row.names = F)
write.csv(sim1_H3N_summary,'Austria sim1_H3N_summary.csv', row.names = F)

figure_H1N <- ggplot(sim1_H1N_summary,aes(x = month,y = mean,ymin = lower,ymax = upper))+
  geom_ribbon(alpha = 0.2,fill = "#0000FF") + 
  geom_line(aes(y = H1N),size = 0.8,color = "black",alpha =0.8) +
  geom_line(aes(y = mean),color = "blue",size = 0.8,alpha =0.9) +
  labs(y=expression(paste('Incidence'%*%10^6)),x="Time(monthly)")+
  annotate("text", x = 2019, y = 4, label = 'H1N1', color='blue')+
  scale_y_continuous(limits = c(0,4.5))+
  theme_bw()  

figure_H3N <- ggplot(sim1_H3N_summary,aes(x = month,y = mean,ymin = lower,ymax = upper))+
  geom_ribbon(alpha = 0.2,fill = "#EE6363") + 
  geom_line(aes(y = H3N),size = 0.8,color = "black",alpha =0.8) +
  geom_line(aes(y = mean),color = "red",size = 0.8,alpha =0.9) +
  labs(y='',x="Time(monthly)")+
  annotate("text", x = 2019, y = 4, label = 'H3N2', color='red')+
  scale_y_continuous(limits = c(0,4.5))+
  theme_bw() 


grid.arrange(figure_H1N,figure_H3N,nrow=1)

##Et------------------------------------------------------------------------------------------
Et <- data.frame(month=inc$month,
                 Et_H3N=cov$y11*exp(-1/best$theta)+cov$y12*exp(-2/best$theta),
                 Et_H1N=cov$y21*exp(-1/best$theta)+cov$y22*exp(-2/best$theta))
write.csv(Et,'AustriaEt.csv', row.names = F)
Et <- read.csv('AustriaEt.csv')
figure_Et <- ggplot(Et, aes(x=month)) +
  geom_line(aes(y=Et_H1N), color='blue') +
  geom_line(aes(y=Et_H3N), color='red') +
  labs(x='Time(monthly)', y='Evolutionary index') +
  scale_x_continuous(breaks = seq(2013,2019,1)) +
  theme_bw()

g1 <- grid.arrange(figure_H1N,figure_H3N,figure_Et,nrow=1)
