# for use if hapmix EM already run on this problem
tmp<-read.csv("hapmix_params.txt")
for (ind in 1:NUMI) # all the same
{
  alpha[[ind]][1]<-tmp$theta[1]
  alpha[[ind]][2]<-1-alpha[[ind]][1]
  lambda[[ind]]<-tmp$lambda[1]
}
Q<-create_Q(alpha,lambda,L,dr,NUMI)
rho[1:2]<-1-exp(-dr*c(mean(tmp$rho1), mean(tmp$rho2)))
theta[1]=mean(tmp$mutation1);theta[2]=mean(tmp$mutation2);theta<-theta/(1+theta) # not a 1->1 compatible model w/ hapmix. 
theta<-rep(0.002169197,2)
#phi.theta=mean(tmp$mutation3)
#theta[1,2]=theta[2,1]=phi.theta/(phi.theta+NUMP)
Mu[ANC!=1,1]<-mean(tmp$miscopy1);Mu[1,1]<-1-mean(tmp$miscopy1)
Mu[ANC!=2,2]<-mean(tmp$miscopy2);Mu[2,2]<-1-mean(tmp$miscopy2)
Mu<-t(t(Mu)/colSums(Mu)) # rescale
# update transitions and emissions
for (ind in 1:NUMI)
  transitions[[ind]]<-s_trans(L,kLL,Q[[ind]],Mu,rho,NL)
mutmat<-fmutmat(theta, L, maxmiss, maxmatch)
source("initProb.R")
