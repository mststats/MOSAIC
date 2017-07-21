path="RESULTS/"
lambda<-alpha<-rho<-theta<-Mu<-chrno.G<-list()
i=1
for (chrno in c.l:c.u)
  {
  l.l=1
  tmp=dir(path,glob2rx(paste0(target,"_", firstind, "_", paste(chrnos[c(1,nchrno)],collapse="-"),"_",NN,"_*EMlog.out")))
  chrno.G[[chrno]]=as.integer(strsplit(tmp,"_")[[1]][4])
  tmp<-read.table(paste0(path,tmp),header=T)
  l.u=l.l+kLL*L-1
  tmp.mu<-as.double(tail(tmp,1)[l.l:l.u])
  Mu[[chrno]]<-t(matrix(tmp.mu,L))
  l.l=l.u+1;l.u=l.u+L
  tmp.rho<-as.double(tail(tmp,1)[l.l:l.u])
  rho[[chrno]]<-tmp.rho
  l.l=l.u+1;l.u=l.u+L
  tmp.alpha<-as.double(tail(tmp,1)[l.l:l.u])
  alpha[[chrno]]<-tmp.alpha
  l.l=l.u+1;l.u=l.l+1
  lambda[[chrno]]<-as.double(tail(tmp,1)[l.l:l.u])
  l.l=l.u+1;l.u=l.u+L
  tmp.theta<-as.double(tail(tmp,1)[l.l:l.u])
  theta[[chrno]]<-tmp.theta
  i=i+1
  }

