output_data_hapmix<-function(ch)
{ 
  options(scipen=10) # make sure no scientific notation is used
  write.table(snps, col.names=F, row.names=F, sep=" ", file=paste0("snpfile.",chrnos[ch]),quote=F)
  write.table(t(hY[which(!is.na(match(label,which(!is.na(match(panels,pops[[1]])))))),]),sep="",file=paste0("panel1genofile.",chrnos[ch]), col.names=F, row.names=F)
  write.table(t(hY[which(!is.na(match(label,which(!is.na(match(panels,pops[[2]])))))),]),sep="",file=paste0("panel2genofile.",chrnos[ch]),col.names=F, row.names=F)
  write.table(t(matrix(hY[!KNOWN,],NUMA)),sep="",file=paste0("admixedhaplofile.",chrnos[ch]),row.names=F,col.names=F) # this is haploid
  if (NUMA>1)
  {
    tmp<-hY[!KNOWN,];tmp<-matrix(tmp[seq(1,nrow(tmp),2),]+tmp[seq(2,nrow(tmp),2),],ncol=S[ch])
    write.table(t(tmp),sep="",file=paste0("admixedgenofile.",chrnos[ch]),row.names=F,col.names=F) # this is proper diploid 
  }
  if (target=="simulated") admixed_truth<-array(true_anc[[ch]],c(L,NUMA,S[ch]))
  if (target!="simulated") admixed_truth<-array(NaN,c(L,NUMA,S[ch]))
  save(admixed_truth, file=paste0("truth.",chrnos[ch],".rdata"))
  write(paste0(":sites:",S[ch]),paste0("rates.",chrnos[ch]))
  write(cbind(snps[,4],rates*100),ncol=S[ch],file=paste0("rates.",chrnos[ch]),append=T) 
  command<-paste0("mv ", "snpfile.",chrnos[ch]," ", "panel*",chrnos[ch]," ", "truth.",chrnos[ch],".rdata ", 
		  "admixed*",chrnos[ch]," ", "rates.",chrnos[ch]," ", "../../HapmixReleasev2/RAINBOW/")
  system(command)
  options(scipen=0) # restore defaults
  }
output_hapmix_initial<-function(ch)
{
  options(scipen=10) # make sure no scientific notation is used
  theta1<-o.theta[1];theta1<-theta1/(1-theta1) 
  theta2<-o.theta[2];theta2<-theta2/(1-theta2) 
  command<-paste("sh create_parfile.sh", (Reduce("+",o.alpha)/NUMI)[1], Reduce("+",o.lambda)/NUMI, # works; don't change!
		 -log(1-o.rho[1])/dr/100*sum(NL[which(!is.na(match(panels,pops[[1]])))]), -log(1-o.rho[2])/dr/100*sum(NL[which(!is.na(match(panels,pops[[2]])))]), # works; don't change!
		 sum(NL[which(!is.na(match(panels,pops[[1]])))])*theta1,sum(NL[which(!is.na(match(panels,pops[[2]])))])*theta2,# works; don't change! 
		 mean(o.theta[1]/(1-o.theta[1]),o.theta[2]/(1-o.theta[2])), 
		 sum(o.Mu[which(!is.na(match(panels,pops[[2]]))),1]), sum(o.Mu[which(!is.na(match(panels,pops[[1]]))),2]), chrnos[ch]) # works; don't change!
  system(command)
  command<-paste0("mv ", "mosaic*",chrnos[ch]," ", "../../HapmixReleasev2/RAINBOW/")
  system(command)
  options(scipen=0) # restore defaults
}

output_hapmix_final<-function(ch)
{
  options(scipen=10) # make sure no scientific notation is used
  # for use after EM
  theta1<-mean(theta[1]*Mu[,1]*NL[1:kLL]);theta1<-theta1/(1-theta1)
  theta2<-mean(theta[2]*Mu[,2]*NL[1:kLL]);theta2<-theta2/(1-theta2)
  theta3<-mean(theta%*%t(Mu)*NL[1:kLL])/NUMP;theta3<-theta3/(1-theta3);
  command<-paste("sh create_parfile.sh", (Reduce("+",alpha)/NUMI)[1], Reduce("+",lambda)/NUMI, 
		 -log(1-rho[1])/dr/100*sum(NL[which(!is.na(match(panels,pops[[1]])))]), -log(1-rho[2])/dr/100*sum(NL[which(!is.na(match(panels,pops[[2]])))]),
		 theta1,theta2,theta3,
		 sum(Mu[which(!is.na(match(panels,pops[[2]]))),1]), sum(Mu[which(!is.na(match(panels,pops[[1]]))),2]), chrnos[ch])
  system(command)
  command<-paste0("mv ", "mosaic*",chrnos[ch]," ", "../../HapmixReleasev2/RAINBOW/")
  system(command)
  options(scipen=0) # restore defaults
}

output_hapmix_final<-function(ch)
{
  options(scipen=10) # make sure no scientific notation is used
  theta1<-theta[1];theta1<-theta1/(1-theta1) 
  theta2<-theta[2];theta2<-theta2/(1-theta2) 
  command<-paste("sh create_parfile.sh", (Reduce("+",alpha)/NUMI)[1], Reduce("+",lambda)/NUMI, # works; don't change!
		 -log(1-rho[1])/dr/100*sum(NL[which(!is.na(match(panels,pops[[1]])))]), -log(1-rho[2])/dr/100*sum(NL[which(!is.na(match(panels,pops[[2]])))]), # works; don't change!
		 sum(NL[which(!is.na(match(panels,pops[[1]])))])*theta1,sum(NL[which(!is.na(match(panels,pops[[2]])))])*theta2,# works; don't change! 
		 mean(theta[1]/(1-theta[1]),theta[2]/(1-theta[2])), 
		 sum(Mu[which(!is.na(match(panels,pops[[2]]))),1]), sum(Mu[which(!is.na(match(panels,pops[[1]]))),2]), chrnos[ch]) # works; don't change!
  system(command)
  command<-paste0("mv ", "mosaic*",chrnos[ch]," ", "../../HapmixReleasev2/RAINBOW/")
  system(command)
  options(scipen=0) # restore defaults
}
