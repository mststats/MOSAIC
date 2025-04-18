# functions to calculate coancestry curves along target genomes based on relative probabilities of pairs of ancestries at various genomic distances
r_create_coancs<-function(t.localanc, gap, MODE="DIP", min.cM=0, max.cM=50,gby=5)
{
  # coarser grid here will mean more averaging / smoothing of curves
  nchrno=length(t.localanc)
  maxconsidered=as.integer(max.cM/gap/100)
  for (ch in 1:nchrno) {
    if (maxconsidered > dim(t.localanc[[ch]])[3]) {
      max.cM=100*gap*dim(t.localanc[[ch]])[3]
      warning("########## trying to fit curves on longer segments than total length of chromosome : changing max.cM to ", max.cM, " cM #########", immediate.=TRUE)
    }
  }
  drange=seq(as.integer(min.cM/gap/100),as.integer(max.cM/gap/100),by=gby) # 50 cM max; drange is in gridpoints
  NUMA=dim(t.localanc[[1]])[2]
  NUMI=max(NUMA/2,1)
  lpop<-dim(t.localanc[[1]])[1]
  if (MODE=="HAP") 
  {
    relprobs<-array(NaN,c(lpop,lpop,NUMA,length(drange))) #store in an array that's lpop x lpop x NUMA x length(drange) 
    ancprobs<-array(NaN,c(lpop,NUMA)) #store in an array that's lpop x NUMA 
  }
  if (MODE=="DIP") 
  {
    relprobs<-array(NaN,c(lpop,lpop,NUMI,length(drange))) #store in an array that's lpop x lpop x NUMI x length(drange) 
    ancprobs<-array(NaN,c(lpop,NUMI)) #store in an array that's lpop x NUMI 
  }

  if (MODE=="HAP"){
    ans=foreach(k=1:NUMA) %dopar%
    {
      kancprobs=rep(0,lpop)
      krelprobs<-array(NaN,c(lpop,lpop,length(drange))) #store in an array that's lpop x lpop x NUMA x length(drange) 
      anclist=list()
      for(ch in 1:nchrno){
	anclist[[ch]]=matrix(nrow=dim(t.localanc[[ch]])[3],ncol=lpop)
	for(i in 1:lpop) anclist[[ch]][,i]=t.localanc[[ch]][i,k,]
      }
      for (i in 1:lpop){
	t1=0
	t2=0
	for(ch in 1:nchrno){
	  t1=t1+sum(anclist[[ch]][,i])
	  t2=t2+length(anclist[[ch]][,i])
	}
	kancprobs[i]=t1/t2
      }
      for (d in 1:length(drange))
      {
	avec=bvec=list()
	for (ch in 1:nchrno) # if we try to go beyond the length of a chromosome this will fail...
	{
	  avec[[ch]]=1:(nrow(anclist[[ch]])-round(drange[d])) # start at left, go to drange[d] from right
	  bvec[[ch]]=(round(drange[d])+1):nrow(anclist[[ch]]) # start at right, go to drange[d] from left 
	}
	for (i in 1:(lpop-1))
	  for (j in (i+1):lpop){
	    #cat(k,i,j,d,"\n")
	    totproda=totlefta=totrighta=0
	    totprodb=totleftb=totrightb=0
	    totlength=0    	
	    for(ch in 1:nchrno){
	      tmpanc=anclist[[ch]]
	      leftvec=tmpanc[,i][avec[[ch]]]
	      rightvec=tmpanc[,j][bvec[[ch]]]
	      leftvecb=tmpanc[,j][avec[[ch]]]
	      rightvecb=tmpanc[,i][bvec[[ch]]]
	      totproda=totproda+sum(leftvec*rightvec)
	      totlefta=totlefta+sum(leftvec)
	      totrighta=totrighta+sum(rightvec)
	      totlength=totlength+length(leftvec)
	      totprodb=totprodb+sum(leftvecb*rightvecb)
	      totleftb=totleftb+sum(leftvecb)
	      totrightb=totrightb+sum(rightvecb)
	    }
	    krelprobs[j,i,d]=krelprobs[i,j,d]=(totproda+totprodb)*totlength/(totlefta*totrighta+totleftb*totrightb)
	  }
	for (i in 1:lpop) # i==j is a special case and therefore faster 
	{
	  totproda=totlefta=totrighta=0
	  totlength=0    	
	  for(ch in 1:nchrno){
	    tmpanc=anclist[[ch]][,i]
	    leftvec=tmpanc[avec[[ch]]]
	    rightvecb=tmpanc[bvec[[ch]]]
	    totproda=totproda+sum(leftvec*rightvecb)
	    totlefta=totlefta+sum(leftvec)
	    totrighta=totrighta+sum(rightvecb)
	    totlength=totlength+length(leftvec)
	  }
	  krelprobs[i,i,d]=(totproda)*totlength/(totlefta*totrighta)
	}
      }
      list(krel=krelprobs, kanc=kancprobs)
    }
      for (k in 1:NUMA)
      {
	relprobs[,,k,]=ans[[k]]$krel
	ancprobs[,k]=ans[[k]]$kanc
      }
  }
  if (MODE=="DIP") { # this only works well when using long chromosomes!
    ans=foreach(ind=1:NUMI) %dopar%
    {
      hap<-c(ind*2-1,ind*2)
      kancprobs=rep(0,lpop)
      krelprobs<-array(NaN,c(lpop,lpop,length(drange))) #store in an array that's lpop x lpop x NUMA x length(drange) 
      anclist=list()
      for(ch in 1:nchrno){
	anclist[[ch]]=matrix(nrow=dim(t.localanc[[ch]])[3],ncol=lpop)
	for(i in 1:lpop) anclist[[ch]][,i]=t.localanc[[ch]][i,hap[1],]+t.localanc[[ch]][i,hap[2],]
      }
      for (i in 1:lpop){
	t1=0
	t2=0
	for(ch in 1:nchrno){
	  t1=t1+sum(anclist[[ch]][,i])
	  t2=t2+length(anclist[[ch]][,i])
	}
	kancprobs[i]=t1/t2
      }
      for (d in 1:length(drange))
      {
	avec=bvec=list()
	for (ch in 1:nchrno) {
	  avec[[ch]]=1:(nrow(anclist[[ch]])-round(drange[d])) # start at left, go to drange[d] from right
	  bvec[[ch]]=(round(drange[d])+1):nrow(anclist[[ch]]) # start at right, go to drange[d] from left 
	}
	for (i in 1:(lpop-1))
	  for (j in (i+1):lpop)
	  {
	    totproda=totlefta=totrighta=0
	    totprodb=totleftb=totrightb=0
	    totlength=0    	
	    for(ch in 1:nchrno){
	      tmpanc=anclist[[ch]]
	      leftvec=tmpanc[,i][avec[[ch]]]
	      rightvec=tmpanc[,j][bvec[[ch]]]
	      leftvecb=tmpanc[,j][avec[[ch]]]
	      rightvecb=tmpanc[,i][bvec[[ch]]]
	      totproda=totproda+sum(leftvec*rightvec)
	      totlefta=totlefta+sum(leftvec)
	      totrighta=totrighta+sum(rightvec)
	      totlength=totlength+length(leftvec)
	      totprodb=totprodb+sum(leftvecb*rightvecb)
	      totleftb=totleftb+sum(leftvecb)
	      totrightb=totrightb+sum(rightvecb)
	    }
	    krelprobs[j,i,d]=krelprobs[i,j,d]=(totproda+totprodb)*totlength/(totlefta*totrighta+totleftb*totrightb)
	  }
	for (i in 1:lpop) # i==j is a special case and therefore faster 
	{
	  totproda=totlefta=totrighta=0
	  totlength=0    	
	  for(ch in 1:nchrno){
	    tmpanc=anclist[[ch]][,i]
	    leftvec=tmpanc[avec[[ch]]]
	    rightvecb=tmpanc[bvec[[ch]]]
	    totproda=totproda+sum(leftvec*rightvecb)
	    totlefta=totlefta+sum(leftvec)
	    totrighta=totrighta+sum(rightvecb)
	    totlength=totlength+length(leftvec)
	  }
	  krelprobs[i,i,d]=(totproda)*totlength/(totlefta*totrighta)
	}
      }
      #cat("ind:", ind, "dim=", dim(krelprobs), "\n")
      list(krel=krelprobs, kanc=kancprobs)
    }
      for (ind in 1:NUMI)
      {
	#cat("relprob dims for ind", ind, "=", dim(relprobs[,,ind,]),dim(ans[[ind]]$krel),"\n")
	relprobs[,,ind,]=ans[[ind]]$krel
	ancprobs[,ind]=ans[[ind]]$kanc
      }
  }
  relprobs[is.na(relprobs)]=1 # remove uninformative ones
  list(relprobs=relprobs,ancprobs=ancprobs,drange=drange)
}
create_coancs<-cmpfun(r_create_coancs,list(optimize=3))
#create_coancs<-r_create_coancs

plot_coanccurves<-function(coancs,gap,lwd=2,cexa=2,k=NULL,popnames=NULL,PLOT=TRUE,targetname=NULL,dd=NULL,min.cM=1,max.cM=NULL,ylab="relative prob.",
			   plotall=(is.null(k)),axisall=FALSE,transalpha=0.5,verbose=FALSE,anc.thresh=0.2,asym=FALSE,samedates=FALSE,ndates=1,optmethod="BFGS")
{
  lpop<-dim(coancs$relprobs)[1]
  if (any(rowMeans(coancs$ancprobs)/lpop<0.05))
    warning("########## minor ancestry proportions small; may be hard to fit estimate coancestry curves ##########",immediate.=TRUE)
  # FLAG should use a vector of ndates, one per pair...
  if (ndates>2) {warning("########## cannot fit more than two dates per pair of ancestries ########## ",immediate.=TRUE);ndates=1}
  if (samedates) {if (ndates==2) warning("########## if samedates specified then cannot fit multiple dates per pair of ancestries ########## ", immediate.=TRUE); ndates=1}
  # plotall indicates whether to plot individual based curves as well as consensus curves
  # axisall indicates whether to use a y-axis limit based on consensus or all curves
  # asym=TRUE allows asymptote to be something other than 1
  # samedates=TRUE forces single estimation of all curves; not sensible for (>2)-way
  if (min.cM!=0) 
  {
    min.cM=which(coancs$drange<(min.cM/gap/100))
    coancs$drange=coancs$drange[-min.cM] # always remove leftmost entries of relprobs i.e. at distance~=0
    coancs$relprobs=array(coancs$relprobs[,,,-min.cM],c(dim(coancs$relprobs)[1],dim(coancs$relprobs)[2],dim(coancs$relprobs)[3],
							dim(coancs$relprobs)[4]-length(min.cM)))
  }
  if (!is.null(max.cM))
  {
    max.cM=which(coancs$drange>(max.cM/gap/100))
    if (length(max.cM)>0) {
      coancs$drange=coancs$drange[-max.cM] # always remove leftmost entries of relprobs i.e. at distance~=0
      coancs$relprobs=array(coancs$relprobs[,,,-max.cM],c(dim(coancs$relprobs)[1],dim(coancs$relprobs)[2],dim(coancs$relprobs)[3],
							dim(coancs$relprobs)[4]-length(max.cM)))
    } else (warning("########## cannot select this distance (you would need to recompute the coancestry curves); using max distance of ", 
			  max(coancs$drange)*100*dr, " ##########", immediate.=TRUE))
  }
  if (samedates& lpop>2)
    warning("########## using samedates for more than 2-way event is not sensible ##########", immediate.=TRUE)
  relcurve=array(0,c(lpop,lpop,dim(coancs$relprobs)[4]))
  kweights=keep=list()
  if (is.null(popnames)) popnames<-1:lpop
  params=array(0,c(lpop,lpop,1+2*ndates))
  if (PLOT)
  {
    if (is.null(dd)) 
    {
      dd=c(NaN,NaN)
      if (lpop==2) {dd[1]=1;dd[2]=3} # dimensions of plots
      if (lpop==3) {dd[1]=2;dd[2]=3} # dimensions of plots
      if (lpop==4) {dd[1]=2;dd[2]=5} # dimensions of plots
      if (lpop==5) {dd[1]=3;dd[2]=5} # dimensions of plots
      if (lpop==6) {dd[1]=3;dd[2]=7} # dimensions of plots
    }
    par(mfrow=c(dd[1],dd[2]), cex.axis=0.7*cexa, cex.lab=cexa, cex.main=cexa, mar=c(cexa,cexa,cexa,0)+0.2)
  }
  for (i in 1:lpop) for (j in 1:lpop)
  {
    if(!is.null(k)) tmpcurve=coancs$relprobs[i,j,k,]
    if (is.null(k))
    { 
      ###only individuals with some of ancestry pair i,j
      #print(coancs$ancprobs) # A X (NUMA or NUMI) marginal probs of being each anc for each target
      zzi=coancs$ancprobs[i,]
      zzj=coancs$ancprobs[j,]
      zzi=zzi^2;zzj=zzj^2
      #if (i!=j) {zzi=zzi^2;zzj=zzj^2} 
      #if (i==j) {zzj=1-zzi;zzi=zzi^3} # 1-p
      # weight by each ind by 
      kweights[[(i-1)*lpop+j]]=rep(0,dim(coancs$relprobs)[3])
      for (kk in 1:length(kweights[[(i-1)*lpop+j]])) # use kk b/c otherwise will create a k we don't want
	kweights[[(i-1)*lpop+j]][kk]=(zzi[kk]*zzj[kk])
      keep[[(i-1)*lpop+j]]=(zzi>mean(zzi)*anc.thresh & zzj>mean(zzj)*anc.thresh)
      if (all(!keep[[(i-1)*lpop+j]])) # if none pass the threshold, use all
	keep[[(i-1)*lpop+j]][]=TRUE
      if (verbose)
	cat("for ancs", i,":",j, "keeping", which(keep[[(i-1)*lpop+j]]), "\n")
      kweights[[(i-1)*lpop+j]][!keep[[(i-1)*lpop+j]]]=0
      #kweights[[(i-1)*lpop+j]][]=0;kweights[[(i-1)*lpop+j]][1]=1 # debug by looking at only one ind / hap
      if (all(kweights[[(i-1)*lpop+j]]==0)) kweights[[(i-1)*lpop+j]][]=1
      kweights[[(i-1)*lpop+j]]=kweights[[(i-1)*lpop+j]]/sum(kweights[[(i-1)*lpop+j]])
      #cat(i,j,kweights[[(i-1)*lpop+j]],"\n")
      tmpcurve=rep(0,dim(coancs$relprobs)[4])
      for (kk in 1:length(kweights[[(i-1)*lpop+j]])) # use kk b/c otherwise will create a k we don't want
	coancs$relprobs[i,j,kk,]=coancs$relprobs[i,j,kk,]*kweights[[(i-1)*lpop+j]][kk]
      if (dim(coancs$relprobs)[3]==1) tmpcurve=coancs$relprobs[i,j,1,] # only one ind / hap
      if (dim(coancs$relprobs)[3]>1) tmpcurve=colSums(coancs$relprobs[i,j,,]) # multiple inds / haps
    }
    relcurve[i,j,]=tmpcurve
  }
  x=array(NaN,c(lpop,lpop,1+2*ndates));x[,,1]=1 # it'll always be approx 1 for x[1]
  for (i in 1:lpop) for (j in 1:lpop)
  {
    #if (diff(range(relcurve[i,j,]))<1e-3)
    #{
    #  print("no change along entire genome: exiting")
    #  return(NULL)
    #}
    # fit an exp to the curve
    closest=relcurve[i,j,which.min(abs(relcurve[i,j,]-1))]
    perc=0.5;someway=(1-perc)*relcurve[i,j,1]+perc*closest # someway b/w closest to 1 and first point on relcurve
    dh=ifelse(relcurve[i,j,1]<1,which((relcurve[i,j,]-someway)>0)[1],which((relcurve[i,j,]-someway)<0)) # first point to cross someway 
    #if (verbose) cat("closest:", closest, perc, "way", someway, "dh:", dh, "at:", gap*100*coancs$drange[dh], "curve:", relcurve[i,j,dh], "\n")
    for (dd in 1:ndates) x[i,j,2*dd]=relcurve[i,j,][1]-x[i,j,1] # at d=0, Y=x[1]+x[2] # match at d=0 to get x[2]
    for (dd in 1:ndates) x[i,j,2*dd+1]=5
    if (!is.na(dh))
      for (dd in 1:ndates)
	x[i,j,2*dd+1]=sqrt(abs(log(((relcurve[i,j,dh]-x[i,j,1])/x[i,j,2*dd]))/coancs$drange[dh]/gap))
    for (dd in 1:ndates) if (is.nan(x[i,j,2*dd]) | is.infinite(x[i,j,2*dd])) x[i,j,2*dd]=0 
      if (is.nan(x[i,j,1]) | is.infinite(x[i,j,1])) x[i,j,1]=1 
    for (dd in 1:ndates) {
      if (is.nan(x[i,j,2*dd+1]) | is.infinite(x[i,j,2*dd+1]) | (x[i,j,2*dd+1]==0)) x[i,j,2*dd+1]=5 
    }
    #print(x[i,j,3])
    if (verbose) cat(i,":",j, "before:", x[i,j,],"\n")
  }
  if (!samedates)
  {
    for (i in 1:lpop) for (j in 1:lpop)
    {
      if (!asym)
      {
	if (ndates==1) { # number of dates per paired mixture
	  mf=function(y) sum((y[2]*exp(-gap*(y[3]^2)*coancs$drange)+y[1]-relcurve[i,j,])^2)
	  fit<-optim(x[i,j,], mf, method=optmethod)
	}
	if (ndates==2) { # number of dates per paired mixture
	  mf=function(y) sum((y[2]*exp(-gap*(y[3]^2)*coancs$drange)+y[4]*exp(-gap*(y[5]^2)*coancs$drange)+y[1]-relcurve[i,j,])^2)
	  fit<-optim(x[i,j,], mf, method=optmethod)
      }
	x[i,j,]=fit$par # b*exp(-lambda*d)+a w/ x=c(a,b,rate)
      }
      if (asym)
      {
	if (ndates==1) {
	  mf=function(y) sum((y[1]*exp(-gap*(y[2]^2)*coancs$drange)+1-relcurve[i,j,])^2)
	  fit<-optim(x[i,j,][2:3], mf, method=optmethod)
	  x[i,j,][2:3]=fit$par # a*exp(-lambda*d)+b w/ x=c(a,b,rate)
	}
	if (ndates==2) {
	  mf=function(y) sum((y[1]*exp(-gap*(y[2]^2)*coancs$drange)+y[3]*exp(-gap*(y[4]^2)*coancs$drange)-relcurve[i,j,])^2)
	  fit<-optim(x[i,j,][2:5], mf, method=optmethod)
	  x[i,j,][2:5]=fit$par # a*exp(-lambda*d)+b w/ x=c(a,b,rate)
	}
      }
    }
  }
  if (samedates)
  {
    x[,,3]=mean(x[,,3])
    if (!asym)
    {
      mf=function(y) 
      {
	ans=0;
	l=1
	for (i in 1:lpop) for (j in i:lpop) 
	{
	  ans=ans+sum((y[l+1]*exp(-gap*(y[length(y)]^2)*coancs$drange)+y[l]-relcurve[i,j,])^2)
	  l=l+2
	}
	ans
      }
      tmp=NULL;for (i in 1:lpop) for (j in i:lpop) tmp=c(tmp,x[i,j,1:2]);tmp=c(tmp,x[1,1,3])
      fit<-optim(tmp, mf, method=optmethod) 
      l=1
      for (i in 1:lpop) for (j in i:lpop) 
      {
	x[i,j,1]=x[j,i,1]=fit$par[l]
	for (dd in 1:ndates)
  	  x[i,j,2*dd]=x[j,i,2*dd]=fit$par[l+1]
	l=l+2
      }
      x[,,3]=fit$par[length(fit$par)]
    }
    if (asym) 
    {
      mf=function(y) 
      {
	ans=0;
	l=1
	for (i in 1:lpop) for (j in i:lpop) 
	{
	  ans=ans+sum((y[l]*exp(-gap*(y[length(y)]^2)*coancs$drange)+1-relcurve[i,j,])^2)
	  l=l+1
	}
	ans
      }
      tmp=NULL;for (i in 1:lpop) for (j in i:lpop) tmp=c(tmp,x[i,j,2]);tmp=c(tmp,x[1,1,3])
      fit<-optim(tmp, mf, method=optmethod)
      l=1
      for (i in 1:lpop) for (j in i:lpop) 
      {
	x[i,j,2]=x[j,i,2]=fit$par[l]
	l=l+1
      }
      x[,,3]=fit$par[length(fit$par)]
    }
  }
  for (i in 1:lpop) for (j in 1:lpop)
  {
    if (verbose) cat(i,j, " after:", x[i,j,], "\n")
    for (dd in 1:ndates)
      x[i,j,][1+dd*2]=x[i,j,][1+dd*2]^2
    params[i,j,]=x[i,j,]
  }
  for (i in 1:lpop) for (j in i:lpop)
  {
    if (PLOT)
    {
      YLIM=range(relcurve[i,j,])
      if (plotall & axisall)
	YLIM=c(min(relcurve[i,j,],coancs$relprobs[i,j,keep[[(i-1)*lpop+j]],]/kweights[[(i-1)*lpop+j]][keep[[(i-1)*lpop+j]]]),
	       max(relcurve[i,j,],coancs$relprobs[i,j,keep[[(i-1)*lpop+j]],]/kweights[[(i-1)*lpop+j]][keep[[(i-1)*lpop+j]]]))
      if (YLIM[1]>1) YLIM[1]=1 # always include 1
      if (YLIM[2]<1) YLIM[2]=1 # always include 1
      if (ndates==1)
        plot(range(gap*coancs$drange*100),YLIM,t='n',xlab="cM",ylab=ylab,main=paste0(targetname, popnames[i],":",popnames[j]," (",round(x[i,j,][3],1),")"))
      if (ndates==2)
        plot(range(gap*coancs$drange*100),YLIM,t='n',xlab="cM",ylab=ylab,main=paste0(targetname, popnames[i],":",popnames[j],
										     " (",round(x[i,j,][3],1),",",round(x[i,j,][5],1),")"))
      if (is.null(k) & plotall) # if plotting group consensus, also show hap / inds in grey
      {
	# don't forget to undo the weighting
	for (kk in which(keep[[(i-1)*lpop+j]]))
	  lines(gap*coancs$drange*100,coancs$relprobs[i,j,kk,]/kweights[[(i-1)*lpop+j]][kk],lwd=lwd/2,col=rgb(190/255,190/255,190/255,kweights[[(i-1)*lpop+j]][kk]^(1-transalpha))) 
	#lines(gap*coancs$drange*100,coancs$relprobs[i,j,kk,]/kweights[[(i-1)*lpop+j]][kk],lwd=lwd/2,col=kk) # debug line
      }
      lines(gap*coancs$drange*100,relcurve[i,j,],t='l',lwd=lwd)
      abline(h=1,col=2)
      if (ndates==1)
        lines(gap*coancs$drange*100 , x[i,j,][2]*exp(-gap*x[i,j,][3]*coancs$drange)+x[i,j,][1], col=3, lwd=lwd) # note x[3] is already exponentiated
      if (ndates==2)
        lines(gap*coancs$drange*100 , x[i,j,][2]*exp(-gap*x[i,j,][3]*coancs$drange)+x[i,j,][4]*exp(-gap*x[i,j,5]*coancs$drange)+x[i,j,][1], col=3, lwd=lwd) 
    }
  }
  return(list(params=params, relcurve=relcurve, gens.matrix=params[,,c(1+2*(1:ndates))], kweights=kweights, keep=keep)) # redundancy here but useful to focus 
}


bootstrap_chromosomes_coanc_curves=function(coancs,gap,localanc,alpha,nsamps=100,min.cM=1,max.cM=50,asym=FALSE,samedates=FALSE,optmethod="BFGS",thresh=1e-4)
{
  NUMA=dim(localanc[[1]])[2]
  NUMI=NUMA/2
  A=dim(localanc[[1]])[1]
  nchrno=length(localanc)
  kgens=rep(NaN,NUMI)
  for (k in 1:NUMI) if (min(alpha[[k]])>thresh) kgens[k]=mean(plot_coanccurves(coancs,gap,k=k,PLOT=FALSE,samedates=samedates,asym=asym,min.cM=min.cM)$params[,,3],na.rm=TRUE)
  G=sapply(localanc,function(x) dim(x)[3])
  boot.gens=list()
  boot.localanc=list()
  for (t.ch in 1:nchrno) 
    boot.localanc[[t.ch]]=array(NaN,c(A,NUMA,G[t.ch]))
  pb<-txtProgressBar(min=0,max=nsamps,style=3)
  for (r in 1:nsamps) ## nsamps bootstrap samples
  {
    setTxtProgressBar(pb, r)
    #bootstrap chromosomes; generate NUMI pseudo-individuals using random chromosomes drawn from all inds w/ replacement; see GT S4.4 for details
    boot.inds=matrix(sample(1:NUMI,NUMI*nchrno,replace=TRUE),NUMI)
    boot.haps=matrix(NaN,NUMA,nchrno);for (t.ch in 1:nchrno) for (ind in 1:NUMI) boot.haps[c(ind*2-1,ind*2),t.ch]=c(boot.inds[ind,t.ch]*2-1,boot.inds[ind,t.ch]*2)
    for (t.ch in 1:nchrno) 
      for (l in 1:A) for (hap in 1:NUMA)
        boot.localanc[[t.ch]][l,hap,]=localanc[[t.ch]][l,boot.haps[hap,t.ch],]
    boot.coancs=create_coancs(boot.localanc,gap,"DIP",max.cM=max.cM)#*mean(unlist(lambda))/100);
    boot.gens[[r]]=plot_coanccurves(boot.coancs,gap,PLOT=FALSE,samedates=samedates,asym=asym,min.cM=min.cM)$gens.matrix
  }  
  close(pb)
  return(list(kgens=kgens,boot.gens=boot.gens))
}

# this function bootstraps over individuals rather than chromosomes in pseudo-individuals
bootstrap_individuals_coanc_curves=function(coancs,gap,alpha,nsamps=100,min.cM=1,max.cM=50,asym=FALSE,samedates=FALSE,optmethod="BFGS",thresh=1e-4)
{
  NUMI=dim(coancs$ancprobs)[2]
  kgens=rep(NaN,NUMI)
  for (k in 1:NUMI) if (min(alpha[[k]])>thresh) kgens[k]=mean(plot_coanccurves(coancs,gap,k=k,PLOT=FALSE,samedates=samedates,asym=asym,min.cM=min.cM)$params[,,3],na.rm=TRUE)
  boot.gens=list()
  pb<-txtProgressBar(min=0,max=nsamps,style=3)
  for (r in 1:nsamps) ## nsamps bootstrap samples
  {
    setTxtProgressBar(pb, r)
    boot.inds=sample(1:NUMI,replace=TRUE)
    tmp.coancs=coancs;tmp.coancs$ancprobs=coancs$ancprobs[,boot.inds];tmp.coancs$relprobs=coancs$relprobs[,,boot.inds,]
    boot.gens[[r]]=plot_coanccurves(tmp.coancs,gap,PLOT=FALSE,samedates=samedates,asym=asym,min.cM=min.cM)$gens.matrix
  }  
  close(pb)
  return(list(kgens=kgens,boot.gens=boot.gens))
}




