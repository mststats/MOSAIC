tract_lengths=function(t.localanc, thresh=0.8) {
  ans=list()
  t.A=dim(t.localanc[[1]])[1]
  for (a in 1:t.A) {
    ans[[a]]=list()
    for (ch in 1:length(t.localanc)) {
      ans[[a]][[ch]]=list()
      anc.calls=t(apply(t.localanc[[ch]][a,,],2,function(x) x>thresh))
      for (i in 1:dim(localanc[[1]])[2]) {
	d=diff(anc.calls[,i])
	segs=cumsum(abs(d))
	ans[[a]][[ch]][[i]]=as.numeric(table(segs))
      }
    }
  }
  return(ans)
}

bin_tracts=function(tract.lengths, nh=20, min.length=NULL, max.length=NULL, mids=NULL) {
  if (is.null(min.length)) min.length=min(unlist(tract.lengths))
  if (is.null(max.length)) max.length=max(unlist(tract.lengths))
  if (is.null(mids)) mids=seq(min.length,max.length,length.out=nh)
  ans=data.frame(mids=mids)
  t.A=length(tract.lengths)
  for (a in 1:t.A)
    ans[,a+1]=sapply(mids,function(x) sum(x==mids[sapply(unlist(tract.lengths[[a]]),function(y) which.min(abs(mids-y)))]))
  colnames(ans)=c("mids",1:t.A)
  return(ans)
}

plot_tract_lengths=function(binned.tracts, t.scale=100, xlab="tract length (cM)",ylab="#tracts") {
  t.A=ncol(binned.tracts)-1 # first column is breakpoints
  plot(t.scale*binned.tracts$mids[c(1,nrow(binned.tracts))], c(1,max(binned.tracts[,-1])), log="y",xlab=xlab,ylab=ylab,t="n") # blank plot
  for (a in 1:t.A) 
    points(t.scale*binned.tracts$mids,binned.tracts[,1+a],pch=20,col=a,t='b') # add one line per ancestry
}
