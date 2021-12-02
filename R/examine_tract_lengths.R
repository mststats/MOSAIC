tract_lengths=function(t.localanc, thresh=0.8) {
  f.lengths=list()
  t.A=dim(t.localanc[[1]])[1]
  for (a in 1:t.A) {
    f.lengths[[a]]=list()
    for (ch in 1:length(t.localanc)) {
      f.lengths[[a]][[ch]]=list()
      anc.calls=t(apply(t.localanc[[ch]][a,,],2,function(x) x>thresh))
      for (i in 1:dim(localanc[[1]])[2]) {
	d=diff(anc.calls[,i])
	segs=cumsum(abs(d))
	f.lengths[[a]][[ch]][[i]]=as.numeric(table(segs))
      }
    }
  }
  return(f.lengths)
}
plot_tract_lengths=function (tract_lengths, t.dr, nh=20) {
  t.A=length(tract_lengths)
  h=list()
  for (a in 1:t.A) 
    h[[a]]=hist(unlist(f.lengths[[a]])*t.dr*100,nh) # multiplying by dr uses the MOSAIC grid-length and 100 to convert to cM
  plot(range(sapply(1:t.A, function(a) range(h[[a]]$mids))),1+range(sapply(1:t.A, function(a) range(h[[a]]$counts))),log="y",xlab="cM",ylab="#tracts",t="n") # blank plot
  for (a in 1:t.A) 
    points(h[[a]]$mids,h[[a]]$counts,pch=20,col=a,t='b') # add one line per ancestry
}
