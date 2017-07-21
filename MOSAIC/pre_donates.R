# of course often returns far more than needed
require(compiler)
r_f_misses=function(donor,ch,ind,subset_match,dw,tw,wg=20) 
{
  ans=rep(0,G[ch])
  tmpg=1:(wg+1)
  haps=c(ind*2-1,ind*2)
  tmp=rep(0,length(tmpg))
  for (g1 in 1:length(tmp))
    tmp[g1]=2*gobs[[ch]][[ind]][tmpg[g1]]-sum(subset_match[[tmpg[g1]]][dw[[tmpg[g1]]][donor],tw[[tmpg[g1]]][haps]+1]) 
  g=1
  tmpleft=1;tmpright=1:wg
  ans[g]=max(sum(tmp[tmpleft])/length(tmpleft),sum(tmp[tmpright])/length(tmpright)) # max of mean over chunk left and mean over chunk right
  for (g in 2:(wg+1))
  {
    # add one point to the right
    tmp=c(tmp,2*gobs[[ch]][[ind]][g+wg]-sum(subset_match[[g+wg]][dw[[g+wg]][donor],tw[[g+wg]][haps]+1]))
    tmpleft=1:g;tmpright=wg:length(tmp)
    ans[g]=max(sum(tmp[tmpleft])/length(tmpleft),sum(tmp[tmpright])/length(tmpright)) # max of mean over chunk left and mean over chunk right
  }
  tmpleft=1:wg;tmpright=wg:(2*wg)
  for (g in g:(G[ch]-wg))
  {
    #remove leftmost point and add one to the right
    tmp=c(tmp[-1],2*gobs[[ch]][[ind]][g+wg]-sum(subset_match[[g+wg]][dw[[g+wg]][donor],tw[[g+wg]][haps]+1])) 
    #tmpleft=1:wg;tmpright=wg:(2*wg) # stays the same here
    ans[g]=max(sum(tmp[tmpleft])/length(tmpleft),sum(tmp[tmpright])/length(tmpright)) # max of mean over chunk left and mean over chunk right
  }
  for (g in (G[ch]-wg+1):G[ch])
  {
    #remove leftmost point only
    tmp=tmp[-1]
    tmpleft=1:wg;tmpright=wg:length(tmp)
    ans[g]=max(sum(tmp[tmpleft])/(wg+1),sum(tmp[tmpright])/(wg+1)) # max of mean over chunk left and mean over chunk right
  }
  ans
}
f_misses=cmpfun(r_f_misses)


r_pre_create_donates<-function(ch,ind,subset_match,dw,tw,thresh=thresh.prop) 
{ 
  nmisses=matrix(0L,G[ch],NUMP)
  for (donor in 1:NUMP)
    nmisses[,donor]=f_misses(donor,ch,ind,subset_match,dw,tw)
  t.ndonors=apply(nmisses,1,function(x) max(min.donors,sum(x<thresh))) # FLAG change to proportion of max matches rather than #misses
  t.donates=apply(nmisses,1,order)
  t.donatesl<-t.donatesr<-matrix(0L,NUMP,G[ch]) # re-sizing and 0s first
  t.donatesl[,2:G[ch]]<-matrix(unlist(tapply(2:G[ch], 2:G[ch], function(g) match(t.donates[,g], t.donates[,g-1])),use.names=F),NUMP)
  t.donatesr[,1:(G[ch]-1)]<-matrix(unlist(tapply(1:(G[ch]-1), 1:(G[ch]-1), function(g) match(t.donates[,g], t.donates[,g+1])),use.names=F),NUMP)
  t.donatesl[is.na(t.donatesl)]<-0 # match() returns NA where no match occurs, replace with 0s
  t.donatesr[is.na(t.donatesr)]<-0 # match() returns NA where no match occurs, replace with 0s  
  return(list(ndonors=t.ndonors,donates=t.donates-1,donatesl=t.donatesl-1,donatesr=t.donatesr-1)) # -1 to convert to cpp indexing
 }
pre_create_donates=cmpfun(r_pre_create_donates,list(optimize=optlevel))


