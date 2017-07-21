r_create_umatch=function(ch.d.w,ch.t.w,g.map,t.G)
{
  ans=list()
  for (g in 1:t.G)
  {
    gm=which(g.map==g)
    if (length(gm)>0)
    {
      tmp=matrix(0,length(ch.d.w$u[[g]]), length(ch.t.w$u[[g]]))
      for (i in 1:nrow(tmp)) for (j in 1:ncol(tmp)) tmp[i,j]=sum(ch.d.w$u[[g]][[i]]==ch.t.w$u[[g]][[j]],na.rm=T)  
      ans[[g]]=tmp
    }
    if (length(gm)==0)
      ans[[g]]=matrix(0)
  }
  ans
}
require(compiler)
create_umatch=cmpfun(r_create_umatch,list(optimize=3))

#d.w=list() # map to unique donor  haps at each gridpoint
#t.w=list() # map to unique target haps at each gridpoint
#umatch=list() # lookup for donor,target to #matches using t.w and d.w
#
#for (ch in 1:nchrno)
#{
#  snps<-read.table(paste0(datasource,"snpfile.",chrnos[ch])) 
#  #### inputs #################################
#  all_rates<-matrix(scan(paste0(datasource,"rates.",chrnos[ch]),skip=1,quiet=T),ncol=2)
#  locs<-as.integer(snps[,4])
#  tmp=match(locs, all_rates[,1])
#  rates=all_rates[tmp,2] # use ones with hap data; some may be nmissing if in snps file but not in rates file
#  # rates are flat in sections so use rate to the left if missing
#  for (l in which(is.na(tmp))) rates[l]=all_rates[which.max(all_rates[all_rates[,1]<locs[l],1]-locs[l]),2]
#  rates<-rates/100 # /100 to move to morgans from centimorgans 
#  rm(all_rates)
#  g.rates<-seq(rates[1],rates[S[ch]],l=G[ch])
#  g.map<-tapply(1:S[ch], 1:S[ch], function(s) which.min((rates[s]-g.rates)^2)) # create map from rates to grid
#  d.w[[ch]]=list(u=list(),w=list());k=1
#  for (k in 1:NUMP) d.w[[ch]]=cpp_unique_haps(Y[k,],k,S[ch],G[ch],g.map-1,max(table(g.map)),d.w[[ch]])
#  t.w[[ch]]=list(u=list(),w=list())
#  for (k in 1:NUMA) t.w[[ch]]=cpp_unique_haps(Y[k+NUMP,],k,S[ch],G[ch],g.map-1,max(table(g.map)),t.w[[ch]])
#  umatch[[ch]]=create_umatch(d.w[[ch]],t.w[[ch]],g.map,G[ch])
#  # don't need to keep the haps any more
#  d.w[[ch]]$u=NULL
#  t.w[[ch]]$u=NULL
#}
#
