# given unique donor and unique target haplotype lists (from compressed_grid.cpp) calculate the number of matches between them along the genome
r_create_umatch=function(ch.d.w,ch.t.w,g.map,t.G)
{
  ans=list()
  for (g in 1:t.G)
  {
    gm=which(g.map==g)
    if (length(gm)>0)
    {
      tmp=matrix(0,length(ch.d.w$u[[g]]), length(ch.t.w$u[[g]]))
      for (i in 1:nrow(tmp)) for (j in 1:ncol(tmp)) tmp[i,j]=sum(ch.d.w$u[[g]][[i]]==ch.t.w$u[[g]][[j]],na.rm=T) # note the na.rm=T for missing data compatibility
      ans[[g]]=tmp
    }
    if (length(gm)==0)
      ans[[g]]=matrix(0)
  }
  ans
}
require(compiler)
create_umatch=cmpfun(r_create_umatch,list(optimize=3))

