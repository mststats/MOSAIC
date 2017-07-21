#source("fast_prethin_donates.R")
require(parallel)
if (verbose)
  cat("\npre thinning donors to below", thresh.prop, "matches: ")
if (NUMA==1) {H=1} else H=2
prethin_ndonors<-prethin_donates<-prethin_donatesl<-prethin_donatesr<-list()
if (HPC!=2)
{
  for (ch in 1:nchrno)
  {
    prethin_ndonors[[ch]]<-list()
    prethin_donates[[ch]]<-prethin_donatesl[[ch]]<-prethin_donatesr[[ch]]<-list()
    # if using all then a vector will suffice for each prethin_donates, prethin_donatesl, prethin_donatesr rather than a matrix with G[ch] columns
    NvecsG=ifelse(max.donors==NUMP, 1, G[ch]) 
    if (HPC==1)
    {
      tmp<-foreach(ind=1:NUMI) %dopar%
      {
        tmp2=pre_create_donates(ch,ind,umatch[[ch]],d.w[[ch]]$w,t.w[[ch]]$w,thresh.prop)
	ans_ndonors=tmp2$ndonors
	ans_donates=ff(tmp2$donates,vmode="integer",dim=c(NUMP,NvecsG),filename=paste0(ffpath,target,"_prethin_donates_",ch,"_",ind,".ff"),overwrite=T);close(ans_donates)
	ans_donatesl=ff(tmp2$donatesl,vmode="integer",dim=c(NUMP,NvecsG),filename=paste0(ffpath,target,"_prethin_donatesl_",ch,"_",ind,".ff"),overwrite=T);close(ans_donatesl)
	ans_donatesr=ff(tmp2$donatesr,vmode="integer",dim=c(NUMP,NvecsG),filename=paste0(ffpath,target,"_prethin_donatesr_",ch,"_",ind,".ff"),overwrite=T);close(ans_donatesr)
	rm(tmp2)
	gc()
	list(ndonors=ans_ndonors,donates=ans_donates,donatesl=ans_donatesl,donatesr=ans_donatesr)
      }
    }
    if (!HPC)
    {
      tmp<-foreach(ind=1:NUMI) %dopar%
      pre_create_donates(ch,ind,umatch[[ch]],d.w[[ch]]$w,t.w[[ch]]$w,thresh.prop)
    }
    for (ind in 1:NUMI)
    {
      prethin_ndonors[[ch]][[ind]]<-tmp[[ind]]$ndonors
      prethin_donates[[ch]][[ind]]<-tmp[[ind]]$donates
      prethin_donatesl[[ch]][[ind]]<-tmp[[ind]]$donatesl
      prethin_donatesr[[ch]][[ind]]<-tmp[[ind]]$donatesr
    }
    rm(tmp)
  }
}
if (HPC==2)
{
  tmp<-foreach(ch_ind=1:(nchrno*NUMI)) %dopar%
  {
    ch=as.integer((ch_ind-0.5)/NUMI)+1
    ind=(ch_ind-1)%%NUMI+1
    NvecsG=ifelse(max.donors==NUMP, 1, G[ch]) 
    tmp2=pre_create_donates(ch,ind,umatch[[ch]],d.w[[ch]]$w,t.w[[ch]]$w,thresh.prop)
    ans_ndonors=tmp2$ndonors
    ans_donates=ff(tmp2$donates,vmode="integer",dim=c(NUMP,NvecsG),filename=paste0(ffpath,target,"_prethin_donates_",ch,"_",ind,".ff"),overwrite=T);close(ans_donates)
    ans_donatesl=ff(tmp2$donatesl,vmode="integer",dim=c(NUMP,NvecsG),filename=paste0(ffpath,target,"_prethin_donatesl_",ch,"_",ind,".ff"),overwrite=T);close(ans_donatesl)
    ans_donatesr=ff(tmp2$donatesr,vmode="integer",dim=c(NUMP,NvecsG),filename=paste0(ffpath,target,"_prethin_donatesr_",ch,"_",ind,".ff"),overwrite=T);close(ans_donatesr)
    rm(tmp2)
    gc()
    list(ndonors=ans_ndonors,donates=ans_donates,donatesl=ans_donatesl,donatesr=ans_donatesr)
  }
  for (ch in 1:nchrno)
  {
    prethin_ndonors[[ch]]=prethin_donates[[ch]]=prethin_donatesl[[ch]]=prethin_donatesr[[ch]]=list()
    for (ind in 1:NUMI)
    {
      ch_ind=(ch-1)*NUMI+ind
      prethin_ndonors[[ch]][[ind]]<-tmp[[ch_ind]]$ndonors
      prethin_donates[[ch]][[ind]]<-tmp[[ch_ind]]$donates
      prethin_donatesl[[ch]][[ind]]<-tmp[[ch_ind]]$donatesl
      prethin_donatesr[[ch]][[ind]]<-tmp[[ch_ind]]$donatesr
    }
  }
  rm(tmp)
}
if (verbose)
  cat("average of", mean(sapply(prethin_ndonors,sapply,mean)), "out of", NUMP, "donors per gridpoint\n")

