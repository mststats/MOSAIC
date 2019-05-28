# Fst calculations computational bottleneck is calculating the frequencies of markers for each population (panel).

# functions to read in panel data and summarise as freqs and counts
summarise_panels=function(panelname, pathin, chrnos)
{
  allp=list();alln=list() 
  for (ch in 1:length(chrnos))
  {
    snps=read.table(paste0(pathin,"snpfile.",chrnos[ch]))
    S=nrow(snps)
    tmp<-scan(paste0(pathin,panelname,"genofile.",chrnos[ch]),what="character",quiet=T) 
    tmp<-strsplit(tmp,"")
    tmpy=matrix(sapply(tmp, as.double), ncol=S)
    ny=nrow(tmpy)
    allp[[ch]]=colMeans(tmpy)
    alln[[ch]]=ny 
  }
  return(list("freqs"=allp,"counts"=alln))
}


# This function will perform these as a one off calculation for each panel and save to file for reuse.
# if panels=NULL then it will run for all panels in sample.names
write_panel_summaries=function(pathout="FREQS/",datasource="HGDP/", chrnos=1:22, panels=NULL) {
  if (!file.exists(pathout))
    dir.create(file.path(getwd(), pathout))
  if (is.null(panels)) 
  {
    panels<-read.table(paste(datasource,"sample.names",sep=""), header=F)
    panels<-as.character(unique(panels[,1]))
  }
  for (panel in panels)
  {
    cat("Saving SNP frequencies of", panel, " \n")
    pdata=summarise_panels(panel, datasource, chrnos) 
    save(pdata, file=paste0(pathout, panel, "_", "freqs.rdata"))
  }
  return(NULL)
}

# similar function that uses estimated local ancestry along the admixed genome
write_admixed_summary=function(target,NL,pathout="FREQS/",datasource="HGDP/",targetdatasource="MOSAIC_RESULTS/",g.loc,t.localanc,chrnos=1:22)
{
  A=dim(t.localanc[[1]])[1]
  if (!file.exists(pathout))
    dir.create(file.path(getwd(), pathout))
  cat("Saving frequencies of simulated data\n")
  ancestral_freqs=maximal_alleles(target,chrnos,g.loc,t.localanc,datasource,targetdatasource) 
  save(ancestral_freqs, file=paste0(pathout, target, "_", A, "way_", sum(NL), "_freqs.rdata"))
  return(NULL)
}
