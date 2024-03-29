# Fst calculations computational bottleneck is calculating the frequencies of markers for each population (panel).

# functions to read in panel data and summarise as freqs and counts
summarise_panels=function(panelname, pathin, chrnos)
{
  allp=list();alln=list() 
  for (ch in 1:length(chrnos))
  {
    tmpfilename=file.path(pathin,paste0(panelname,"genofile.",chrnos[ch]))
    tmp<-scan(tmpfilename,what="character",quiet=TRUE,nlines=1)
    N2<-nchar(tmp)
    tmpy=matrix(suppressWarnings(as.integer(as.matrix(laf_open_fwf(tmpfilename, column_widths=rep(1,N2),column_types=rep("character",N2))[,]))),ncol=N2)
    allp[[ch]]=rowMeans(tmpy,na.rm=TRUE)
    alln[[ch]]=N2
  }
  return(list("freqs"=allp,"counts"=alln))
}


# This function will perform these as a one off calculation for each panel and save to file for reuse.
# if panels=NULL then it will run for all panels in sample.names
write_panel_summaries=function(pathout="MOSAIC_RESULTS/FREQS", datasource, chrnos=1:22, panels=NULL) {
  if (!dir.exists(pathout))
    dir.create(file.path(pathout))
  if (is.null(panels)) 
  {
    panels<-read.table(file.path(datasource,"sample.names"), header=FALSE)
    panels<-as.character(unique(panels[,1]))
  }
  for (panel in panels)
  {
    cat("Saving SNP frequencies of", panel, " \n")
    pdata=summarise_panels(panel, datasource, chrnos) 
    save(pdata, file=file.path(pathout, paste0(panel, "_", "freqs.rdata")))
  }
  return(NULL)
}

# similar function that uses estimated local ancestry along the admixed genome
write_admixed_summary=function(target,NL,pathout="MOSAIC_RESULTS/FREQS",datasource,targetdatasource="MOSAIC_RESULTS",g.loc,t.localanc,chrnos=1:22)
{
  A=dim(t.localanc[[1]])[1]
  if (!dir.exists(pathout))
    dir.create(file.path(pathout))
  cat("Saving frequencies of simulated data\n")
  ancestral_freqs=maximal_alleles(target,chrnos,g.loc,t.localanc,datasource,targetdatasource) 
  save(ancestral_freqs, file=file.path(pathout, paste0(target, "_", A, "way_", sum(NL), "_freqs.rdata")))
  return(NULL)
}
