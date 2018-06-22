# Fst calculations computational bottleneck is calculating the frequencies of markers for each population (panel).
# This function will perform these as a one off calculation for each panel and save to file for reuse.
# if panels=NULL then it will run for all panels in sample.names
write_panel_summaries=function(pathout="FREQS/",datasource="HGDP/", chrnos=1:22, panels=NULL) {
  if (is.null(panels)) 
  {
    panels<-read.table(paste(datasource,"sample.names",sep=""), header=F)
    panels<-as.character(unique(panels[,1]))
  }
  for (panel in panels)
  {
    cat("Looking at ", panel, " \n")
    pdata=summarise_panels(panel, datasource, chrnos)
    save(pdata, file=paste0(pathout, panel, "_", "freqs.rdata"))
  }
  return(NULL)
}

# similar function that uses estimated local ancestry along the admixed genome
# after loading results from a MOSAIC run target, NL, etc are in global memory
write_admixed_summary=function(pathout="FREQS/",datasource="HGDP/",t.localanc,chrnos=1:22)
{
  ancestral_freqs=maximal_alleles(target,chrnos,t.localanc,datasource,datasource) 
  save(ancestral_freqs, file=paste0(pathout, target, "_", L, "way_", sum(NL), "_freqs.rdata"))
  return(NULL)
}
