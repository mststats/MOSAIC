# Fst calculations computational bottleneck is calculating the frequencies of markers for each population (panel).
# This function will perform these as a one off calculation for each panel and save to file for reuse.
source("fst.R")
write_panel_summaries=function(pathout="FREQS/",datasource="HGDP/", chrnos=1:22) {
  panels<-read.table(paste(datasource,"sample.names",sep=""), header=F)
  panels<-as.character(unique(panels[,1]))
  for (panel in panels)
  {
    cat("Looking at ", panel, " \n")
    pdata=summarise_panels(panel, datasource, chrnos)
    save(pdata, file=paste0(pathout, panel, "_", "freqs.rdata"))
  }
  return(NULL)
}
