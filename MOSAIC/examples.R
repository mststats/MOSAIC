if (target!="simulated")
{
  if (datasource=="HGDP/" & target=="Hazara") 
    pops<-list(EA=c("HanNchina", "Mongola", "Daur", "Oroqen", "Hezhen", "Xibo", "Tu", "Yi"),#, "Uzbekistani", "BantuKenya", "Mozabite"),
               NE=c("Pathan", "Iranian", "Balochi", "Makrani", "Cypriot", "Armenian", "Turkish", "Polish"),#, "Chuvash", "Mandenka"), 
               AD=target)
  if (datasource=="HGDP/" & target=="NorthAfrican") 
    pops<-list(European=c("French", "English", "Scottish", "Welsh", "GermanyAustria", "Ireland", 
			  "Greek", "NorthItalian", "Turkish", "Spanish"), #"Sardinian", "Tuscan", "WestSicilian", "EastSicilian", "SouthItalian", "Cypriot", 
	       African=c("Yoruba", "BantuKenya", "BantuSouthAfrica", "Mandenka", "BiakaPygmy", "Hadza", "MbutiPygmy"), #"Sandawe", 
	       AD=target)
  if (datasource=="spanish/") 
    pops<-list(EU=c("France_17", "Germany-Belgium_1", "Germany-Hungary_7", "Ireland_8", "Neatherlands-Germany_27", "Poland_19", "Switzerland_11", "Switzerland_14"),
	       `NA`=c("Egypt_85", "Libya-Algeria_82", "NorthAfrica.M-A-L_89", "NorthMorocco_92", "SS.Kenya-LWK_40", "Sub-saharan.Nigeria_43", "Tunisia_83", "WesternSahara_90"),   
	       AD=target)
  if (datasource=="chromopainter/") 
    pops<-list(panels=panels,AD=target)
  panels=unlist(pops)
  kLL=length(panels)-1
}

if (target=="simulated")
{
  require(gtools) # for rdirichlet
  sim.alpha<-sim.lambda<-list()
  NUMI<-max(1,NUMA/2)
  if (!exists("mean.sim.alpha")) mean.sim.alpha=c(rdirichlet(1, rep(8,L)))
  #mean.sim.alpha<-c(0.9,0.1) # manual choice for L=2
  #mean.sim.alpha<-c(0.1,0.1,0.8) # manual choice for L=3
  for (ind in 1:NUMI)
  {
    sim.alpha[[ind]]=c(rdirichlet(1, mean.sim.alpha*L*4))
    #sim.alpha[[ind]]=c(0.2,0.8) # leads to difficulties
    sim.lambda[[ind]]<-o.lambda
  }
  #sim.alpha[[1]]=c(1,0,0);sim.alpha[[2]]=c(0,1,0)#;sim.alpha[[3]]=c(0,0,1); # one pure from each anc=pop; useful when theta and rho not updated
  #sim.alpha[[1]]=c(1,1,1);sim.alpha[[2]]=c(1,1,1); # two that are equal parts each anc=pop; useful when theta and rho not updated
  if (!exists("fewer_ancs")) fewer_ancs=NULL
  if (!is.null(fewer_ancs)) 
  {
    for (ind in fewer_ancs)
    {
      remanc=sample(1:L)[1:sample(1:(L-1),1)] # sample a subset of the L to remove from this individual, but never all of them
      sim.alpha[[ind]][remanc]=0
      if (L>2)
      {
	# never have ancs 1 and 2 both present
	if (sim.alpha[[ind]][1]) sim.alpha[[ind]][2]=0; if (sim.alpha[[ind]][2]) sim.alpha[[ind]][1]=0 
      }
      sim.alpha[[ind]]=sim.alpha[[ind]]/sum(sim.alpha[[ind]])
    } # remove first anc from first ind
    #for (ind in (1:NUMI)[-fewer_ancs]) # make those with this first anc only have a little
    #  {sim.alpha[[ind]][1]=0.1;sim.alpha[[ind]]=sim.alpha[[ind]]/sum(sim.alpha[[ind]])} # remove first anc from first ind
  }
  if (!exists("prop.missing")) prop.missing=0
  # takes in multipanels and simulates admixed individuals using the first panel in each anc
  if (!is.null(ANC)) # if set to TRUE then use one of these examples, depending on choice of L
  {
    if (ANC[1]==T) # if just set to TRUE but mixing panels not supplied directly
    {
      if (L==2)
      {
	# pops<-list(Eu=c("Brahui", "Balochi", "Burusho", "Armenian", "Iranian", "Lezgin","Georgian"),
	#            As=c("Hezhen", "HanNchina", "Han", "Tu", "Oroqen", "Daur", "Xibo", "Tujia", "Yakut"))
	# 
	#mixers=c("Sindhi","Mongola")
	pops<-list(Eu=c("Spanish","Norwegian","Polish","English","Welsh","Ireland","Scottish"),
		   Af=c("Mandenka","Hadza","BantuKenya","BantuSouthAfrica","Sandawe","BiakaPygmy","MbutiPygmy"))
	mixers=c("French","Yoruba")
      }
      if (L==3)
      {
	pops<-list(Eu=c("English", "Scottish", "Spanish", "Polish", "Melanesian", "Turkish", "Norwegian", "Hungarian"),
		   #As=c("Han", "Yi", "Daur", "She", "Maya"),# "HanNchina", "Hezhen", "Oroqen", "Tu", # easier
		   SEA=c("Sindhi", "Indian", "Makrani", "IndianJew", "Yemeni"), # harder 
		   Af=c("Mandenka", "BantuKenya", "Sandawe", "BantuSouthAfrica", "Pima", "Surui", "MbutiPygmy", "Hadza")
		   )
	mixers=c("French", "Pathan", "Yoruba") # harder
	#mixers=c("French", "Japanese", "Yoruba") # easier
      }
      if (L==4)
      {
	pops<-list(Eu=c("English", "Scottish", "Polish", "Spanish"),
		   As=c("Han", "Yi", "Daur", "She"),
		   Af=c("Mandenka","BantuKenya", "Sandawe"),
		   Sa=c("Maya", "Surui"))
	mixers=c("French", "Japanese", "Yoruba", "Pima")
      }
    }
    if (ANC[1]!=T) # something supplied
    { 
      mixers=ANC
      # check that the supplied panels to mix are indeed members of panels and that length(ANC)==L
      refs=sample(panels[!(panels%in%mixers)]);tmp=sample(1:L, length(refs), replace=T); # random order of unused panels, then split into approx equal groups
      pops=list()
      for (l in 1:L) pops[[l]]=c(ANC[l],refs[tmp==l]) # first panel is always the panel used to simulate the admixed individuals, remaining panels taken at random
    }
  }
  if (is.null(ANC))
  {
    ANC=mixers=sample(panels,L)
    if (verbose)
      cat("no ancestral panels specified so using:", ANC, "\n")
    refs=sample(panels[!(panels%in%mixers)]);tmp=sample(1:L, length(refs), replace=T); # random order of unused panels, then split into approx equal groups
    pops=list()
    for (l in 1:L) pops[[l]]=c(ANC[l],refs[tmp==l]) # # first panel is always the panel used to simulate the admixed individuals, remaining panels taken at random
  }
  pops[[L+1]]="simulated";names(pops)[L+1]="target"
  panels<-unlist(pops)
  kLL=length(panels)-1 # remove simulated as there is no such panel
  panels=panels[1:kLL]
  if (verbose) cat("Creating ", NUMI, " simulated ", L, "-way admixed target individuals with ", kLL, " panels\n", sep="")
  LL=kLL+1
  ANC<-NULL
  for (i in 1:L)
    ANC[i]=mixers[i] 
}
