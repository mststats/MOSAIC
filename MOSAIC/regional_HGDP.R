# as per Table S15 of GlobeTrotter paper
Americas=c("Colombian","Karitiana","Maya","Pima","Surui")
Africa=c("BantuKenya","BantuSouthAfrica","BiakaPygmy","Egyptian","Ethiopian","EthiopianJew","Hadza","Mandenka","MbutiPygmy","Sandawe","SanNamibia","SanKhomani","Yoruba")
SouthEastAsia=c("Cambodian","Dai","Han","HanNchina","Japanese","Lahu","Miao","Myanmar","Naxi","She","Tu","Tujia","Xibo","Yi")
Oceania=c("Melanesian","Papuan")

# these are the regional analyses
MiddleEastNorthAfrica=list("recipients"=c("Bedouin","Egyptian","Iranian","Jordanian","Moroccan","Mozabite","Palestinian","Syrian","Tunisian","UAE"),
		 "masked"=c(Americas,Oceania))
Ethiopian=list("recipients"=c("Ethiopian","EthiopianJew"),
		 "masked"=c(Americas,Oceania))
Mediterranean=list("recipients"=c("EastSicilian","Greek","Sardinian","SouthItalian","Spanish","WestSicilian"),
		 "masked"=c(Americas,Oceania,MiddleEastNorthAfrica$recipients))
CentralAsia=list("recipients"=c("Balochi","Brahui","Burusho","Hazara","Kalash","Makrani","Pathan","Sindhi","Uygur","Uzbekistani"),
		 "masked"=c(Americas,Oceania))
San=list("recipients"=c("SanKhomani","SanNamibia"),
		 "masked"=c(Americas,Oceania))
EastEuropeI=list("recipients"=c("Belorussian","Bulgarian","Chuvash","Finnish","Hungarian","Lithuanian","Polish","Romanian","Russian"),
		 "masked"=c(Americas,Oceania,Africa,SouthEastAsia))
EastEuropeII=list("recipients"=c("Belorussian","Bulgarian","Chuvash","Finnish","Hungarian","Lithuanian","Romanian","Russian"),
		 "masked"=c(Americas,Oceania,Africa,SouthEastAsia))

file.create(file="regional.txt")
regions=list(MiddleEastNorthAfrica,Ethiopian,CentralAsia,San,EastEuropeI,EastEuropeII)
names(regions)=c("MiddleEastNorthAfrica","Ethiopian","CentralAsia","San","EastEuropeI","EastEuropeII")

for (r in 1:length(regions))
  for (i in 1:length(regions[[r]]$recipients))
    write(c(names(regions)[r],regions[[r]]$recipients[i], regions[[r]]$recipients[-i], regions[[r]]$masked), file="regional.txt", append=T, ncol=length(unlist(regions[[r]]))+1)
