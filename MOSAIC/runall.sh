#targets="Hazara Chuvash Kalash Russian SanKhomani Maya Pima Spanish SanNamibia Colombian Adygei Armenian Balochi BantuKenya BantuSouthAfrica Basque Bedouin 
#targets="Belorussian BiakaPygmy Brahui Bulgarian Burusho Cambodian 
#targets="Cypriot Dai Daur Druze EastSicilian Egyptian English Ethiopian EthiopianJew Finnish French Georgian 
# GermanyAustria Greek Hadza Han HanNchina Hezhen Hungarian IndianJew Indian Iranian Ireland Japanese Jordanian Karitiana Lahu Lezgin Lithuanian 
# Makrani Mandenka MbutiPygmy Melanesian Miao Mongola Moroccan Mozabite Myanmar Naxi NorthItalian Norwegian Orcadian 
#Oroqen Palestinian Papuan Pathan Polish Romanian Sandawe Sardinian Saudi Scottish She Sindhi SouthItalian Surui Syrian Tu Tujia Tunisian Turkish 
# Tuscan UAE Uygur Uzbekistani Welsh WestSicilian Xibo Yakut Yemeni Yi Yoruba"
targets="Burusho Mozabite Russian Yakut French Japanese Sandawe Sardinian SanKhomani Han Spanish Druze Bedouin Palestinian"

firstind=1
datasource="HGDP/";NUMA=92
for L in {2..3}
do
  for target in $targets
  do
    Rscript run.R $target $datasource $L $firstind $NUMA > LOGS/"$target"_"$firstind"_"$L".out 2> LOGS/"$target"_"$firstind"_"$L".error 
  done 
done 

#datasource=HGDP_PEL;NUMA=32
#L=3;targets="PEL" 
#for target in $targets
#do
# Rscript run.R $target $datasource $L $firstind $NUMA > LOGS/"$target"_"$firstind"_"$L".out 2> LOGS/"$target"_"$firstind"_"$L".error 
#done

#"MALAWI" "MOSSI" "SEMI-BANTU" "BANTU" "MANDINKA" "FULA" "JOLA" "WOLLOF" "MANJAGO" "SERERE" "SEREHULE" "NAMKAM" "KASEM" "AKANS" "CHONYI" "GIRIAMA" "KAUMA" "KAMBE" "GBR" "FIN" "CHS" "CDX" "IBS" "PEL" "KHV" "CEU" "YRI" "CHB" "JPT" "LWK" "TSI" "GIH" "MKK" "MALINKE" "BAMBARA" "WABONDEI" "WASAMBAA" "MZIGUA" "KHOMANI" "KARRETJIE" "KHWE" "GUIGHANAKGAL" "JUHOANSI" "NAMA" "XUN" "SEBANTU" "SWBANTU" "AMAXHOSA" "JUHOAN" "XUNV" "TYGRAY" "SUDANESE" "OROMO" "WOLAYTA" "AFAR" "ESOMALI" "ARIBLACKSMITHIBD" "GUMUZ" "SOMALI" "ARICULTIVATOR" "ARIBLACKSMITH" "ANUAK" "AMHARA" "ARICULTIVATORIBD" "GUMUZIBD" "SUDANESEIBD" "AMHARAIBD"

