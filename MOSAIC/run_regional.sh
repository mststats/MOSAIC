L=2
firstind=1
datasource="HGDP/";NUMA=100
MC=16
Rscript regional_HGDP.R # create the regional targets and maskings
K=`wc -l regional.txt | cut -d " " -f 1`
for k in $(eval echo {1.."$K"})
do
 Rscript run.R $k $datasource $L $firstind $NUMA $MC > LOGS/"$k"_"$firstind"_"$L".out 2> LOGS/"$k"_"$firstind"_"$L".error 
done 

