#!/bin/bash
#$ -cwd -N mosaic_run
#$ -P myers.prjc -q short.qc -j yes -V -pe shmem 8
MC=8 # don't forget to match to shmem above
target=SpainPopn_2
NUMA=16
datasource=spanish/
L=2
firstind=1
echo $target $L-way started at: `date`
Rscript run.R $target $datasource $L $firstind $NUMA $MC > LOGS/"$target"_"$firstind"_"$L".out 2> LOGS/"$target"_"$firstind"_"$L".error 
echo $target $L-way finished at: `date`
