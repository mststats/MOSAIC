HapMixMODE="BOTH"
PLOTMODE="BAR"
OUTHAPMIX=T
# run MOSAIC
ptm=proc.time()
source("run.R")
mosaic.time=proc.time() - ptm

# run HapMix equivalent
setwd("../../HapmixReleasev2/")
EM=T;TRUTH=T
ptm=proc.time()
runcode=paste0("bash run_mosaic_example.sh ", chrnos[1], " ", chrnos[length(chrnos)], " '", HapMixMODE, "' ", NUMA, " ", EM)
system(runcode)
hapmix.time=proc.time() - ptm

cat("MOSAIC took ", mosaic.time[3], "seconds \n")
cat("HapMix took ", hapmix.time[3], "seconds \n")

# compare the results
setwd("../Rainbow/sandbox")
cexa=1.5;PNG=F
source("read_hapmix_results.R")
save.image("hapalso.RData")
#source("plot_hapmix_compare.R")

