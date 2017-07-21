# of course often returns far more than needed
#wg=as.integer(0.5/mean(rho)) # half the expected chunk width
#wg=5;thresh.misses=5
max.donors=NUMP-1
source("all_donates.R")
f_max=f_misses(max.donors,ch,ind,umatch[[ch]],d.w[[ch]]$w,t.w[[ch]]$w)
f_mid=f_misses(100,ch,ind,umatch[[ch]],d.w[[ch]]$w,t.w[[ch]]$w)
f_min=f_misses(1,ch,ind,umatch[[ch]],d.w[[ch]]$w,t.w[[ch]]$w)
plot(density(f_min),xlim=c(0,max(f_max)),col=3);lines(density(f_mid));lines(density(f_max),col=2)
legend("topright", c("best donor", "100th donor", "worst donor"), lty=rep(1,3), col=c(3,1,2))
