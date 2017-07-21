source("ancunaware.R")
for (ind in 1:NUMI) for (ch in 1:nchrno) flips[[ind]][[ch]][]=F # this required if not reloading but working with current session OR if RPE>0
getnoancgfbs=T;o.LOG=F;PLOT=F;get_switches=F;source("noanc.R");Mu=a.Mu;rho=a.rho;theta=a.theta;Q=a.Q;lambda=a.lambda;alpha=a.alpha
unphased_localanc=get_ancunaware_localanc()  # works off noanc_gfbs
source("coancestry.R");if (verbose) cat("calculating coancestry curves\n"); coancs=create_coancs(unphased_localanc,dr,"DIP")

source("cleanup.R")
save(file=paste0(resultsdir,target,"_", L, "way_", firstind, "-", firstind+NUMI-1, "_", paste(chrnos[c(1,nchrno)],collapse="-"),"_",NN,"_",GpcM,"_",a.prop.don,"_",a.max.donors,".RData"), 
     target, phase.error.locs, o.Mu, o.lambda, o.theta, o.alpha, o.Q, o.rho, Mu, lambda, theta, alpha, Q, rho, L, NUMA, nchrno, chrnos, g.loc, tol, dr, NL, kLL, RPE, acoancs, coancs)
