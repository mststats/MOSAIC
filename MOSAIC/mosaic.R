# run EM algorithm for total iterations
#if (verbose) pb<-txtProgressBar(min=1,max=ITER,style=3)
for (ITER in 1:total)
{
  old.Mu<-Mu; old.PI<-PI; old.lambda<-lambda; old.alpha<-alpha; old.rho<-rho; old.theta<-theta
  old.cloglike<-cloglike
  source("EM_updates.R") 
  source("initProb.R")
  source("klikelihood.R") # E-step: extra work here as fors will be calculated next iteration of E.n above
  cat(round(100*ITER/total), "%: ", cloglike, "(", cloglike-old.cloglike, ")", "\n")
  if (!is.na(old.cloglike)) 
  {
    if ((old.cloglike - cloglike)>1e-3)
    {
      Mu<-old.Mu; PI<-old.PI; lambda<-old.lambda; alpha<-old.alpha; rho<-old.rho; theta<-old.theta
      source("initProb.R")
      cloglike<-old.cloglike
      warning("loglikelihood has decreased; abandoning EM", immediate.=T)
      break
    }
    if ((cloglike - old.cloglike)< eps) 
      {cat("EM iterations have converged\n");break;}
  }
  if (LOG) 
  {
    runtime<-as.numeric(Sys.time());diff.time<-runtime-old.runtime;old.runtime<-runtime;
    writelog("EM")
  }
  #if (verbose) setTxtProgressBar(pb, m)
}
#if (verbose) close(pb)
