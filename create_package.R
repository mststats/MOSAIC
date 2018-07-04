require(Rcpp)

Rcpp.package.skeleton(name="MOSAIC",example_code=F,cpp_files=dir(pattern="*.cpp"),force=T,
		      maintainer="michael.salter-townshend@ucd.ie",author="Michael Salter-Townshend",
		      code_files=c("grid.R","compressed_grid.R","calc_r2.R","transitions.R","ancunaware.R","coancestry.R",
				   "plot_funcs.R","fst.R","admix.R","donates.R","create_logfile.R","init_Mu.R","intermediate_calcs.R",
				   "mix_hmm.R","cleanup.R","phase_funcs.R","initProb.R","EM_updates.R","all_donates.R","noanc.R","klikelihood.R",
				   "plot_localanc.R","log_funcs.R","summarise_panels.R","read_panels.R","localanc.R","setup.R","run.R"))

