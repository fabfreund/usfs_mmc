This folder contains the data analysis from the main manuscript. It runs the main inference tools (from the folder main_sim_inf_tool/) which
optimize over the pseudolikelihood (Eq. 2 in the manuscript) of a grid of parameter values (coalescent model parameter alpha/Psi, growth parameter g, misorientation parameter e). From this, it then performs the testing procedure based on pseudolikelihood ratios/BFs (Eq. 3).

There are usually three variants of scripts and folder, indicated by names including beta, psi or smallpsi. This indicates whether we infer under the Beta coalescent, the coarser Psi coalescent grid or the finder Psi coalescent grid which only covers small Psi values (indicated by smallpsi)   

 * The folder sfs_data (so not a subfolder) contains all input data, the SFS of all data sets considered
 * The optimisation scripts are ..._dataanalysis.R. These scripts are running a loop analysing all data sets. 
 * We have precomputed the expected SFS for all grid values from the manuscript (subfolders phi_tables_...)
 * The output of the main inference tools are a) the full likelihood grid for all parameter combinations (subfolders likelihood_grids_...) and, for convenience, also b) the ML estimates (subfolders ML_est_...)
 * From this output and the expected SFS grids for all three grids, one can compute the model selection test and the Cramer's V values by running post_check_gofmeasures.R
 * Its output is in subfolder res_summaries/, both as a semicolon delimited text file and a .RData object. This is the full result table for the main analysis!
 * Then, we can generate the supplemental figures 1 and 2 by running fits_pic.R. For each species, the figures show SFS, neutralized SFS and the expected SFS under best fitting parameters for all three grid and model choices.
 
