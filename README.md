This repository contains analysis and simulation pipelines as well as the analysis result tables/data objects from our recent preprint

"Interpreting the pervasive observation of U-shaped Site Frequency Spectra"

available on biorxiv: https://www.biorxiv.org/content/10.1101/2022.04.12.488084v2

Please cite this paper if you intend to work with our code. 

This repository is still in beta state, if there are questions/comments please write to fafreund_popgen*at*gmx.net (replace *at* with @) 

The subfolders contain:

 * datafilter_main: An example script how we extracted SNPs, the SFS and reconstructed the (putative) ancestral alleles from an alignment of ingroup sample and an outgroup sample
 * example_fasta: Contains an example alignment (fasta format, see subfolder README for the data sources)
 * sfs_data: Contains the SFS data and related data from which we performed the main text analyses (code and results in main_dataanalysis/, data sources are referenced in the pdf within the folder)
 * main_dataanalysis: Contains the full main text data analysis from our paper, including all results  
 * main_sim_inf_tool: This contains the main tool (written in C++, needs to be compiled before usage, see README within the folder)   
 * sim_study: Full pipeline for the simulation study part from the paper (we simulated SFS data under different coalescent models and then ran our pseudolikelihood inference method on it to see how well we retrieve the true model and its parameters)
 * supp_analysis: Contains all supplementary analyses scripts 
	* population structure scans
	* sliding windows nucleotide diversity for the data samples
	* comparison of ingroup and between-group nucleotide diversity for sample and outgroup to estimate bias through incomplete lineage sorting (ILS)
	* adjusted pseudolikelihood model selection and parameter estimation when taking into account transition/transversion information for ancestral allele misidentification
         
