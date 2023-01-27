The two subfolders contain three subfolders each

 * MMC-CoalescentSimulator/: Allows to simulate data under a specific MMC coalescent model with exponential growth and ancestral allele misorientation
 * MMC-Phi-LookUpTable/: Essentially computes the expected site frequency spectrum (SFS) for a grid of equidistant choices of MMC coalescent parameter and exponential growth parameter. 
	* "Essentially" means that the executable actually computes the expected sum of lengths of all the branches from the genealogical tree on which mutations cause exactly 1,2,...,n-1 sequences to inherit this mutation (these are the E(T_i)'s from Eq. 2 in the manuscript )
        * The output is a data table/matrix, with each row containing information for one parameter combination   
        * The first column records the total expected length E(T_tot) of the genealogical tree (again see Eq. 2 from the manuscript), and the following columns show E(T_i)/E(T_tot) (which equals the expected ith entry of the SFS divided by the total number of segregating sites) 
 * MMC-MaxLikelihoodInference-GridSearch: Computes and maximises the pseudolikelihood function (Eq.2 in the manuscript) across a expected SFS grid (generated via the executable from  MMC-Phi-LookUpTable/)

To get executables, simply run make in each of the three sub-subfolders.

Then, run the executables with --help to get a short documentation on how to use them 
