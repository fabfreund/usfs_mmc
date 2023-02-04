The folder contains an example of our supplementary data analysis - tailored to one relatively small data set (provided in the folder ../example_fasta/). 
Other data sets were analysed analogously, with one exception this only meant different input data formats (e.g. vcf instead of fasta) and/or bioinformatic adjustments (running scripts on chunks of data to circumvent memory overflow, using more instance of parallel computing). 
The exception is the sliding windows measurement of nucleotide diversity \pi, which was performed directly with vcftools for vcf input files.

