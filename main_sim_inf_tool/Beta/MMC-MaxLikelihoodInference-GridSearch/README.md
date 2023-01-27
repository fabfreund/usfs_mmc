# README #

This repository contains all files related to the Maximum Likelihood Estimation Program for estimating 'psi' (i.e., the multiple merger (MMC) parameter) and/or the exponential growth rate. 
There are two different folders: MaxLikelihood-MMC and Scripts.
The main simulation program, written in C++, is contained in MaxLikelihood-MMC. It can be easily compiled using the provided Makefile (by simply typing 'make'). Note that compilation requires the Gnu Scientific Library, the GMP Library, the MPFR library, and the Boost Library to be installed.
Scripts contains a shell script to run the Maximum Likelihood Estimation Program.

### Wishlist ###

* Write own log-binomial function for coalescent rates to increase numerical stability (analogous to HybridLambda) ::DONE BUT CURRENTLY UNUSED::
* Check how to integrate with R [RCPP](http://www.rcpp.org/)
* Make inits member of class. Should probably change for each data set (treated as fixed now).
* Add more explanation (comments)
* Add output file path as command line argument
* Provide pre-compiled binaries for multiple platforms
* Add lambda/beta coalescent (this requires some mathematical proofs though)
* change int to unsigned int
* efficient and stable way to compute harmonics numbers
* best place to compute the harmonic numbers and other constants for summary stats to avoid redundancies
* check

### Who do I talk to? ###

* If you have questions or comments please contact sebastian.matuszewski[ÄT]epfl.ch
