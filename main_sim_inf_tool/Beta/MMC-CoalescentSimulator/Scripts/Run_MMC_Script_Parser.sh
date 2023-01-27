#!/bin/sh

#************************* Run_MMC_Script ********************************************
# This script runs the MMC program and passes the command line arguments.
# For further details please refer to the manual.
#*************************************************************************************


_pathToData="/Path/To/Data"
_pathToRename="/Path/To/Rename"

_alphaVector=( 1.25 1.5 1.75 2)			# alpha (sweepstake parameter)
_rhoVector=( 0 1 10 100 )				# Exponential growth parameter
_initialPopSizeVector=( 10000 )			# Population size (when samples were taken)
_muVector=( 0.0001 )					# Population-wide mutation rate
_noSamplesVector=( 20 50 100 200 )		# Sample size

_pathToProgram="/Path/To/Program"				# Pass the path to the compiled MMC program
_nameOfProgram="MMC-CoalescentSimulator.out"																		# Pass the the name to the compiled MMC program

_noReplicates=10000				# Number of simulations
_inference=1

_SFS=0
_sequences=0
_cumulBranchLength=1
_relBranchLength=0
_rates=0
_sumStats=0
_trees=0
_seed=-1
_moran=0
_newick=""
_outputPath=""


	#********************
	#	PARSING INPUT
	#********************

	if [ $_SFS == 1 ]
	then
		_SFS="-SFS"
	else
		_SFS=""
	fi

	if [ $_inference == 1 ]
	then
		_inference="-inference"
	else
		_inference=""
	fi

	if [ $_sequences == 1 ]
	then
		_sequences="-seq"
	else
		_sequences=""
	fi

	if [ $_cumulBranchLength == 1 ]
	then
		_cumulBranchLength="-cumulBrL"
	else
		_cumulBranchLength=""
	fi

	if [ $_relBranchLength == 1 ]
	then
		_relBranchLength="-relBrL"
	else
		_relBranchLength=""
	fi

	if [ $_rates == 1 ]
	then
		_rates="-rates"
	else
		_rates=""
	fi

	if [ $_sumStats == 1 ]
	then
		_sumStats="-sumStats"
	else
		_sumStats=""
	fi

	if [ $_trees == 1 ]
	then
		_trees="-trees"
	else
		_trees=""
	fi

	if [ $_moran == 1 ]
	then
		_moran="-moran"
	else
		_moran=""
	fi

	if [ -n "$_newick" ]
	then
		_newick="-newick $_newick"
	fi

	if [ -n "$_outputPath" ]
	then
		_outputPathFinal=$_outputPath
		_outputPath="-out $_outputPath/"
	else
		_outputPath="-out $_pathToProgram/"
		_outputPathFinal=$_pathToProgram
	fi

	if [ $_seed != -1 ]
	then
		_seed="-seed $_seed"
	else
		_seed=""
	fi


for (( i=1; i<=${#_alphaVector[@]}; i++ ))
	do
	_alpha=${_alphaVector[$((i - 1))]}									# Alpha (sweepstake parameter)

	for (( j=1; j<=${#_rhoVector[@]}; j++ ))
		do
		_rho=${_rhoVector[$((j - 1))]}								# Exponential growth parameter
		for (( n=1; n<=${#_noSamplesVector[@]}; n++ ))
		do
		_noSamples=${_noSamplesVector[$((n - 1))]}
			for (( k=1; k<=${#_initialPopSizeVector[@]}; k++ ))
				do
				_initialPopSize=${_initialPopSizeVector[$((k - 1))]}	# Population size (when samples were taken)
				for (( m=1; m<=${#_muVector[@]}; m++ ))
					do
					_mu=${_muVector[$((m - 1))]}					# Population-wide mutation rate

					#***************
					#	EXECUTION
					#***************

					$_pathToProgram/$_nameOfProgram $_inference -s $_noSamples -alpha $_alpha -mu $_mu -num $_noReplicates -N $_initialPopSize -rho $_rho $_SFS $_sumStats $_sequences $_rates $_cumulBranchLength $_relBranchLength $_trees $_moran $_newick $_seed $_outputPath

					#***********************
					#	PROCESSING OUTPUT
					#***********************

					perl $_pathToRename/rename.pl 's/(.*_)\d*_\d*_\d*_\d*_(.*txt)/$1$2/g' $_outputPathFinal/*.txt
					perl $_pathToRename/rename.pl 's/(.*_)\d*_\d*_\d*_\d*_(.*tree)/$1$2/g' $_outputPathFinal/*.tree

					eval "mkdir -p $_pathToData/CoalescentDataNewBrLength/noSamples-$_noSamples/Alpha-$_alpha/Rho-$_rho/relativeBranchLengths"
					mv $_outputPathFinal/*Coalescent*.txt $_pathToData/CoalescentDataBrLength/noSamples-$_noSamples/Alpha-$_alpha/Rho-$_rho/relativeBranchLengths

					_fileName=Alpha-$_alpha-Rho-$_rho-noSamples-$_noSamples-branchLengths.txt
					tail -n+11 $_pathToData/CoalescentDataBrLength/noSamples-$_noSamples/Alpha-$_alpha/Rho-$_rho/relativeBranchLengths/*.txt > $_pathToData/CoalescentDataBrLength/noSamples-$_noSamples/Alpha-$_alpha/Rho-$_rho/relativeBranchLengths/$_fileName
					rm $_pathToData/CoalescentDataBrLength/noSamples-$_noSamples/Alpha-$_alpha/Rho-$_rho/relativeBranchLengths/*cumulativeBranchLength.txt
					done
				done
			done
		done
	done

exit 0
