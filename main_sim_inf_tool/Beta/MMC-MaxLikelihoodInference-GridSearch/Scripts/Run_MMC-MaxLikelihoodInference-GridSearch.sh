#!/bin/sh

#********************************** Run_MMC-MaxLikelihoodInference-GridSearch *************************************
#  This script runs the MMC-MaxLikelihoodInference-GridSearch program and passes the command line arguments.
#  For further details please refer to the manual.
#******************************************************************************************************************


_pathToRename="/Path/To/Rename"
_nameOfProgram="/Path/To/Program/MMC-MaxLikelihoodInference-GridSearch.out"				# Pass the path (including the name) to the compiled MMC program

_pathToDataSFSFile="/Path/To/SFS"
_pathToDataPhiFile="/Path/To/Phi"

_pathToResultsFile="/Path/To/Results"

_alphaVector=( 1 1.25 1.5 1.75 2 )
_rhoVector=( 0 1 10 100 )
_noSegSitesVector=( 100 1000 10000 )
_sampleSizeVector=( 20 50 100 200 )
_lumpVector=( 0 5 10 15 )

_misroeintation

for (( m=1; m<=${#_lumpVector[@]}; m++ ))
	do
	for (( i=1; i<=${#_alphaVector[@]}; i++ ))
		do
		for (( j=1; j<=${#_rhoVector[@]}; j++ ))
			do
			for (( n=1; n<=${#_noSegSitesVector[@]}; n++ ))
				do
				for (( k=1; k<=${#_sampleSizeVector[@]}; k++ ))
				do
					_Alpha=${_alphaVector[$((i - 1))]}
					_Rho=${_rhoVector[$((j - 1))]}
					_noSegSites=${_noSegSitesVector[$((n - 1))]}
					_noSamples=${_sampleSizeVector[$((k - 1))]}
					_lump=${_lumpVector[$((m - 1))]}

					_minAlpha=1
					_maxAlpha=2
					_noStepsAlpha=100

					_minRho=0
					_maxRho=1024
					_noStepsRho=1024

					_minMisIdent=0
					_maxMisIdent=0.15
					_noStepsMisIdent=15

					_specPath=noSamples-$_noSamples/Alpha-$_Alpha/Rho-$_Rho/segSites-$_noSegSites
					_fileNameSFS="SFS-Alpha-$_Alpha-rho-$_Rho-noSamples-$_noSamples-segSites-$_noSegSites.txt"
					_fileNamePhi=MinAlpha-$_minAlpha-MaxAlpha-$_maxAlpha-NoStepAlpha-$_noStepsAlpha-MinRho-$_minRho-MaxRho-$_maxRho-NoStepRho-$_noStepsRho-noSamples-$_noSamples\_expectedPhi.txt

					_folded=0
					_printGrid=1

					#********************
					#	PARSING INPUT
					#********************

					if [ $_lump == 0 ]
					then
						_lump=""
					else
						_lump="-lump $_lump"
					fi

					if [ $_folded == 0 ]
					then
						_folded=""
					else
						_folded="-folded"
					fi

					if [ $_printGrid == 0 ]
					then
						_printGrid=""
					else
						_printGrid="-printGrid"
					fi

					if [ $(echo "$_minAlpha == -1" | bc -l) == 1 ]
					then
						_minAlpha=""
					else
						_minAlpha="-minAlpha $_minAlpha"
					fi

					if [ $(echo "$_maxAlpha == -1" | bc -l) == 1 ]
					then
						_maxAlpha=""
					else
						_maxAlpha="-maxAlpha $_maxAlpha"
					fi

					if [ $_noStepsAlpha == -1 ]
					then
						_noStepsAlpha=""
					else
						_noStepsAlpha="-noStepsAlpha $_noStepsAlpha"
					fi

					if [ $(echo "$_minRho == -1" | bc -l) == 1 ]
					then
						_minRho=""
					else
						_minRho="-minRho $_minRho"
					fi

					if [ $(echo "$_maxRho == -1" | bc -l) == 1 ]
					then
						_maxRho=""
					else
						_maxRho="-maxRho $_maxRho"
					fi

					if [ $_noStepsRho == -1 ]
					then
						_noStepsRho=""
					else
						_noStepsRho="-noStepsRho $_noStepsRho"
					fi


					if [ $(echo "$_minMisIdent == -1" | bc -l) == 1 ]
					then
						_minMisIdent=""
					else
						_minMisIdent="-minMisIdent $_minMisIdent"
					fi

					if [ $(echo "$_maxMisIdent == -1" | bc -l) == 1 ]
					then
						_maxMisIdent=""
					else
						_maxMisIdent="-maxMisIdent $_maxMisIdent"
					fi

					if [ $_noStepsMisIdent == -1 ]
					then
						_noStepsMisIdent=""
					else
						_noStepsMisIdent="-noStepsMisIdent $_noStepsMisIdent"
					fi

					#***************
					#	EXECUTION
					#***************

					$_nameOfProgram -SFS $_pathToDataSFSFile/$_specPath/$_fileNameSFS -Phi $_pathToDataPhiFile/$_fileNamePhi $_lump $_folded $_printGrid $_minAlpha $_maxAlpha $_noStepsAlpha $_minRho $_maxRho $_noStepsRho $_minMisIdent $_maxMisIdent $_noStepsMisIdent

					#***********************
					#	PROCESSING OUTPUT
					#***********************

					perl $_pathToRename/rename.pl 's/(.*_)\d*_\d*_\d*_\d*_(.*txt)/$1$2/g' $_pathToDataSFSFile/$_specPath/*.txt
					mv $_pathToDataSFSFile/$_specPath/*GridML*.txt $_pathToResultsFile/$_specPath/

					done
				done
			done
		done
	done
exit 0
