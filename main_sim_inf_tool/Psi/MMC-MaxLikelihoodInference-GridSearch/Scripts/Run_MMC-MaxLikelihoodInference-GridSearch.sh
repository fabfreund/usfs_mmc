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

_psiVector=( 0 0.15 0.3 0.45 0.6 0.75 0.9 )
_rhoVector=( 0 1 10 100 )
_noSegSitesVector=( 100 1000 10000 )
_sampleSizeVector=( 20 50 100 200 )
_lumpVector=( 0 5 10 15 )

for (( m=1; m<=${#_lumpVector[@]}; m++ ))
	do
	for (( i=1; i<=${#_psiVector[@]}; i++ ))
		do
		for (( j=1; j<=${#_rhoVector[@]}; j++ ))
			do
			for (( n=1; n<=${#_noSegSitesVector[@]}; n++ ))
				do
				for (( k=1; k<=${#_sampleSizeVector[@]}; k++ ))
				do
					_Psi=${_psiVector[$((i - 1))]}
					_Rho=${_rhoVector[$((j - 1))]}
					_noSegSites=${_noSegSitesVector[$((n - 1))]}
					_noSamples=${_sampleSizeVector[$((k - 1))]}
					_lump=${_lumpVector[$((m - 1))]}

					_minPsi=0
					_maxPsi=1
					_noStepsPsi=100

					_minRho=0
					_maxRho=1024
					_noStepsRho=1024

					_specPath=noSamples-$_noSamples/Psi-$_Psi/Rho-$_Rho/segSites-$_noSegSites
					_fileNameSFS="SFS-Psi-$_Psi-rho-$_Rho-noSamples-$_noSamples-segSites-$_noSegSites.txt"
					_fileNamePhi=MinPsi-$_minPsi-MaxPsi-$_maxPsi-NoStepPsi-$_noStepsPsi-MinRho-$_minRho-MaxRho-$_maxRho-NoStepRho-$_noStepsRho-noSamples-$_noSamples\_expectedPhi.txt

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

					if [ $(echo "$_minPsi == 0" | bc -l) == 1 ]
					then
						_minPsi=""
					else
						_minPsi="-minPsi $_minPsi"
					fi

					if [ $(echo "$_maxPsi == 0" | bc -l) == 1 ]
					then
						_maxPsi=""
					else
						_maxPsi="-maxPsi $_maxPsi"
					fi

					if [ $_noStepsPsi == 0 ]
					then
						_noStepsPsi=""
					else
						_noStepsPsi="-noStepsPsi $_noStepsPsi"
					fi

					if [ $(echo "$_minRho == 0" | bc -l) == 1 ]
					then
						_minRho=""
					else
						_minRho="-minRho $_minRho"
					fi

					if [ $(echo "$_maxRho == 0" | bc -l) == 1 ]
					then
						_maxRho=""
					else
						_maxRho="-maxRho $_maxRho"
					fi

					if [ $_noStepsRho == 0 ]
					then
						_noStepsRho=""
					else
						_noStepsRho="-noStepsRho $_noStepsRho"
					fi

					#***************
					#	EXECUTION
					#***************

					$_nameOfProgram -SFS $_pathToDataSFSFile/$_specPath/$_fileNameSFS -Phi $_pathToDataPhiFile/$_fileNamePhi $_lump $_folded $_printGrid $_minPsi $_maxPsi $_noStepsPsi $_minRho $_maxRho $_noStepsRho

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
