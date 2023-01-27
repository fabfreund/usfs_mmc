#!/bin/sh

#************************* Run_MMC-Phi-LookUpTable_Script *********************************
#  This script runs the MMC-Phi-LookUpTable program and passes the command line arguments.
#  For further details please refer to the manual.
#******************************************************************************************


_pathToOutput="/Users/matu/Documents/MMC-Marcel/MultipleMerger-Project/Data/PhiTables/PsiCoalescent"

_pathToProgram="/Users/matu/Documents/MMC-Marcel/MultipleMerger-Project/Psi/MMC-Phi-LookUpTable"
_nameOfProgram="MMC-Phi-LookUpTable.out"				# Pass the name of the compiled MMC program


_minPsi=0.74
_maxPsi=0.78
_noStepsPsi=4

_minRho=0
_maxRho=100
_noStepsRho=0

_sampleSize=79
_precisionPhi=150


	#********************
	#	PARSING INPUT
	#********************


	if [ $(echo "$_minPsi == 0" | bc -l) == 1 ]
	then
		_minPsi=""
	else
		_minPsi="-MinPsi $_minPsi"
	fi

	if [ $(echo "$_maxPsi == 0" | bc -l) == 1 ]
	then
		_maxPsi=""
	else
		_maxPsi="-MaxPsi $_maxPsi"
	fi

	if [ $_noStepsPsi == -1 ]
	then
		_noStepsPsi=""
	else
		_noStepsPsi="-NoStepsPsi $_noStepsPsi"
	fi

	if [ $(echo "$_minRho == 0" | bc -l) == 1 ]
	then
		_minRho=""
	else
		_minRho="-MinRho $_minRho"
	fi

	if [ $(echo "$_maxRho == 0" | bc -l) == 1 ]
	then
		_maxRho=""
	else
		_maxRho="-MaxRho $_maxRho"
	fi

	if [ $_noStepsRho == 0 ]
	then
		_noStepsRho=""
	else
		_noStepsRho="-NoStepsRho $_noStepsRho"
	fi

	if [ $_sampleSize == 0 ]
	then
		_sampleSize=""
	else
		_sampleSize="-sampleSize $_sampleSize"
	fi

	if [ $(echo "$_precisionPhi == 0" | bc -l) == 1 ]
	then
		_precisionPhi=""
	else
		_precisionPhi="-precPhi $_precisionPhi"
	fi

	if [ -n "$_pathToOutput" ]
	then
		_pathToOutput="-OUT $_pathToOutput"
	else
		_pathToOutput=""
	fi

	#***************
	#	EXECUTION
	#***************

	$_pathToProgram/$_nameOfProgram $_pathToOutput/ $_minPsi $_maxPsi $_noStepsPsi $_minRho $_maxRho $_noStepsRho $_sampleSize $_precisionPhi


exit 0
