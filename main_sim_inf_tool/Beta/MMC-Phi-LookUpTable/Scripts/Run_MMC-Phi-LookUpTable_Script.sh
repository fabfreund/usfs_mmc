#!/bin/sh

#************************* Run_MMC-Phi-LookUpTable_Script *********************************
#  This script runs the MMC-Phi-LookUpTable program and passes the command line arguments.
#  For further details please refer to the manual.
#******************************************************************************************


_pathToOutput="/Path/To/Output"

_pathToProgram="/Path/To/Program"
_nameOfProgram="MMC-Phi-LookUpTable.out"				# Pass the name of the compiled MMC program


_minAlpha=0
_maxAlpha=0
_noStepsAlpha=-1

_minRho=0
_maxRho=0
_noStepsRho=0

_sampleSize=0
_precisionPhi=0


	#********************
	#	PARSING INPUT
	#********************


	if [ $(echo "$_minAlpha == 0" | bc -l) == 1 ]
	then
		_minAlpha=""
	else
		_minAlpha="-MinAlpha $_minAlpha"
	fi

	if [ $(echo "$_maxAlpha == 0" | bc -l) == 1 ]
	then
		_maxAlpha=""
	else
		_maxAlpha="-MaxAlpha $_maxAlpha"
	fi

	if [ $_noStepsAlpha == -1 ]
	then
		_noStepsAlpha=""
	else
		_noStepsAlpha="-NoStepsAlpha $_noStepsAlpha"
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

	$_pathToProgram/$_nameOfProgram $_pathToOutput/ $_minAlpha $_maxAlpha $_noStepsAlpha $_minRho $_maxRho $_noStepsRho $_sampleSize $_precisionPhi


exit 0
