/*
 * MMC-MaxLikelihood-Grid is used to estimate the MMC parameter alpha and/or the
 * exponential growth parameter from SFS data.
 *
 * Copyright (C) 2016 Sebastian Matuszewski & Marcel Hildebrandt
 *
 * This file is part of MMC-MaxLikelihood-Grid.
 *
 * MMC-MaxLikelihood-Grid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "myIOParser.hpp"

#define MAX_DATE 20

#ifndef VERSION
	#define VERSION "v1.0"
#endif


#if defined _WIN32
	#define NAVIGATE "\\"
#endif

#if defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__))
	#define  NAVIGATE "/"
#endif


//*********************
//** Output Switches **
//*********************

char outfile_name[255];
std::string outfile_prefix;
std::string dateconst;

//**********************
//*** Output Streams ***
//**********************

std::ofstream _offMLEstimate;
std::ofstream _offL2AbsEstimate;
std::ofstream _offL1AbsEstimate;
std::ofstream _offGrid;

long basePrecision = std::cout.precision();
long extendedPrecision = 20;

////////////////
// parameters //
////////////////

int sampleSizeSFS = 1, sampleSizePhi = 1, noDataSets = 0;

double minAlpha = 1.0, maxAlpha = 2.0, minRho = 0.0, maxRho = 5.0, minMisIdent = 0, maxMisIdent = 0.15, misIdentSim = -1;
int noStepsAlpha = 100, noStepsRho = 100, noStepsMisIdent = 15;

int lumping = 0;
int noIncongruentSites = -1;
bool lhoodMisIdent = true;
bool folded = false;


std::vector<std::vector<double> > SFS;							// vector of SFS data sets
std::vector<std::vector<std::vector<double> > > Phi;		// vector of Phi data
std::vector<std::vector<double> > T_Total;					// vector of T_Total data
std::string SFSInputFile = "false";
std::string PhiInputFile = "false";
bool print_Grid = false;

void readcommandLineArguments(int argc, const char *argv[])
{
	if (argc == 1) // if input is incomplete stop
	{
		std::cout << "Error: Insufficient number of arguments passed.\n";
		printHelpMessage();
		exit(1);
	}
	
	int argc_i = 1;
	while (argc_i < argc)
	{
		std::string argv_i(argv[argc_i]);
		if ( argv_i == "-h" || argv_i == "-help" ){ printHelpMessage(); exit(0);}
		else if ( argv_i == "-minAlpha" ){ minAlpha = readNextInput<double>(argc, argc_i, argv); }
		else if ( argv_i == "-maxAlpha" ){ maxAlpha = readNextInput<double>(argc, argc_i, argv); }
		else if ( argv_i == "-minRho" ){ minRho = readNextInput<double>(argc, argc_i, argv); }
		else if ( argv_i == "-maxRho" ){ maxRho = readNextInput<double>(argc, argc_i, argv); }
		else if ( argv_i == "-noIncongruentSites" ){ noIncongruentSites = readNextInput<int>(argc, argc_i, argv); }
		else if ( argv_i == "-minMisIdent" ){ minMisIdent = readNextInput<double>(argc, argc_i, argv); }
		else if ( argv_i == "-maxMisIdent" ){ maxMisIdent = readNextInput<double>(argc, argc_i, argv); }
		else if ( argv_i == "-noStepsAlpha" ){ noStepsAlpha = readNextInput<int>(argc, argc_i, argv); }
		else if ( argv_i == "-noStepsRho" ){ noStepsRho = readNextInput<int>(argc, argc_i, argv); }
		else if ( argv_i == "-noStepsMisIdent" ){ noStepsMisIdent = readNextInput<int>(argc, argc_i, argv); }
		else if ( argv_i == "-lump" ){ lumping = readNextInput<int>(argc, argc_i, argv); }
		else if ( argv_i == "-SFS" ){ readNextStringto(SFSInputFile , argc_i, argc,  argv); }
		else if ( argv_i == "-Phi" ){ readNextStringto(PhiInputFile , argc_i, argc,  argv); }
		else if ( argv_i == "-folded" ){ folded = true; }
		else if ( argv_i == "-printGrid" ){ print_Grid = true; }
		
		else{ throw (std::string("Unknown flag: ") + argv_i); }
		argc_i++;
	}
	
	dateconst = get_date();
}

template < class T > T readNextInput(int argc, int& argc_i, const char *argv[])
{
	++argc_i;
	if (argc_i >= argc) throw (std::string( "Not enough parameters when parsing options: ") + argv[argc_i-1]);
	
	char c;
	T input;
	std::stringstream ss( argv[argc_i] );
	ss >> input;
	if (ss.fail() || ss.get(c)) throw ( std::string( "Failed to parse option: ") + argv[argc_i]);
	return input;
}

void readNextStringto(std::string &readto , int& argc_i, int argc, char const *argv[])
{
	argc_i++;
	if (argc_i >= argc) throw (std::string("Not enough parameters when parsing options: ") + argv[argc_i-1]);
	readto = std::string(argv[argc_i]);
	if ( readto[0] == '-' ) throw (std::string("Not enough parameters when parsing options: ") + argv[argc_i-1]);
}

void printHelpMessage()
{
	std::cout << std::endl << "MaxLikelihood-MMC-Grid " << VERSION << std::endl << std::endl;
	std::cout << "Usage:" << std::endl;
	std::cout << std::setw(20) << "-h/-help"         << "  --  " << "Help. List the following content." << std::endl;
	std::cout << std::setw(20) << "-minAlpha FLT"           << "  --  " << "Lower boundary for multiple merger parameter alpha. By default 1 is used." << std::endl;
	std::cout << std::setw(20) << "-maxAlpha FLT"           << "  --  " << "Upper boundary for multiple merger parameter alpha. By default 2 is used (corresponding to the Kingman case)." << std::endl;
	std::cout << std::setw(20) << "-minRho FLT"           << "  --  " << "Lower boundary for exponential growth rate rho. By default 0 is used (no growth)." << std::endl;
	std::cout << std::setw(20) << "-maxRho FLT"           << "  --  " << "Upper boundary for exponential growth rate rho. By default 5 is used." << std::endl;
	std::cout << std::setw(20) << "-noInCongruentSites INT"           << "  --  " << "Number of incongruent Sites that have been eliminated before calculating the SFS. By default -1 is used. This option implies that misidentification will not be estimated." << std::endl;
	std::cout << std::setw(20) << "-minMisIdent FLT"           << "  --  " << "Lower boundary for misidentification rate MisIdent. By default 0 is used (no growth)." << std::endl;
	std::cout << std::setw(20) << "-maxMisIdent FLT"           << "  --  " << "Upper boundary for misidentification rate MisIdent. By default 0.15 is used." << std::endl;
	std::cout << std::setw(20) << "-noStepsAlpha INT"           << "  --  " << "Number of equally spaced steps between minAlpha and maxAlpha in grid search algorithm. By default 100 is used." << std::endl;
	std::cout << std::setw(20) << "-noStepsRho INT"           << "  --  " << "Number of equally spaced steps between minRho and maxRho in grid search algorithm. By default 100 is used." << std::endl;
	std::cout << std::setw(20) << "-noStepsMisIdent INT"           << "  --  " << "Number of equally spaced steps between minMisIdent and maxMisIdent in grid search algorithm. By default 15 is used." << std::endl;
	std::cout << std::setw(20) << "-lump INT"           << "  --  " << "Frequency class cutoff above which all SFS classes are lumped into a single class. By default the unlumped SFS is used." << std::endl;
	std::cout << std::setw(20) << "-folded"           << "  --  " << "Whether SFS should be folded. By default the unfolded SFS is used." << std::endl;
	std::cout << std::setw(20) << "-printGrid"           << "  --  " << "Whether Grid should be printed. By default it is not printed." << std::endl;
	std::cout << std::setw(20) << "-SFS STR"                << "  --  " << "MANDATORY: Specify the FILE NAME of SFS input data."  << std::endl;
	std::cout << std::setw(20) << "-Phi STR"                << "  --  " << "MANDATORY: Specify the FILE NAME of Phi input data."  << std::endl;

	std::cout << std::endl;
	
	std::cout << "Examples:" << std::endl	<< std::endl;
	std::cout << "MMC-MaxLikelihood-Grid.out -SFS /path/To/SFSData -Phi /path/To/PhiData " << std::endl;
	std::cout << "MMC-MaxLikelihood-Grid.out -SFS /path/To/SFSData -Phi /path/To/PhiData -folded -lumped 10" << std::endl;

	std::cout << std::endl;
	
}


//******************************************************
//**  INITIALIZING FUNCTIONS  **************************
//******************************************************


void checkArgvInput()
{
	if(SFSInputFile == "false")
	{
		std::cout << "Error: No SFS data provided. See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	if(PhiInputFile == "false")
	{
		std::cout << "Error: No Phi data provided. See ReadMe file for input arguments.\n";
		exit(1);
	}

	if(minAlpha < 1)
	{
		std::cout << "Error: minAlpha = " << minAlpha << ", but needs to be [1,2]. See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	if(maxAlpha > 2)
	{
		std::cout << "Error: maxAlpha = " << maxAlpha << ", but needs to be [1,2]. See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	if(minAlpha > maxAlpha)
	{
		std::cout << "Error: minAlpha > maxAlpha. See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	if(minRho > maxRho)
	{
		std::cout << "Error: minRho > maxRho. See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	if(minMisIdent < 0)
	{
		std::cout << "Error: minMisIdent = " << minMisIdent << ", but needs to be [0, 1]. See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	if(maxMisIdent > 1)
	{
		std::cout << "Error: maxMisIdent = " << maxMisIdent << ", but needs to be [0, 1]. See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	if(minMisIdent > maxMisIdent)
	{
		std::cout << "Error: minMisIdent > maxMisIdent. See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	if(lumping < 0)
	{
		std::cout << "Error: lumping = " << lumping << ", but needs to be non-negative. See ReadMe file for input arguments.\n";
		exit(1);
	}
	

	if(noIncongruentSites == -1)
	{
		if (minMisIdent == 0)
		{
			std::cout << "Note: Misidentification will not be estimated and not be part of the likelihood calculation.\n";
			maxMisIdent = minMisIdent;
			noStepsMisIdent = 0;
			lhoodMisIdent = false;
		}
		else
		{
			std::cout << "Note: Misidentification will not be estimated but be considered fixed at " << minMisIdent << ".\n";
			maxMisIdent = minMisIdent;
			noStepsMisIdent = 0;
		}
	}
	
}

//***checkdataSFS************************************
//  Pre-scans the data to infer the number of
//	data sets and the number of samples
//
void checkdataSFS()
{
	
	char delimLine = '\n';			// line delimiting symbol
	std::string delimiter = "\t";	// this is the delimiter symbol which separates the different columns
	
	std::ifstream csvread;
	const char * c = SFSInputFile.c_str();
	csvread.open(c, std::ifstream::in);
	if(csvread.is_open())
	{
		std::string s;
		int firstRowOnly = 1;
		while( getline(csvread, s, delimLine) ) // every line (i.e., ended by a newline character is read and saved in 's') ...
		{
			s+= "#";							// a '#' is added at the end of the line
			
			size_t pos = 0;
			noDataSets++;						// counts the number of rows (that equals the number of data sets)
			
			if (firstRowOnly)												// to count the number of columns it is enough to do that once
			{
				firstRowOnly = 0;
				while ((pos = s.find('#')) != std::string::npos)					// until end of line has been reached (marked by a '#')...
				{
					std::string val = s.substr(0, pos);
					size_t subpos = 0;
					std::string token;
					while ((subpos = val.find(delimiter)) != std::string::npos)	// ... tokenize string into values seperated by 'delimiter' ...
					{
						sampleSizeSFS++;											// ... and count the number of columns
						token = val.substr(0, subpos);
						val.erase(0, subpos + delimiter.length());
					}
					s = s.substr(pos + 1);
				}
			}
		}
		csvread.close();
		
		if (sampleSizeSFS <= 0)
		{
			std::cerr << "Poorly formatted data file. Less than two samples." << std::endl;
			exit(1);
		}
		
		if (noDataSets <= 0)
		{
			std::cerr << "Poorly formatted data file. No data set found." << std::endl;
			exit(1);
		}
		
	}
	else
	{
		std::cerr << "Could not read data file!" << std::endl;
		exit(1);
	}
	
	if(lumping > sampleSizeSFS)
	{
		std::cerr << "Lumped class bigger than SFS entries!" << std::endl;
		exit(1);
	}
	
	SFS = std::vector<std::vector<double> > (noDataSets, std::vector<double>(sampleSizeSFS, 0));
}


//***checkdataPhi************************************
//  Pre-scans the data to infer the number of
//	data sets and the number of samples
//
void checkdataPhi()
{
	
	char delimLine = '\n';			// line delimiting symbol
	std::string delimiter = "\t";	// this is the delimiter symbol which separates the different columns
	
	std::ifstream csvread;
	const char * c = PhiInputFile.c_str();
	csvread.open(c, std::ifstream::in);
	
	if(csvread.is_open())
	{
		std::string s;
		getline(csvread, s, delimLine);
		s+= "#";							// a '#' is added at the end of the line
		
		size_t pos = 0;
		while ((pos = s.find('#')) != std::string::npos)					// until end of line has been reached (marked by a '#')...
		{
			std::string val = s.substr(0, pos);
			size_t subpos = 0;
			std::string token;
			while ((subpos = val.find(delimiter)) != std::string::npos)	// ... tokenize string into values seperated by 'delimiter' ...
			{
				sampleSizePhi++;											// ... and count the number of columns
				token = val.substr(0, subpos);
				val.erase(0, subpos + delimiter.length());
			}
			s = s.substr(pos + 1);
		}
		
		csvread.close();
		sampleSizePhi--;
		
		if (folded == true)
		{
			if ( (sampleSizePhi + 1)/2 != sampleSizeSFS )
			{
				std::cerr << "Number of phi entries does not match number of SFS entries." << std::endl;
				exit(1);
			}
		}
		else
		{
			if (sampleSizePhi != sampleSizeSFS)
			{
				std::cerr << "Number of phi entries does not match number of SFS entries." << std::endl;
				exit(1);
			}
		}
		
	}
	else
	{
		std::cerr << "Could not read data file!" << std::endl;
		exit(1);
	}
	
	T_Total = std::vector<std::vector<double > > (noStepsAlpha + 1, std::vector<double>(noStepsRho + 1, 0.));
	Phi = std::vector<std::vector<std::vector<double > > >(noStepsAlpha + 1, std::vector<std::vector<double > >(noStepsRho + 1, std::vector<double> (sampleSizeSFS, 0)));
}


//***readSFSFile************************************
//  Reads SFS Data file
//
void readSFSFile()
{
	int _tempRow = 0;
	int _tempCol;
	char delimLine = '\n';			// line delimiting symbol
	std::string delimiter = "\t";	// this is the delimiter symbol which separates the different columns
	std::ifstream csvread;
	
	const char * c = SFSInputFile.c_str();
	csvread.open(c, std::ifstream::in);
	if(csvread.is_open())
	{
		std::string s;
		while( getline(csvread, s, delimLine) ) // every line (i.e., ended by a newline character is read and saved in 's')
		{
			s+= "#";
			size_t pos = 0;
			
			while ((pos = s.find('#')) != std::string::npos) // until end of line has been reached (marked by a '#')...
			{
				std::string val = s.substr(0, pos);
				size_t subpos = 0;
				_tempCol = 0;
				std::string token;
				
				while ((subpos = val.find(delimiter)) != std::string::npos) // ... tokenize string into values seperated by 'delimiter'
				{
					token = val.substr(0, subpos);
					SFS[_tempRow][_tempCol] = atoi(token.c_str());
					_tempCol++;	
					val.erase(0, subpos + delimiter.length());
				}
				SFS[_tempRow][_tempCol] = atoi(val.c_str());
				s = s.substr(pos + 1);
			}
			_tempRow++;
		}
		csvread.close();
	}
	else
	{
		std::cerr << "Could not read data file!" << std::endl;
		exit(1);
	}
}


//***readPhiFile************************************
//  Reads Phi Data file
//
void readPhiFile()
{
	int _tempRho = 1;
	int _tempAlpha = 0;
	char delimLine = '\n';			// line delimiting symbol
	std::string delimiter = "\t";	// this is the delimiter symbol which separates the different columns
	std::ifstream csvread;
	std::vector<std::vector<std::vector<double> > > _initPhi;
	_initPhi = std::vector<std::vector<std::vector<double > > >(noStepsAlpha + 1, std::vector<std::vector<double > >(noStepsRho + 1, std::vector<double> (sampleSizePhi, 0)));
	
	int _index;
	const char * c = PhiInputFile.c_str();
	csvread.open(c, std::ifstream::in);
	if(csvread.is_open())
	{
		std::string s;
		while( getline(csvread, s, delimLine) ) // every line (i.e., ended by a newline character is read and saved in 's')
		{

			if (_tempRho % (noStepsRho + 2) == 0)
			{
				_tempRho = 1;
				_tempAlpha++;
			}
			s+= "#";
			size_t pos = 0;
			
			while ((pos = s.find('#')) != std::string::npos) // until end of line has been reached (marked by a '#')...
			{
				std::string val = s.substr(0, pos);
				size_t subpos = 0;
				_index = 0;
				std::string token;
				
				while ((subpos = val.find(delimiter)) != std::string::npos) // ... tokenize string into values seperated by 'delimiter'
				{
					token = val.substr(0, subpos);
					
					
					if(_index == 0)
					{
						//std::cout << "\t" << _tempAlpha << "\t" << _tempRho << "\t" << atof(token.c_str()) << std::endl;
						T_Total[_tempAlpha][_tempRho - 1] = atof(token.c_str());
					}
					else
					{
						_initPhi[_tempAlpha][_tempRho - 1][_index - 1] = atof(token.c_str());
					}
					_index++;
					
					val.erase(0, subpos + delimiter.length());
				}
				_initPhi[_tempAlpha][_tempRho - 1][_index - 1] = atof(val.c_str());
				s = s.substr(pos + 1);
			}
			_tempRho++;
		}
		csvread.close();
	}
	else
	{
		std::cerr << "Could not read grid file!" << std::endl;
		exit(1);
	}
	
	if(folded == true)
	{
		for (int i = 0; i < noStepsAlpha + 1; ++i)
		{
			for (int j = 0; j < noStepsRho + 1; ++j)
			{
				_index = 0;
				for (int k = 0; k < sampleSizePhi; ++k)
				{
					if(k >= sampleSizeSFS)
					{
						_index--;
					}
					else
					{
						_index++;
					}
					Phi[i][j][_index - 1] += _initPhi[i][j][k];
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < noStepsAlpha + 1; ++i)
		{
			for (int j = 0; j < noStepsRho + 1; ++j)
			{
					Phi[i][j] = _initPhi[i][j];
			}
		}
	}
}


			//////////////////////////////////////////////////
		   /////////            OUTPUT            ///////////
		  //////////////////////////////////////////////////


//**initiateoutput*******************************************
//	Initates all output files
//
void initiateoutput()
{

	char _MLEstimate_file[255];
	strcpy(_MLEstimate_file, outfile_name);
	strcat (_MLEstimate_file, "GridML_");
	strcat (_MLEstimate_file, const_cast<char*> ( dateconst.c_str() ));
	strcat (_MLEstimate_file, "_MLEstimates.txt");
	_offMLEstimate.open(_MLEstimate_file);
	printparameters(& _offMLEstimate);
	_offMLEstimate << "#dataSet" << "\tAlpha" << "\tRho" << "\tMisIdent" << "\tLogL" << "\tThetaW";
	
	_offMLEstimate << std::endl;

	char _L2AbsEstimate_file[255];
	strcpy(_L2AbsEstimate_file, outfile_name);
	strcat (_L2AbsEstimate_file, "GridML_");
	strcat (_L2AbsEstimate_file, const_cast<char*> ( dateconst.c_str() ));
	strcat (_L2AbsEstimate_file, "_L2AbsEstimates.txt");
	_offL2AbsEstimate.open(_L2AbsEstimate_file);
	printparameters(& _offL2AbsEstimate);
	_offL2AbsEstimate << "#dataSet" << "\tAlpha" << "\tRho" << "\tMisIdent" << "\tMinL2Abs" << "\tThetaW";
	
	_offL2AbsEstimate << std::endl;
	
	char _L1AbsEstimate_file[255];
	strcpy(_L1AbsEstimate_file, outfile_name);
	strcat (_L1AbsEstimate_file, "GridML_");
	strcat (_L1AbsEstimate_file, const_cast<char*> ( dateconst.c_str() ));
	strcat (_L1AbsEstimate_file, "_L1AbsEstimates.txt");
	_offL1AbsEstimate.open(_L1AbsEstimate_file);
	printparameters(& _offL1AbsEstimate);
	_offL1AbsEstimate << "#dataSet" << "\tAlpha" << "\tRho" << "\tMisIdent" << "\tMinL1Abs" << "\tThetaW";
	
	_offL1AbsEstimate << std::endl;

	if (print_Grid == true)
	{
		char _Grid_file[255];
		strcpy(_Grid_file, outfile_name);
		strcat (_Grid_file, "GridML_");
		strcat (_Grid_file, const_cast<char*> ( dateconst.c_str() ));
		strcat (_Grid_file, "_Grid.txt");
		_offGrid.open(_Grid_file);
		printparameters(& _offGrid);
		
		_offGrid << "#dataSet" << "\tAlpha" << "\tRho" << "\tMisIdent" << "\tLogL" << "\tL2_Abs" << "\tL1_Abs";
		_offGrid << std::endl;
	}
}


//***printparameters************************************
//  Prints parameters to output file
//
void printparameters(std::ofstream * outfile)
{
	
	*outfile << "# Joined estimation of psi, rho and misindent.";
	*outfile << "\n# VERSION: " << VERSION;
	*outfile << "\n# SFS Input: " << SFSInputFile;
	*outfile << "\n# Phi Input: " << PhiInputFile;
	*outfile << "\n# sampleSize: " << sampleSizeSFS;
	*outfile << "\n# noDataSets: " << noDataSets;
	
	*outfile << "\n# minAlpha: " << minAlpha;
	*outfile << "\n# maxAlpha: " << maxAlpha;
	*outfile << "\n# noStepsAlpha: " << noStepsAlpha;

	*outfile << "\n# minRho: " << minRho;
	*outfile << "\n# maxRho: " << maxRho;
	*outfile << "\n# noStepsRho: " << noStepsRho;
	
	*outfile << "\n# noIncongruentSites: " << noIncongruentSites;
	*outfile << "\n# minMisIdent: " << minMisIdent;
	*outfile << "\n# maxMisIdent: " << maxMisIdent;
	*outfile << "\n# noStepsMisIdent: " << noStepsMisIdent;
	
	*outfile << "\n# lumping: " << lumping;
	*outfile << "\n# folded: " << folded;

	*outfile << std::endl << std::endl;

}

void printMLEstimate(myEstimationMethod myEstimate, int dataSet)
{

	_offMLEstimate << dataSet + 1 << "\t" << myEstimate.getAlphaML() << "\t" << myEstimate.getRhoML() << "\t" << myEstimate.getMisIdentML() << "\t";
	_offMLEstimate.precision(extendedPrecision);
	_offMLEstimate << myEstimate.getMaxLogLikelihood();
	_offMLEstimate.precision(basePrecision);
	_offMLEstimate << "\t" << myEstimate.getThetaML();
	_offMLEstimate << std::endl;
}

void printL2AbsEstimate(myEstimationMethod myEstimate, int dataSet)
{

	_offL2AbsEstimate << dataSet + 1 << "\t" << myEstimate.getAlphal2_abs() << "\t" << myEstimate.getRhol2_abs() << "\t" << myEstimate.getMisIdentl2_abs() << "\t";
	_offL2AbsEstimate.precision(extendedPrecision);
	_offL2AbsEstimate << myEstimate.getMinl2_abs();
	_offL2AbsEstimate.precision(basePrecision);
	_offL2AbsEstimate << "\t" << myEstimate.getThetal2_abs();
	_offL2AbsEstimate << std::endl;
}


void printL1AbsEstimate(myEstimationMethod myEstimate, int dataSet)
{

	_offL1AbsEstimate << dataSet + 1 << "\t" << myEstimate.getAlphal1_abs() << "\t" << myEstimate.getRhol1_abs() << "\t" << myEstimate.getMisIdentl1_abs() << "\t";
	_offL1AbsEstimate.precision(extendedPrecision);
	_offL1AbsEstimate << myEstimate.getMinl1_abs();
	_offL1AbsEstimate.precision(basePrecision);
	_offL1AbsEstimate << "\t" << myEstimate.getThetal1_abs();
	_offL1AbsEstimate << std::endl;
}


void printGrid(myEstimationMethod myEstimate, int dataSet)
{

	for(int j = 0; j < (myEstimate.getAlphaRange().size()); j++)
	{
		for(int i = 0; i < (myEstimate.getRhoRange().size()); i++)
		{
			for(int h = 0; h < (myEstimate.getMisIdentRange().size()); h++)
			{
				_offGrid << dataSet + 1 << "\t" << (myEstimate.getAlphaRange())[j] << "\t" << (myEstimate.getRhoRange())[i] << "\t" << (myEstimate.getMisIdentRange())[h] << "\t";
				_offGrid.precision(extendedPrecision);
				_offGrid << (myEstimate.getLLGrid())[j][i][h] << "\t" << (myEstimate.getml2Grid_abs())[j][i][h] << "\t" << (myEstimate.getml1Grid_abs())[j][i][h] << std::endl;
				_offGrid.precision(basePrecision);
			}
		}
	}
}

//***get_date********************************************
// Gets current date and time for naming the outputfile uniquely
// Format is dd_mm_yyyy_hh_minmin_secsec
//
std::string get_date()
{
	time_t now;
	char the_date[MAX_DATE];
	
	now = time(NULL);
	
	if (now != -1)
	{
		strftime(the_date, MAX_DATE, "%d_%m_%Y_%H%M%S", gmtime(&now));
	}
	
	return (std::string (the_date));
}

//***convert_NumberToString*****************************
// Converts number to string
// Needed for setting output file name
//
template <typename _typeName>
std::string convert_NumberToString ( _typeName Number )
{
	std::string _temp;
	std::ostringstream _convert;
	_convert << Number;
	_temp = _convert.str();
	
	return (_temp);
}

//***convert_StringToChar*****************************
// Converts string to char
// Needed for setting output file name
//
char * convert_StringToChar (std::string _string)
{
	return (const_cast<char*> (_string.c_str()));
}

char * set_filename (char outfile_Name[])
{
	unsigned long int found = SFSInputFile.find_last_of("/\\");
	unsigned long int found2 = SFSInputFile.find_last_of(".");
	
	outfile_prefix = SFSInputFile.substr(0, found);
	std::string datafileName = SFSInputFile.substr(found+1,found2-found-1);
	
	strcpy (outfile_Name, convert_StringToChar(outfile_prefix));
	strcat (outfile_Name, NAVIGATE);
	strcat (outfile_Name, datafileName.c_str());
	strcat (outfile_Name, "-Lump-");
	strcat (outfile_Name, convert_StringToChar (convert_NumberToString(lumping)));
	strcat (outfile_Name, "-");
	
	return (outfile_Name);
}

