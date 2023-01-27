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


#ifndef myEstimationMethod_hpp
#define myEstimationMethod_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>

class myEstimationMethod
{
	double mMaxLogLikelihood;
	double mRhoML;
	double mAlphaML;
	double mMisIdentML;
    double mMinl2_abs;
    double mAlphal2_abs;
    double mRhol2_abs;
	double mMisIdentl2_abs;
	
    double mMinl1_abs;
    double mAlphal1_abs;
	double mRhol1_abs;
	double mMisIdentl1_abs;
	
	double mThetaML;
	double mThetal2_abs;
	double mThetal1_abs;
	
	bool mSuccess;
	unsigned long mNoSteps;
	std::vector<double> mRhoRange;
	std::vector<double> mAlphaRange;
	std::vector<double> mMisIdentRange;
	std::vector<double> mLogLikelihood;
	
	std::vector<std::vector<std::vector<bool> > > mSuccessGrid;
	std::vector<std::vector<std::vector<double> > > mLLGrid;
    std::vector<std::vector<std::vector<double> > > ml1Grid_abs;
    std::vector<std::vector<std::vector<double> > > ml2Grid_abs;
	std::vector<std::vector<std::vector<double> > > mGrid_Alpha;
	std::vector<std::vector<std::vector<double> > > mGrid_Rho;
	std::vector<std::vector<std::vector<double> > > mGrid_MisIdent;
	
public:

	///////////////////
	//  constructor //
	/////////////////
	
	// :: Grid Search ::
	myEstimationMethod(std::vector<std::vector<std::vector<double>>> &Phi,
					   std::vector<std::vector<double>> &T_TotalGrid,
					   std::vector<double> SFS, double MinAlpha,
					   double MaxAlpha,
					   int noStepsAlpha,
					   double MinRho,
					   double MaxRho,
					   int noStepsRho,
					   double minMisIdent,
					   double maxMisIdent,
					   int noStepsMisIdent,
					   int lumping,
					   int noIncongruentSites,
					   bool lhoodMisIdent
					   );
	
	//////////////////////
	//  help functions //
	////////////////////
	
	/////////////
	// getter //
	///////////
	
	unsigned long getNoSteps();
	double getMaxLogLikelihood();
	std::vector<double> getRhoRange();
	std::vector<double> getAlphaRange();
	std::vector<double> getMisIdentRange();
	double getRhoML();
	double getAlphaML();
	double getMisIdentML();
	bool getSuccess();
	std::vector<std::vector<std::vector<bool> > > getGrid_Success();
	std::vector<std::vector<std::vector<double> > > getLLGrid();
	std::vector<std::vector<std::vector<double> > > getGrid_Alpha();
	std::vector<std::vector<std::vector<double> > > getGrid_Rho();
	std::vector<std::vector<std::vector<double> > > getGrid_MisIdent();
	
	double getMinl2_abs();
	double getAlphal2_abs();
	double getMinl1_abs();
	double getAlphal1_abs();
	
	double getRhol2_abs();
	double getRhol1_abs();
	
	double getMisIdentl2_abs();
	double getMisIdentl1_abs();
	
	std::vector<std::vector<std::vector <bool > > > getSuccessGrid();
	std::vector<std::vector<std::vector<double> > > getml2Grid_abs();
	std::vector<std::vector<std::vector<double> > > getml1Grid_abs();
	
	std::vector<double> getLogLikelihood();
	
	double getThetaML();
	double getThetal1_abs();
	double getThetal2_abs();
	
};

double computeLogLikelihood(
							std::vector<double> &_passPhi,
							int sampleSize,
							std::vector<double> &SFS,
							int segSites,
							double misIdent,
							int cutoff,
							int noIncongruentSites,
							bool lhoodMisIdent
							);
double computeDistance(
					   std::vector<double> &expected,
					   int sampleSize,
					   std::vector<double> &observed,
					   int segSites,
					   double misIdent,
					   int cutoff,
					   double scaleParameter
					   );
double computeConstantFactor(
							 std::vector<double> &observed,
							 unsigned int segSites,
							 int cutoff
							 );
double computeLogChoose(int n, int k);
extern bool folded;


#endif /* myEstimationMethod_hpp */

