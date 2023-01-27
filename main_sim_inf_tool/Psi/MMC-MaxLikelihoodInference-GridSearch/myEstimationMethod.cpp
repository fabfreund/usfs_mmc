/*
 * MMC-MaxLikelihood-Grid is used to estimate the MMC parameter psi and/or the
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


#include "myEstimationMethod.hpp"

////////////////////
//   grid search  //
////////////////////

myEstimationMethod::myEstimationMethod(
									   std::vector<std::vector<std::vector<double>>> &Phi,
									   std::vector<std::vector<double>> &T_TotalGrid,
									   std::vector<double> SFS,
									   double MinPsi,
									   double MaxPsi,
									   int noStepsPsi,
									   double MinRho,
									   double MaxRho,
									   int noStepsRho,
									   double minMisIdent,
									   double maxMisIdent,
									   int noStepsMisIdent,
									   int lumping,
									   int noIncongruentSites,
									   bool lhoodMisIdent
									   ):
mLLGrid(noStepsPsi + 1, std::vector<std::vector<double>>(noStepsRho + 1, std::vector<double>(noStepsMisIdent + 1, -std::numeric_limits<double>::infinity()))),
ml1Grid_abs(noStepsPsi + 1, std::vector<std::vector<double>>(noStepsRho + 1, std::vector<double>(noStepsMisIdent + 1, std::numeric_limits<double>::infinity()))),
ml2Grid_abs(noStepsPsi + 1, std::vector<std::vector<double>>(noStepsRho + 1, std::vector<double>(noStepsMisIdent + 1, std::numeric_limits<double>::infinity()))),
mMaxLogLikelihood (- std::numeric_limits<double>::infinity()),
mRhoRange(noStepsRho + 1),
mPsiRange(noStepsPsi + 1),
mMisIdentRange(noStepsMisIdent + 1),
mMinl1_abs(std::numeric_limits<double>::infinity()),
mMinl2_abs(std::numeric_limits<double>::infinity()),
mThetaML(NAN),
mThetal2_abs(NAN),
mThetal1_abs(NAN)
{
	int sampleSize = (int) SFS.size() + 1;
	
	if (lumping != 0)
	{
		while (SFS.size() > lumping)
		{
			SFS[lumping-1] += SFS.back();
			SFS.pop_back();
		}
	}
	
	if (noStepsRho == 0)
	{
		mRhoRange[0] = MinRho;
	}
	else
	{
		for (int i = 0; i < noStepsRho + 1; ++i)
		{

			mRhoRange[i] = MinRho + i*(MaxRho - MinRho)/(noStepsRho);
		}
	}
	
	if (noStepsPsi == 0)
	{
		mPsiRange[0] = MinPsi;
	}
	else
	{
		for (int j = 0; j < noStepsPsi + 1; ++j)
		{
			mPsiRange[j] = MinPsi + j*(MaxPsi - MinPsi)/(noStepsPsi);
		}
	}
	
	if (noStepsMisIdent == 0)
	{
		mMisIdentRange[0] = minMisIdent;
	}
	else
	{
		for (int j = 0; j < noStepsMisIdent + 1; ++j)
		{
			mMisIdentRange[j] = minMisIdent + j*(maxMisIdent - minMisIdent)/(noStepsMisIdent);
		}
	}
	
    int cutoff;
    if (lumping == 0)
    {
        cutoff = sampleSize-1;
    }
    else
    {
        cutoff = lumping;
    }
	
    unsigned int segSites = 0;
    std::vector<double> relBranchLength_observed(cutoff);
	
	for (int j = 0; j < cutoff; ++j)
    {
        segSites += SFS[j];
        relBranchLength_observed[j] = SFS[j];
    }
	
	if (noIncongruentSites == -1)
	{
		noIncongruentSites = int(std::round(2 * minMisIdent * segSites));
	}
	
	if (lumping == 0)
	{
		if (folded == true)
		{
			cutoff = (sampleSize) / 2;
		}
		else
		{
			cutoff = sampleSize - 1;
		}
	}
	else
	{
		if (folded == true)
		{
			cutoff = (sampleSize) / 2;
			cutoff = std::min(cutoff, lumping);
		}
		else
		{
			cutoff = lumping;
		}
	}

    std::transform(relBranchLength_observed.begin(), relBranchLength_observed.end(), relBranchLength_observed.begin(), std::bind1st(std::multiplies<double>(),1.0/((double) segSites)));
	
	std::vector<double> _Phi;
	std::vector<double> _passPhi(cutoff, 0.);
	std::vector<double> _SFS;
	std::vector<double> _passSFS(cutoff, 0.);
	double tTotal, theta;
	int index;
	
	double constantLogTerm = computeConstantFactor(SFS, segSites, cutoff);
	
	for (int k = 0; k < noStepsPsi + 1; ++k)
	{
		for (int l = 0; l < noStepsRho + 1; ++l)
		{

            tTotal = T_TotalGrid[k][l];
            theta = segSites / tTotal;
            _Phi = Phi[k][l];
            _SFS = _Phi;
			std::fill(_passPhi.begin(), _passPhi.end(), 0.);
			std::fill(_passSFS.begin(), _passSFS.end(), 0.);
            std::transform(_SFS.begin(), _SFS.end(), _SFS.begin(), std::bind1st(std::multiplies<double>(), segSites));
			
			for (int j = 0; j < sampleSize-1; ++j)
			{
				index = std::min(cutoff-1, j);
				_passPhi[index] += _Phi[j];
				_passSFS[index] += _SFS[j];
			}
			
			for (int m = 0; m < noStepsMisIdent + 1; m++)
			{
			
				mLLGrid[k][l][m] = constantLogTerm + computeLogLikelihood(_passPhi, sampleSize, SFS, segSites, mMisIdentRange[m], cutoff, noIncongruentSites, lhoodMisIdent);
				ml1Grid_abs[k][l][m] = computeDistance(_passSFS, sampleSize, SFS, segSites, mMisIdentRange[m], cutoff, 1.0);
				ml2Grid_abs[k][l][m] = computeDistance(_passSFS, sampleSize, SFS, segSites, mMisIdentRange[m], cutoff, 2.0);
				
				if (mMaxLogLikelihood < mLLGrid[k][l][m])
				{
					mRhoML = mRhoRange[l];
					mPsiML = mPsiRange[k];
					mMisIdentML = mMisIdentRange[m];
					mMaxLogLikelihood = mLLGrid[k][l][m];
					mNoSteps = 1;
					mThetaML = theta;
				}
				
				if (mMinl1_abs > ml1Grid_abs[k][l][m])
				{
					mRhol1_abs = mRhoRange[l];
					mPsil1_abs = mPsiRange[k];
					mMisIdentl1_abs = mMisIdentRange[m];
					mMinl1_abs = ml1Grid_abs[k][l][m];
					mNoSteps = 1;
					mThetal1_abs = theta;
				}
				if (mMinl2_abs > ml2Grid_abs[k][l][m])
				{
					mRhol2_abs = mRhoRange[l];
					mPsil2_abs = mPsiRange[k];
					mMisIdentl2_abs = mMisIdentRange[m];
					mMinl2_abs = ml2Grid_abs[k][l][m];
					mNoSteps = 1;
					mThetal2_abs = theta;
				}
			}
			
		}
	}
}


/////////////////////////
//    help functions   //
/////////////////////////

double computeLogLikelihood(
							std::vector<double> &_passPhi,
							int sampleSize,
							std::vector<double> &SFS,
							int segSites,
							double misIdent,
							int cutoff,
							int noIncongruentSites,
							bool lhoodMisIdent
							)
{
	//Will not cope with cutoff properly
	double _logLikelihood = 0;
	double fracIncongruentSites;
	for (int k = 0; k < cutoff; ++k)
	{
		_logLikelihood += SFS[k] * std::log(_passPhi[k] * (1 - misIdent) + _passPhi[cutoff - 1 - k] * misIdent);
	}
	
	if (lhoodMisIdent == true && misIdent > 0)
	{
		fracIncongruentSites = 2 * misIdent / (1 + 2 * misIdent);
		_logLikelihood += computeLogChoose(segSites + noIncongruentSites, noIncongruentSites) +  noIncongruentSites * std::log(fracIncongruentSites) + segSites * std::log(1 - fracIncongruentSites);
	}
	
	return (_logLikelihood);
}


double computeDistance(std::vector<double> &expected, int sampleSize, std::vector<double> &observed, int segSites, double misIdent, int cutoff, double scaleFactor)
{
	//Will not cope with cutoff properly
	double distance = 0;
    for (int k = 0; k < cutoff; ++k)
    {
		distance += std::pow((std::abs(observed[k] - (expected[k] * (1 - misIdent) + expected[cutoff - 1 - k] * misIdent))), scaleFactor);
    }
    return (std::pow(distance, 1.0 / scaleFactor));
}

double computeConstantFactor(std::vector<double> &observed, unsigned int segSites, int cutoff)
{
	double _constF = segSites * std::log(segSites) - segSites;
	double _temp;
	for(int i = 0; i < cutoff; i++)
	{
		_temp = 0;
		
		for(int j = 1; j <= observed[i]; j++)
		{
			_temp += std::log(j);
		}
		_constF -= _temp;
	}
	return (_constF);
}

double computeLogChoose(int n, int k)
{
	double _result = 0;
	
	for(int i = 1; i <=k; i++)
	{
		_result += std::log(n+1-i) - std::log(i);
	}
	
	return (_result);
}

/////////////
// getter //
///////////

double myEstimationMethod::getMaxLogLikelihood()
{
	return (this->mMaxLogLikelihood);
}

std::vector<double> myEstimationMethod::getRhoRange()
{
	return (this->mRhoRange);
}

std::vector<double> myEstimationMethod::getPsiRange()
{
	return (this->mPsiRange);
}

std::vector<double> myEstimationMethod::getMisIdentRange()
{
	return (this->mMisIdentRange);
}

std::vector<std::vector<std::vector<double> > > myEstimationMethod::getLLGrid()
{
	return (this->mLLGrid);
}

std::vector<std::vector<std::vector<double> > > myEstimationMethod::getGrid_Psi()
{
	return (this->mGrid_Psi);
}

std::vector<std::vector<std::vector<bool> > > myEstimationMethod::getGrid_Success()
{
	return (this->mSuccessGrid);
}

std::vector<std::vector<std::vector<double> > > myEstimationMethod::getGrid_Rho()
{
	return (this->mGrid_Rho);
}

std::vector<std::vector<std::vector<double> > > myEstimationMethod::getGrid_MisIdent()
{
	return (this->mGrid_MisIdent);
}

bool myEstimationMethod::getSuccess()
{
	return (this->mSuccess);
}

double myEstimationMethod::getRhoML()
{
	return (this->mRhoML);
}

double myEstimationMethod::getPsiML()
{
	return (this->mPsiML);
}

double myEstimationMethod::getMisIdentML()
{
	return (this->mMisIdentML);
}

unsigned long myEstimationMethod::getNoSteps()
{
	return (this->mNoSteps);
}

double myEstimationMethod::getMinl2_abs()
{
	return (this->mMinl2_abs);
}

double myEstimationMethod::getMinl1_abs()
{
	return (this->mMinl1_abs);
}

double myEstimationMethod::getPsil2_abs()
{
	return (this->mPsil2_abs);
}

double myEstimationMethod::getPsil1_abs()
{
	return (this->mPsil1_abs);
}

double myEstimationMethod::getRhol2_abs()
{
	return (this->mRhol2_abs);
}

double myEstimationMethod::getRhol1_abs()
{
	return (this->mRhol1_abs);
}

double myEstimationMethod::getMisIdentl2_abs()
{
	return (this->mMisIdentl2_abs);
}

double myEstimationMethod::getMisIdentl1_abs()
{
	return (this->mMisIdentl1_abs);
}

std::vector<std::vector<std::vector<bool > > > myEstimationMethod::getSuccessGrid()
{
	return (this->mSuccessGrid);
}

std::vector<std::vector<std::vector<double> > > myEstimationMethod::getml2Grid_abs()
{
	return (this->ml2Grid_abs);
}

std::vector<std::vector<std::vector<double> > > myEstimationMethod::getml1Grid_abs()
{
	return (this->ml1Grid_abs);
}

std::vector<double> myEstimationMethod::getLogLikelihood()
{
	return (this->mLogLikelihood);
}

double myEstimationMethod::getThetaML()
{
	return (this->mThetaML);
}

double myEstimationMethod::getThetal1_abs()
{
	return (this->mThetal1_abs);
}
double myEstimationMethod::getThetal2_abs()
{
	return (this->mThetal2_abs);
}
