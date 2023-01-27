/*
 * MMC-CoalescentSimulator is used to simulate gene trees under the beta- and
 * kingman coalescent process with exponential growth.
 *
 * Copyright (C) 2016 Sebastian Matuszewski & Marcel Hildebrandt
 *
 * This file is part of MMC_Growth.
 *
 * MMC_Growth is free software: you can redistribute it and/or modify
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


#include "myTree.hpp"

////////////////////
// help functions //
////////////////////

int  removeRandomElement(std::vector<int>& v);
double Beta(double x, double y);
std::vector<long double> computeCoalescenceRates(double alpha, int i);
long double logbinomial(long double N, long double k);
long double n_choose_k(int N, int k);
void ReadNewick(int& current, std::string Newick, std::vector<myNode> &AllNodes, int & SampleSize);

////////////////////
//   constructor  //
////////////////////



///////////////////////
//   Beta-Coalescent  //
///////////////////////


myTree::myTree(int SampleSize, double alpha, double mu, long double rho, unsigned int N)
:mBranchlengthsCumulative(SampleSize-1,0.0), mSFS(SampleSize-1)
{
	
	double time = 0;                                //initial time
	mSampleSize = SampleSize;
	mTimeScale = 1.;
	std::vector<myNode> AllNodes (SampleSize);
	std::vector<int> ActiveNodes (SampleSize);      //keeps tracks of the indices of the active nodes (indices correspond to location in AllNodes)
	for (int i = 0; i < SampleSize; ++i)
	{
		AllNodes[i] = myNode(0, i,1);               // create leaves with time 0, label i, 1 descendent
		AllNodes[i].addDescendent(i);               // add itsself as Descendent
		ActiveNodes[i] = i;                         //activate leaves
	}

	double theta;
	if (alpha == 2)
	{
	 theta = 2*mu*N;
	}
	else
	{
		theta = 2*mu*std::pow(N, alpha);
	}

	
	int i = SampleSize;                             //number of ac active nodes
	int lable = i ;                                 //label for the next node
	while (i > 1)
	{
		
		////////////////////////
		// coalescent events  //
		////////////////////////
		
		std::vector< long double> lambda = computeCoalescenceRates(alpha, i);
		double dt = randexp(lambda[0]);              //sample the waiting time until the next coalescensce event
		lambda.erase(lambda.begin());
		int x = randdiscrete(lambda) + 2; time += dt;                                 //coalescensce time
		AllNodes.push_back(myNode(time, lable,0));
		mMerger.push_back(std::vector<double>(4));
		mMerger.back()[0] = i;
		mMerger.back()[1] = x;
		mMerger.back()[2] = time;
		mMerger.back()[3] = N;
		
		for (int j = 0; j < x; ++j)
		{
			int random_index = removeRandomElement(ActiveNodes); //pick one active node at random and deactivate/ remove from active set
			AllNodes.back().addChild(AllNodes[random_index]);
		}
		ActiveNodes.push_back(lable); //add new node to active set
		++lable;
		i = int(ActiveNodes.size()); //update size of active nodes
	}
	mTotalBranchlength = 0.;
	mNoSegSites = 0;
	mSequences.resize(mSampleSize);
	for (int j = 0; j < AllNodes.size(); ++j)
	{
		AllNodes[j].sortDescendent();
		if (j < mSampleSize)
		{
			mSequences[j].push_back(0);
		}
		
	}
	mTreeWalk(int(AllNodes.size()-1), AllNodes, rho, theta, alpha);
	mNewick += ";";
	
	if (rho > 0)
	{
		mTMRC = std::log(1. + time*rho)/(rho); //transformation/ time scaling according to coalescent model (\tilde\rho)
	}
	else
	{
		mTMRC = time;
	}
}


//////////////
//  Newick  //
//////////////

myTree::myTree(std::string Newick, double alpha, double mu, unsigned int N, long double rho)
{
	mTimeScale = 1.;
	std::vector<myNode> AllNodes ;
	int current = 0;
	int SampleSize = 0;
	ReadNewick(current, Newick, AllNodes, SampleSize);
	mSampleSize = SampleSize;
	mBranchlengthsCumulative.resize(SampleSize-1);
	mSFS.resize(SampleSize-1);
	mSequences.resize(mSampleSize);
	mTotalBranchlength = 0.;
	mNoSegSites = 0;
	for (int j = 0; j < AllNodes.size(); ++j)
	{
		AllNodes[j].sortDescendent();
		if (j < mSampleSize)
		{
			mSequences[j].push_back(0);
		}
	}
	
	double theta;
	if (alpha == 2)
	{
		theta = 2*mu*N;
	}
	else
	{
		theta = 2*mu*std::pow(N, alpha);
	}

	mTreeWalk(0, AllNodes, rho, theta, alpha);
	mNewick += ";";
	if (rho > 0)
	{
			mTMRC = std::log(AllNodes[0].getTime()*rho+1.)/(rho); //transformation/ time scaling
	}
	else
	{
		mTMRC = AllNodes[0].getTime();
	}
    
}


////////////////////
// help functions //
////////////////////

void ReadNewick(int& current, std::string Newick, std::vector<myNode> &AllNodes, int& sampleSize)
{
	if (Newick[current] == 40 ) //opening bracket
	{
		AllNodes.push_back(myNode());
		AllNodes.back().setLable((int) AllNodes.size()-1);
		int i = (int) AllNodes.size()-1;
		++current;
		while (Newick[current] != 41)                              //check if closing bracket
		{
			if  ((Newick[current] >= 48 &&  Newick[current] <= 57) || (Newick[current] >= 65 &&  Newick[current] <= 90) || (Newick[current] >= 97 &&  Newick[current] <= 122)) //check if digit
			{
				sampleSize++;
				AllNodes.push_back(myNode(0, (int) AllNodes.size(),1));
				AllNodes.back().addDescendent((int) AllNodes.size()-1);
				AllNodes[i].addChild(AllNodes.back());
				while(Newick[current] != 58) //as long as next element is not a colon iterate
					++current;
				++current;
				std::string time_string;
				while(Newick[current] != 41 && Newick[current] != 44) //as long as next element is not a colon iterate
				{
					time_string += Newick[current];
					++current;
				}
				--current;
				double time = atof(time_string.c_str());
				AllNodes[i].setTime(time);
			}
			if (Newick[current] == 40) //opening bracket
			{
				int next_index = (int) AllNodes.size();
				ReadNewick(current, Newick, AllNodes, sampleSize);
				AllNodes[i].addChild(AllNodes[next_index]);
				while(Newick[current] != 58)    //as long as next element is not a colon iterate
					++current;
				++current;
				std::string time_string;
				while(Newick[current] != 41 && Newick[current] != 44) //as long as next element is not a colon iterate
				{
					time_string += Newick[current];
					++current;
				}
				--current;
				double time = atof(time_string.c_str())+AllNodes[next_index].getTime();
				AllNodes[i].setTime(time);
			}
			++current;
		}
	}
}


int removeRandomElement(std::vector<int>& v)
{
	int n = int(v.size());
	int i = randint(n);
	int x = v[i];
	v.erase(v.begin() +i);
	return x;
}

double Beta(double x, double y)
{
	return std::exp(std::log(std::tgamma(x))+std::log(std::tgamma(y))-std::log(std::tgamma(x+y)));
}

std::vector<long double> computeCoalescenceRates(double alpha, int i)
{
	std::vector<long double> lambda (i);
	if (alpha == 2)
	{
		lambda[1] = n_choose_k(i, 2);
		lambda[0] = lambda[1];
		for (int x = 3; x <= i; ++x)
		{
			lambda[x-1] = 0.;
		}
	}
	else
	{
		lambda[1] = std::exp( std::log(n_choose_k(i, 2)) + std::log(Beta(2.0 - alpha, i - 2.0 + alpha)) - std::log(Beta(2.0 - alpha, alpha)) );
		lambda[0] = lambda[1];
		for (int x = 3; x <= i; ++x)
		{
			lambda[x-1] = std::exp( std::log(n_choose_k(i, x)) + std::log(Beta(x - alpha, i - x + alpha)) - std::log(Beta(2.0 - alpha, alpha)) );
			lambda[0] += lambda[x-1];
		}
	}
	return lambda;
}


void myTree::mTreeWalk(int node_lable, std::vector<myNode>& AllNodes, long double rho, double theta, double alpha)
{
	if (AllNodes[node_lable].getChildren().size() == 0) // node is a leave
	{
		mNewick += "A_";
		mNewick+= convert_NumberToStringTree(AllNodes[node_lable].getLable());
	}
	else
	{   mNewick+= "(";
		int i = 0;
		while (i < AllNodes[node_lable].getChildren().size())
		{
			int IndexChild = AllNodes[node_lable].getChildren()[i]; // location of the i-th child in AllNodes
			mTreeWalk(IndexChild, AllNodes, rho, theta, alpha);
			mNewick+= ":";
			double branch_length;
			if (rho > 0)
			{
				branch_length= std::log(1+AllNodes[node_lable].getTime()*rho)/(rho) - std::log(1+AllNodes[(AllNodes[node_lable].getChildren()[i])].getTime()*rho)/(rho);
			}
			else
			{
				branch_length = AllNodes[node_lable].getTime()-AllNodes[(AllNodes[node_lable].getChildren()[i])].getTime();
			}
			
			double MutationRate = branch_length*theta/2.;
			int noMutations = randpoisson(MutationRate);
			mSFS[AllNodes[IndexChild].getNumberOfDescendents()-1] += noMutations;
			mNoSegSites += noMutations;
			if (noMutations > 0)
			{
				std::vector<bool> zeros(noMutations,0);
				std::vector<bool> ones (noMutations,1);
				
				int counterDescendents = 0;
				for (int k = 0; k < mSampleSize; ++k)
				{
					if (AllNodes[IndexChild].getDescendents()[counterDescendents] == k)
					{
						mSequences[k].insert(mSequences[k].end(),ones.begin(),ones.end());
						++counterDescendents;
					}
					else
					{
						mSequences[k].insert(mSequences[k].end(),zeros.begin(),zeros.end());
					}
				}
			}
			branch_length/=mTimeScale;
			mBranchlengths.push_back(branch_length);
			mTotalBranchlength += branch_length;
			mNumberOfDescendents.push_back(AllNodes[IndexChild].getNumberOfDescendents());
			mBranchlengthsCumulative[AllNodes[IndexChild].getNumberOfDescendents()-1] += branch_length;
			mNewick+= convert_NumberToStringTree(branch_length);
			if (i < AllNodes[node_lable].getChildren().size()-1)
			{
				mNewick+= ",";
			}
			++i;
		}
		mNewick+= ")";
	}
}


/////////////
// getter //
///////////


std::vector<double> myTree::getBranchlengths()
{
	return this->mBranchlengths;
}


std::vector<int> myTree::getNumberOfDescendents()
{
	return this->mNumberOfDescendents;
}

std::string myTree::getNewick()
{
	return this->mNewick;
}


int myTree::getSampleSize()
{
	return this->mSampleSize;
}

std::vector<int> myTree::getSFS()
{
	return this->mSFS;
}

std::vector<double> myTree::getBranchlengthsCumulative()
{
	return this->mBranchlengthsCumulative;
}

double myTree::getTMRCA()
{
	return this->mTMRC;
}

double myTree::getTotalBranchlength()
{
	return this->mTotalBranchlength;
}

int myTree::getNoSegSites()
{
	return this->mNoSegSites;
}

std::vector<std::vector<bool> > myTree::getSequence()
{
	return this->mSequences;
}

std::vector<std::vector<double> > myTree::getMerger()
{
	return this->mMerger;
}

long double n_choose_k( int N,  int k)
{
	if(N < k || N < 0 || k < 0)
		return 0;
	else
	{
		return boost::math::binomial_coefficient<long double>(N, k);
	}
}


// Alternative way to calculate the log-binomial (taken from Hybrid-Lambda)
long double logbinomial(long double N, long double k)
{
	assert(N >= 0.0L);
	assert(k >= 0.0L);
	assert(k <= N);
	
	long double x = 0.0L;
	long double bi = 0.0L;
	
	while (k > x)
	{
		bi += logl(N-x);
		bi -= logl(k-x);
		x++;
	}
	return bi;
}


//***convert_NumberToString*****************************
// Converts number to string
// Needed for setting output file name
//
template <typename _typeName>
std::string convert_NumberToStringTree ( _typeName Number )
{
	std::string _temp;
	std::ostringstream _convert;
	_convert << Number;
	_temp = _convert.str();
	
	return (_temp);
}

